library(Seurat)
library(dplyr)
library(stringr)
library(Matrix)
library(ggplot2)
utils = snakemake@params[["utils"]]
source(utils)
options(Seurat.object.assay.version = "v5")
#Set the directory containing Pelka related files

pelka_path = snakemake@params[["pelka_path"]]
pelka = list.files(pelka_path,full.names=T)
#Set the directory containing Pelka related files
h5 = pelka[str_detect(pelka,'h5')]
#h5 = paste0(pelka_path,h5)
#Load Pelka metadata and clusters
meta = pelka[str_detect(pelka,'metatables')]
meta = read.delim(meta,sep=',')

obs = pelka[str_detect(pelka,'cluster')]
obs = read.delim(obs,sep=',')

#Keep only epithelial and macrophage cells (including those derived from normal adjacent tissue)
all_cells = obs[(obs$clMidwayPr %in% c('Epi','Macro')),]$sampleID
message("Considering ",length(all_cells)," epithelial/macrophage cells from Pelka dataset")
rownames(meta) = meta$cellID
rownames(obs) = obs$sampleID
meta = meta[all_cells,]
obs = obs[all_cells,]
meta$origCluster = obs$clMidwayPr
rm(obs)

#Divide normal and tumour epithelial cells in two separate groups
#meta[meta$origCluster == 'Epi' & meta$SPECIMEN_TYPE == 'T',]$origCluster = 'TumourEpi'
meta = meta[meta$SPECIMEN_TYPE == 'T',]
#Unify patient ID information (remove tumour/normal specification)
meta$PatientID = sapply(meta$PatientTypeID,function(x) strsplit(x,'_')[[1]][[1]])


#Load Pelka dataset
data = Read10X_h5(filename = h5, use.names = TRUE, unique.features = TRUE)
bcs = colnames(data)

#Filter Pelka dataset
common_bcs = intersect(rownames(meta),bcs)
message(length(all_cells),"/",length(common_bcs)," single cells matching with metadata in Pelka dataset")
data = data[,common_bcs]
meta = meta[common_bcs,]
genes = rownames(data)



#Load NCG7 gene list that will be used to filter the count matrix and update the gene symbols
ncg = snakemake@params[['ncg']]
ncg = read.delim(ncg)
#Update gene symbols in the dataset when necessary
new_symbols = sapply(genes,function(x){
	is_new = x %in% unique(ncg$NCG_symbol)
	is_old = !(x %in% unique(ncg$NCG_symbol)) & (x %in% unique(ncg$symbol))
	if(is_new == F & is_old == T){
		x = unique(ncg[which(ncg$symbol == x),]$NCG_symbol)
	} else {
		x = x
	}
	return(x)
}
)

#Filter matrix to keep only the selected NCG genes
data = data[names(new_symbols),]
#Update gene symbols
rownames(data) = new_symbols
message("Keeping ",nrow(data)," genes included in the NCG7 catalogue")
data=data[!duplicated(rownames(data)),]

#Save matrix as Seurat objects
seurat = CreateSeuratObject(counts = data,meta.data = meta,)
rm(data)
seurat = NormalizeData(seurat)

data  = FetchData(seurat,c('CD74',colnames(seurat@meta.data)))
data$origCluster = factor(data$origCluster,levels = c('Epi','Macro'))

rm(seurat)

#Summarise CD74 expression at the patient level, with/without including cells with no CD74 expression
#If some patients have only <10 cells available for epithelial or macrophage compartment, don't summarise CD74 median expression for those cells
patient_count = data %>% group_by(PatientID,origCluster) %>% count(.)
patient_count = patient_count %>% mutate(label = paste0(PatientID,'_',origCluster))

patient_data = data %>% 
group_by(PatientID,MMRStatus,origCluster) %>% 
summarise(avg_CD74 = median(CD74)) %>% 
ungroup() %>%
mutate(label = paste0(PatientID,'_',origCluster)) %>%
filter(label %in% patient_count$label)

#Paired Wilcoxon test
patient_data = patient_data %>% group_by(PatientID) %>% arrange(origCluster,.by_group=T)
pw1 = wilcox.test(data=patient_data,avg_CD74~origCluster)$p.value
pw1 = format(pw1,scientific=T,digits=3)
pw1 = paste0('p = ',pw1)

#Patient level (ALL)
g1 = ggplot(patient_data, aes(x=origCluster,y=avg_CD74,fill=origCluster))+
theme_classic()+
geom_boxplot(outlier.shape = NA,width=0.6,alpha = 0.2)+
geom_point(aes(color= MMRStatus))+
#geom_line(aes(group = PatientID),color='grey',alpha=0.5)+
#stat_summary(fun='median',fill='white',geom='label',aes(label=round(after_stat(y),2)))+
scale_fill_manual(values=c('Epi' = 'blue','Macro'='forestgreen'))+
scale_colour_manual(values = c('MMRd' = '#8C2F89','MMRp' = '#F29FC4'))+
ylim(0,10)+scale_y_continuous(breaks = c(0,round(max(patient_data$avg_CD74),2)))+
geom_text(aes(label=paste0('n = ',after_stat(count))), y=7, stat='count', colour="black", size=4)+
annotate("text", x = 1.5, y = 8.5, label = pw1)+
xlab('')+ylab('Normalized CD74 expression')+labs(fill='Cell Type')+
guides(color='none',fill='none')+
ggtitle('CD74 expression in Epi/Macrophages Pelka (ALL)')

pdf(snakemake@output[["plots"]],width=6,height=5)
g1
dev.off()



