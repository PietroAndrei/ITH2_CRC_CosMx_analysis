library(Seurat)
library(sf)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
utils = snakemake@params[["utils"]]
source(utils)
meta = readRDS(snakemake@input[["meta"]])
#Clean cell type annotation after InSituType refinement step
meta[str_detect(meta$final_anno,'Monocytes'),]$final_anno = 'Monocytes_Macrophages'
meta[str_detect(meta$final_anno,'Epithelial'),]$final_anno = 'Epithelial'
meta[str_detect(meta$final_anno,'Endothelial'),]$final_anno = 'Endothelial'
meta[str_detect(meta$final_anno,'Fibroblasts'),]$final_anno = 'Fibroblasts'
meta[str_detect(meta$final_anno,'TCD8'),]$final_anno = 'TCD8'
#Compute distance between each Mono/Macro cell and their closest CD8 T cell
mac = get_cosmx_ct_dists(meta,ct.1='Monocytes_Macrophages',ct.2='TCD8',class='final_anno',plot=F)
#Try different thresholds to define TAM-T and TAM-Other groups
mact = mac %>% filter(TCD8_dist <= 20) 
maco = mac %>% filter(TCD8_dist > 20) 
mac$group_20um = ifelse(mac$cell_ID %in% rownames(mact),'TAM-T','TAM-Other')
#Compute distance between each Epithelial cell and their closest CD8 T cell
epi = get_cosmx_ct_dists(meta,ct.1='Epithelial',ct.2='TCD8',class='final_anno',plot=F)
#Try different thresholds to define EpiT and EpiOnly groups
epit = epi %>% filter(TCD8_dist <= 20) 
epio = epi %>% filter(TCD8_dist > 40)
epi$group_20um = ifelse(epi$cell_ID %in% rownames(epit),'EpiT',ifelse(epi$cell_ID %in% rownames(epio),'EpiOnly','Other'))

#Load seurat dataset and normalize protein expression values
seurat = readRDS(snakemake@input[["seurat"]])
seurat@reductions[names(seurat@reductions)] = NULL
Idents(seurat) = 'cell_ID'
seurat = NormalizeData(seurat)
mat = LayerData(seurat,'data')
markers = rownames(mat)
#subset matrix to keep only TAMs and epithelial cells, and add tam/epi grouping scheme columns
mac_mat = Matrix::t(mat[,rownames(mac)])
mac_mat = as.data.frame(as.matrix(mac_mat))
mac_mat$macType = mac$group_20um
epi_mat = Matrix::t(mat[,rownames(epi)])
epi_mat = as.data.frame(as.matrix(epi_mat))
epi_mat$epiType = epi$group_20um

#Calculate Wilcoxon test pvalue for the two TAM-T/TAM-Other comparisons
#20um
mact_pws = sapply(markers,function(x){
mact_pw = wilcox.test(mac_mat[rownames(mact),x],mac_mat[rownames(maco),x])$p.value
mact_pw = format(mact_pw,scientific=T,digits=3)
mact_pw = paste0('p = ',mact_pw)
return(mact_pw)}
)
names(mact_pws) = markers

mac_mat$macType = factor(mac_mat$macType,levels=c('TAM-T','TAM-Other'))

#Calculate Wilcoxon test pvalue for the Epi-T/Epi-Only comparisons
#20um
epit_pws = sapply(markers,function(x){
epit_pw = wilcox.test(epi_mat[rownames(epit),x],epi_mat[rownames(epio),x])$p.value
epit_pw = format(epit_pw,scientific=T,digits=3)
epit_pw = paste0('p = ',epit_pw)
return(epit_pw)}
)
names(epit_pws) = markers
epi_mat$epiType = factor(epi_mat$epiType,levels=c('EpiT','EpiOnly','Other'))

mac_plots = lapply(markers,
		   function(x){
			   ggplot(mac_mat,aes(x=macType,y=.data[[x]],fill=macType))+
			   theme_classic()+
			   geom_boxplot(width=0.5,outlier.shape =NA)+geom_point(size=0.04)+
			   scale_fill_manual(values = c("TAM-T" = "#83C65D", "TAM-Other" = "#3B4E36"),guide='none')+
			   ylab(paste0(x,' Normalized protein expression'))+xlab('')+
			   ggtitle('TAM next to T (CD8 T cell distance <= 20um) vs Other TAM')+
			   scale_y_continuous(
					      expand=c(0.04,0),breaks = c(min(mac_mat[,x]),max(mac_mat[,x])),
					      labels=c(round(min(mac_mat[,x]),1),round(max(mac_mat[,x]),1))
					      )+
		   theme(plot.title = element_text(size=6))+
		   annotate("text", x = 1.5, y = 8.5, label = mact_pws[x])+
		   scale_x_discrete(labels = c(paste0('TAM-T\n(',nrow(mact),')'),paste0('TAM-Other\n(',nrow(maco),')')))
		   }
		   )



pdf(snakemake@output[["tam_plots"]],width=5,height=5)
for(p in mac_plots){
	print(p)
}
dev.off()

epi_plots = lapply(markers,
		   function(x){
			   ggplot(epi_mat %>% filter(epiType != 'Other'),aes(x=epiType,y=.data[[x]],fill=epiType))+
theme_classic()+
geom_boxplot(width=0.5,outlier.shape=NA)+geom_point(size=0.04)+
scale_fill_manual(values = c("EpiT" = "cyan", "EpiOnly" = "blue"),guide='none')+
ylab(paste0(x,' Normalized protein expression'))+xlab('')+
ggtitle('Epi next to T (CD8 T cell distance <= 20um) vs EpiOnly (>40um)')+
scale_y_continuous(expand=c(0.04,0),breaks = c(min(epi_mat[,x]),max(epi_mat[,x])),
		   labels=c(round(min(epi_mat[,x]),1),round(max(epi_mat[,x]),1))
		   )+
theme(plot.title = element_text(size=6))+
annotate("text", x = 1.5, y = 8.5, label = epit_pws[x])+
scale_x_discrete(labels = c(paste0('Epi-T\n(',nrow(epit),')'),paste0('Epi-Only\n(',nrow(epio),')')))
		   }
		   )


pdf(snakemake@output[["epi_plots"]],width=5,height=5)
for(p in epi_plots){
	print(p)
}
dev.off()
