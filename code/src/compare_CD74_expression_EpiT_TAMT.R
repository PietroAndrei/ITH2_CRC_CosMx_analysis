library(Seurat)
library(sf)
library(dplyr)
library(stringr)
library(ggplot2)
utils = snakemake@params[["utils"]]
source(utils)
meta_path = snakemake@params[["meta_path"]]

samples = c('CR21','CR36','CR48','UH20')
seurat = list.files('Seurat',pattern='Cosmx',full.names=T)
meta = list.files(meta_path,pattern='Cosmx',full.names=T)
#Aggregate CD74 expression levels from all samples/cells
all_cd74 = lapply(samples,function(x){
		  seurat = seurat[str_detect(seurat,x)]
		  meta = meta[str_detect(meta,x)]
		  seurat=readRDS(seurat)
		  meta=readRDS(meta)
		  seurat = seurat[,WhichCells(seurat,cells = rownames(meta))]
		  seurat@meta.data = as.data.frame(meta)
		  Idents(seurat) = 'cell_ID'
		  seurat = NormalizeData(seurat)
		  data = FetchData(seurat,c('CD74',colnames(seurat@meta.data)))
		  rm(seurat)
		  data$final_anno = as.character(data$final_anno)
		  data = st_as_sf(data)
		  #Identify Epi-T/TAM-T cell groups
		  nb = st_intersects(data)
		  nb = lapply(seq_len(length(nb)),function(x) nb[[x]][nb[[x]] != x])
		  attr(nb,'class') = c('sgbp','list')
		  nb = as.data.frame(nb)
		  nb$row.type = as.character(data[nb$row.id,]$final_anno)
		  nb$col.type = as.character(data[nb$col.id,]$final_anno)
		  nb$row.cellID = data[nb$row.id,]$cell_ID
		  nb$col.cellID = data[nb$col.id,]$cell_ID
		  tcells_ie = data %>% filter(final_anno == 'T/NK',nb_epi >= 3)
		  nb_epit = nb %>% filter(col.type == 'Epi',row.cellID %in% tcells_ie$cell_ID)
		  epit = data %>% filter(cell_ID %in% unique(nb_epit$col.cellID))
		  tcells_stroma = data %>% filter(final_anno == 'T/NK',nb_epi == 0)
		  data = data %>% mutate(Epi_groups = ifelse(cell_ID %in% epit$cell_ID,'Epi-T',ifelse(epiType == 'epiOnly' & nb_epi >=4,'Epi-Only','')))
		  data = data %>% mutate(Lymph_groups = ifelse(cell_ID %in% tcells_ie$cell_ID,'IE-T',ifelse(cell_ID %in% tcells_stroma$cell_ID,'Ts-T','')))
		  data = data %>% mutate(Mye_groups = ifelse(final_anno == 'Myeloid' & nb_tnk >=1,'TAM-T',ifelse(final_anno=='Myeloid' & nb_tnk == 0,'TAM-Other','')))
		  return(data)
})

all_cd74 = purrr::reduce(all_cd74,rbind)


all_cd74 = all_cd74 %>% filter(Epi_groups == 'Epi-T' | Mye_groups == 'TAM-T')

#Compare CD74 expression in Epi-T vs TAM-T
pw=wilcox.test(CD74~final_anno,data=all_cd74)$p.value
if(pw == 0){
	pw = paste0('Wilcoxon, p < 2.2e-16')
}

all_cd74[all_cd74$final_anno == 'Epi',]$final_anno = 'Epi-T'
all_cd74[all_cd74$final_anno == 'Myeloid',]$final_anno = 'TAM-T'
#Plot the comparison result as a box plot
pdf(snakemake@output[["plot"]],width=5,height = 5)
ggplot(all_cd74,aes(x=final_anno,y=CD74,color=final_anno))+
theme_classic()+
geom_boxplot(outlier.shape = NA,width=0.6)+
geom_point(size=0.25)+
scale_y_continuous(breaks=c(0,round(max(all_cd74$CD74),2)))+
scale_color_manual(values=c('Epi-T'='blue','TAM-T'='forestgreen'))+
annotate('text',x=1.5,y=9.5,label=pw)+
xlab('')+ylab('Normalized CD74 expression')+labs(fill='Cell Type')+
geom_text(aes(label=paste0('n = ',after_stat(count))), y=8.5, stat='count', colour="black", size=4)+
guides(color='none')+
ggtitle('CD74 expression: Epi-T vs TAM-T')
dev.off()
