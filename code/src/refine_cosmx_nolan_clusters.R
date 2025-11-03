library(Seurat)
library(presto)
library(InSituType)
library(sf)
library(dplyr)
library(stringr)
library(ggplot2)
utils = snakemake@params[["utils"]]
source(utils)

seurat = readRDS(snakemake@input[["seurat"]])
seurat@reductions[names(seurat@reductions)] = NULL
Idents(seurat) = 'cell_ID'
seurat = NormalizeData(seurat) 
neg = round(seurat$nCount_negprobes/10)

immunof = seurat@meta.data[,str_detect(colnames(seurat@meta.data),'Mean\\.|Area.um2|AspectRatio')]
cohort <- fastCohorting(immunof,
			gaussian_transform = TRUE) 

#For InSituType, limit the genes to the top 150 most variable ones
mat = Matrix::t(LayerData(seurat,layer='counts'))
#Round protein intensity to counts
mat = round(mat)

markers = readRDS(snakemake@input[["markers"]])
anno = readRDS(snakemake@input[["anno"]])
markers_filt = markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% top_n(15) %>% ungroup() %>% arrange(as.character(cluster)) %>% as.data.frame()

#Define cell type identity for each cluster based on their protein markers
labels = data.frame(cluster = unique(sort(as.character(anno$InSitu_clusters))),celltype=c('immune','TCD8_CD11b','Bcell','Monocytes_Macrophages','unknown','Neutrophils','Monocytes_Macrophages','immune','immune','Bcell','Monocytes_Macrophages','Epithelial','Epithelial','Epithelial','TCD8_CD56','Endothelial','proliferating','Fibroblasts','TCD4','Fibroblasts'),n=as.numeric(table(anno$InSitu_clusters)))
unsup = readRDS(snakemake@input[["unsup"]])
celltypes = labels$celltype
names(celltypes) = labels$cluster

newclusts <- refineClusters(logliks = unsup$logliks,
			    merges = celltypes) 

anno_aggr = newclusts$clust

set.seed(42)
newclusts = refineClusters(logliks = newclusts$logliks,
			   subcluster = list(immune = 5:10,unknown=5:10,proliferating=5:10),
			   counts =mat,
			   cohort = cohort,
			   neg = neg)
rm(mat)

seurat$refined_anno = newclusts$clust
Idents(seurat) = 'refined_anno'

fc.func = function(x=NULL,pscount=1,base=2){return(log(x = rowMeans(x = expm1(x = x)) + pscount, base = base))}
new_markers = FindAllMarkers(seurat,logfc.threshold = log2(1.5),mean.fxn = fc.func)
saveRDS(new_markers,snakemake@output[["new_markers"]])
saveRDS(seurat@meta.data %>% select(cell_ID,refined_anno),snakemake@output[["new_anno"]])

#Annotate the new unknown/immune/prolif subclusters
unk_clusts = paste0('unknown_',seq(1,10,1))
immune_clusts = paste0('immune_',seq(1,10,1))
prolif_clusts = paste0('proliferating_',seq(1,10,1))

unk_cells = c('unknown','unknown','Epithelial','unknown','Endothelial','TCD8','Epithelial','Fibroblasts','Monocytes_macrophages','Fibroblasts')
prolif_cells = c('Epithelial','Epithelial','Epithelial','prolif_immune','Epithelial','prolif_immune','Epithelial','prolif_immune','prolif_immune','Epithelial')
immune_cells = c('immune','immune','immune','immune','CD8','Monocytes_Macrophages','immune','immune','immune','DC')

#unk_cells = c('Endothelial','TCD4','TCD8','unknown','Epithelial','Epithelial','Epithelial','Monocytes_Macrophages','Fibroblasts','Monocytes_Macrophages')
#immune_cells = c('DC','immune','DC','Tcell','Monocytes_Macrophages','Monocytes_Macrophages','immune','immune','immune','Monocytes_Macrophages')
#prolif_cells = c('prolif_immune','Epithelial','prolif_immune','Epithelial','prolif_immune','prolif_immune','prolif_immune','prolif_immune','Epithelial','Epithelial')

names(unk_cells) = unk_clusts
names(immune_cells) = immune_clusts
names(prolif_cells) = prolif_clusts

#Finalize celltype annotation
finalclusts = refineClusters(logliks = newclusts$logliks,
			     merges = c(unk_cells,immune_cells,prolif_cells)
			     )

saveRDS(finalclusts,snakemake@output[["final_unsup"]])
seurat$final_anno = finalclusts$clust

meta_sf = seurat@meta.data
meta_sf$x_slide_um = meta_sf$x_slide_mm*1000
meta_sf$y_slide_um = meta_sf$y_slide_mm*1000
meta_sf = make_geometry_df(meta_sf,xcoord='x_slide_um',ycoord='y_slide_um')
meta_sf$final_anno = as.character(meta_sf$final_anno)
saveRDS(meta_sf,snakemake@output[["meta_sf"]])

all_celltypes = unique(meta_sf$final_anno)
colscheme = c24[1:length(all_celltypes)]
names(colscheme) = all_celltypes

meta_list = split(meta_sf,meta_sf$fov)
plots = snakemake@output[["plots"]]
plot_list = lapply(meta_list,function(x){
		   g = x %>% sf::st_as_sf() %>%
		   ggplot()+geom_sf(aes(fill=final_anno),lwd=0.01,shape=21,size=2)+
		   scale_fill_manual(values=colscheme)+labs(fill='Cell Type')
		   g = g + theme(
				 panel.background = element_rect(fill = "black",
								 colour = "black",
								 linewidth  = 0.5, linetype = "solid"),
				 panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
								 colour = "black"), 
				 panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
								 colour = "black")
				 )+
		   guides(color = guide_legend(override.aes = list(size = 7.5)))+
		   theme(legend.title = element_text(size=12,face='bold',color='black'),
			 legend.text = element_text(face='bold',color='black'))+ggtitle(paste0('CRC -',unique(x$fov)))
		   g})

pdf(plots,width=9,height=9)
for(p in plot_list){
	print(p)
}
dev.off()



