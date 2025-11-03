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
seurat = FindVariableFeatures(seurat)
neg = round(seurat$nCount_negprobes/10)

immunof = seurat@meta.data[,str_detect(colnames(seurat@meta.data),'Mean\\.|Area.um2|AspectRatio')]
cohort <- fastCohorting(immunof,
			gaussian_transform = TRUE) 

#For InSituType, limit the genes to the top 150 most variable ones
mat = Matrix::t(LayerData(seurat,layer='counts'))
#Round protein intensity to counts
mat = round(mat)
#Perform InSituType unsupervised clustering (NB* expr matrix must have cell ids in rows an gene ids in columns)
unsup <- insitutype(
		    x = mat,
		    neg = neg,
		    cohort = cohort,
		    # Enter your own per-cell background estimates here if you have them;
		    # otherwise insitutype will use the negprobes to estimate background for you.
		    bg = NULL,
		    n_starts = 10,
		    max_iters = 30,
		    # condensed to save time. n_clusts = 5:20 would be more optimal
		    n_clusts = 10:20
		    )

seurat = NormalizeData(seurat) %>% 
ScaleData() %>% 
RunPCA(npcs = 30) %>% 
RunUMAP(reduction = 'pca',dims=1:25)
seurat$InSitu_clusters = unsup$clust
Idents(seurat) = 'InSitu_clusters'


fc.func = function(x=NULL,pscount=1,base=2){return(log(x = rowMeans(x = expm1(x = x)) + pscount, base = base))}
markers = FindAllMarkers(seurat,only.pos = T,logfc.threshold = log2(1.5),mean.fxn = fc.func)
saveRDS(markers,snakemake@output[["markers"]])

meta_sf = seurat@meta.data
meta_sf$x_slide_um = meta_sf$x_slide_mm*1000
meta_sf$y_slide_um = meta_sf$y_slide_mm*1000
meta_sf = make_geometry_df(meta_sf,xcoord='x_slide_um',ycoord='y_slide_um')


meta_list = split(meta_sf,meta_sf$fov)
plots = snakemake@output[["plots"]]
plot_list = lapply(meta_list,function(x){
		   g = x %>% sf::st_as_sf() %>%
		   ggplot()+geom_sf(aes(fill=InSitu_clusters),lwd=0.01,shape=21,size=2)+
		   scale_fill_manual(values=c24)+labs(fill='Cell Type')
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


cell_anno = snakemake@output[["cell_anno"]]
saveRDS(seurat@meta.data %>% select(cell_ID,InSitu_clusters),cell_anno)

InSitu_res = snakemake@output[["InSitu_unsup_res"]]
saveRDS(unsup,InSitu_res)
