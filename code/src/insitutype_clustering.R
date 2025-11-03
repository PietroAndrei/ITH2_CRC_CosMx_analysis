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
seurat = filter_cosmx_cells(seurat)
seurat = FindVariableFeatures(seurat,nfeatures=150)
var = VariableFeatures(seurat)
Idents(seurat) = 'cell_ID'
neg = seurat$nCount_negprobes/10

immunof = seurat@meta.data[,str_detect(colnames(seurat@meta.data),'Mean\\.')]
cohort <- fastCohorting(immunof,
			gaussian_transform = TRUE) 

#For InSituType, limit the genes to the top 150 most variable ones
mat = Matrix::t(LayerData(seurat,layer='counts'))

#Perform InSituType unsupervised clustering (NB* expr matrix must have cell ids in rows an gene ids in columns)
unsup <- insitutype(
		    x = mat,
		    neg = neg,
		    cohort = cohort,
		    # Enter your own per-cell background estimates here if you have them;
		    # otherwise insitutype will use the negprobes to estimate background for you.
		    bg = NULL,
		    n_starts = 7,
		    max_iters = 20,
		    # condensed to save time. n_clusts = 5:20 would be more optimal
		    n_clusts = 10:14
		    )

seurat = NormalizeData(seurat) %>% 
ScaleData() %>% 
RunPCA(npcs = 30) %>% 
RunUMAP(reduction = 'pca',dims=1:25)
seurat$InSitu_clusters = unsup$clust
Idents(seurat) = 'InSitu_clusters'


fc.func = function(x=NULL,pscount=1,base=2){return(log(x = rowMeans(x = expm1(x = x)) + pscount, base = base))}
markers = FindAllMarkers(seurat,only.pos = T,logfc.threshold = log2(1.5),min.diff.pct = 0.1,mean.fxn = fc.func)
saveRDS(markers,snakemake@output[["markers"]])


plots = snakemake@output[["plots"]]
pdf(plots,width=13,height=9)
DimPlot(seurat,cols=c24,raster=F)+theme_void()
cluster_fov_plot(seurat@meta.data,class='InSitu_clusters',colors=c24)
dev.off()

cell_anno = snakemake@output[["cell_anno"]]
saveRDS(seurat@meta.data %>% select(cell_ID,InSitu_clusters),cell_anno)

InSitu_res = snakemake@output[["InSitu_unsup_res"]]
saveRDS(unsup,InSitu_res)
