library(Seurat)
library(presto)
library(InSituType)
library(sf)
library(dplyr)
library(stringr)
library(ggplot2)
utils = snakemake@params[["utils"]]
source(utils)

#Use InSituType to recluster TNK cells separately, then update the global cell annotation
#NB*!! 'to_delete' parameter in refineClusters function does not remove cells, it will just reassign them to the closest cluster. Remove bad clusters in Seurat

sample = snakemake@wildcards[["sample"]]
seurat = readRDS(snakemake@input[["seurat"]])
seurat = filter_cosmx_cells(seurat)
seurat = FindVariableFeatures(seurat,nfeatures=1000)
Idents(seurat) = 'cell_ID'
neg = seurat$nCount_negprobes/10

immunof = seurat@meta.data[,str_detect(colnames(seurat@meta.data),'Mean\\.')]
cohort <- fastCohorting(immunof,
			gaussian_transform = TRUE) 

#(NB* expr matrix must have cell ids in rows an gene ids in columns)
mat = Matrix::t(LayerData(seurat,layer='counts'))

unsup_res = snakemake@input[["unsup_res"]]
unsup_res = readRDS(unsup_res)

anno_raw = unsup_res$clust


if(sample == 'CR36'){
newclusts <- refineClusters(logliks = unsup_res$logliks,
			    merges = c('c' = 'Epi','n' = 'Epi','l'='Epi','m'='Epi','k'='Epi',
				       'b' = 'TNK','e' = 'TNK','d'='Plasma/B','f'='Plasma/B',
				       'i' = 'Myeloid','g' = 'Myeloid',
				       'a' = 'FibroEndoMuscle','j' = 'FibroEndoMuscle','h'='FibroEndoMuscle')
			    )
anno_aggr = newclusts$clust
newclusts = refineClusters(logliks = newclusts$logliks,
			   subcluster = list(TNK = 2:4),
			   counts =mat,
			   cohort = cohort,
			   neg = neg)
}
if(sample == 'CR21'){
newclusts <- refineClusters(logliks = unsup_res$logliks,
			    merges = c('m' = 'Epi','k'='Epi','j'='TNK','e'='Epi','a'='FibroEndoMuscle',
				       'h' = 'FibroEndoMuscle','g'='Myeloid','c' = 'Epi','b'='Epi','i'='FibroEndoMuscle',
				       'd'='FibroEndoMuscle','f'='Plasma/B','l'='Epi')
			    )
anno_aggr = newclusts$clust
newclusts = refineClusters(logliks = newclusts$logliks,
			   subcluster = list(TNK = 2:4),
			   counts =mat,
			   cohort = cohort,
			   neg = neg)
}
if(sample == 'CR48'){
newclusts <- refineClusters(logliks = unsup_res$logliks,
			    merges = c('a' = 'FibroEndoMuscle','k'='Epi','m'='Myeloid','n'='Epi','e'='Epi',
				       'j'='TNK','i'='Epi','g'='Myeloid','d'='FibroEndoMuscle','c'='FibroEndoMuscle',
				       'l' = 'FibroEndoMuscle','f'='Epi','h'='FibroEndoMuscle','b'='Plasma/B')
			    )
anno_aggr = newclusts$clust
#Try to subset Lymphoid cells
newclusts = refineClusters(logliks = newclusts$logliks,
			   subcluster = list(TNK = 2:4),
			   counts =mat,
			   cohort = cohort,
			   neg = neg)
}
if(sample == 'UH20'){
newclusts <- refineClusters(logliks = unsup_res$logliks,
			    merges = c('d' = 'Epi','f'='Myeloid','c'='FibroEndoMuscle','m'='Epi','k'='Myeloid',
				       'l' = 'FibroEndoMuscle','n'='Epi','j'='FibroEndoMuscle','g'='Plasma/B',
				       'i'='TNK','b'='Myeloid','e'='FibroEndoMuscle','a'='Epi','h'='Epi')
			    )
anno_aggr = newclusts$clust
newclusts <- refineClusters(logliks = newclusts$logliks,
			    subcluster = list(TNK = 2:4),
			    counts =mat,
			    cohort = cohort,
			    neg = neg)
}

rm(mat)
all_celltypes = c('Epi','FibroEndoMuscle','Myeloid','T/NK','T_Other','Plasma/B','Mast')
colscheme = c24[1:length(all_celltypes)]
names(colscheme) = all_celltypes
#Prepare plots for the raw and aggregate cell typing
seurat = NormalizeData(seurat) %>% 
ScaleData() %>% 
RunPCA(npcs = 30) %>% 
RunUMAP(reduction = 'pca',dims=1:25)
seurat$InSitu_clusters_raw = anno_raw
seurat$InSitu_clusters_aggr = anno_aggr

all_cells_anno = seurat@meta.data[,c('cell_ID','InSitu_clusters_raw','InSitu_clusters_aggr')]

seurat$InSitu_clusters_raw = factor(seurat$InSitu_clusters_raw,levels = sort(unique(seurat$InSitu_clusters_raw)))
seurat$InSitu_clusters_aggr = factor(seurat$InSitu_clusters_aggr,levels = sort(unique(seurat$InSitu_clusters_aggr)))
p1 = DimPlot(seurat,cols=c24,group.by='InSitu_clusters_raw',raster=F)+theme_void()
p2 = DimPlot(seurat,cols=c24,group.by='InSitu_clusters_aggr',raster=F)+theme_void()

#Add final cell annotation as final cell identity. If some cells have been removed in the final refinement step, adjust PCA and UMAP
seurat$final_anno = newclusts$clust
if(sample == 'CR48'){
	seurat@meta.data[seurat@meta.data$final_anno %in% c('TNK_1','TNK_2','TNK_3'),]$final_anno = 'T/NK'
	seurat = seurat[,seurat$final_anno != 'TNK_4']
}
if(sample == 'CR21'){
	seurat@meta.data[seurat@meta.data$final_anno == 'TNK_1',]$final_anno = 'Myeloid'
	seurat@meta.data[seurat@meta.data$final_anno == 'TNK_3',]$final_anno = 'T/NK'
	seurat@meta.data[seurat@meta.data$final_anno == 'TNK_4',]$final_anno = 'Mast'
	seurat = seurat[,seurat$final_anno != 'TNK_2']
}
if(sample == 'CR36'){
	seurat@meta.data[seurat@meta.data$final_anno %in% c('TNK_1','TNK_2','TNK_3'),]$final_anno = 'T/NK'
	seurat@meta.data[seurat@meta.data$final_anno == 'TNK_4',]$final_anno = 'T_Other'
}
if(sample == 'UH20'){
	seurat@meta.data[seurat@meta.data$final_anno %in% c('TNK_1','TNK_2','TNK_3'),]$final_anno = 'T/NK'
	seurat = seurat[,seurat$final_anno != 'TNK_4']
}

message(length(anno_raw) - ncol(seurat)," cells have been removed after T/NK cell subclustering analysis")
seurat$final_anno = factor(seurat$final_anno,levels = all_celltypes)
Idents(seurat) = 'final_anno'

if((length(anno_raw) - ncol(seurat)) > 0){
	seurat = seurat %>%
	ScaleData() %>% 
	RunPCA(npcs = 30) %>% 
	RunUMAP(reduction = 'pca',dims=1:25)
}

p3 = DimPlot(seurat,cols=c24,raster=F)+theme_void()

fc.func = function(x=NULL,pscount=1,base=2){return(log(x = rowMeans(x = expm1(x = x)) + pscount, base = base))}
markers = FindAllMarkers(seurat,only.pos = T,logfc.threshold = log2(1.5),min.diff.pct = 0.1,mean.fxn = fc.func)
saveRDS(markers,snakemake@output[["markers"]])


plots = snakemake@output[["plots"]]
pdf(plots,width=13,height=9)
p1+labs(color='Cell Type')+theme(plot.title = element_blank())
p2+labs(color='Cell Type')+theme(plot.title = element_blank())
p3+labs(color='Cell Type')+theme(plot.title = element_blank())
cluster_fov_plot(seurat@meta.data,class='final_anno',colors=colscheme)
dev.off()
rm(p1)
rm(p2)
rm(p3)

cell_anno = snakemake@output[["cell_anno"]]
final_cells_anno = seurat@meta.data[,c('cell_ID','final_anno')]
final_cells_anno = merge(x=all_cells_anno,y=final_cells_anno,by='cell_ID',all.x=T)
saveRDS(final_cells_anno,cell_anno)

InSitu_res = snakemake@output[["InSitu_refined_res"]]
saveRDS(newclusts,InSitu_res)
rm(newclusts)
rm(unsup_res)
#Define EpiT/EpiOnly epithelial groups on the basis of the touching cells around them
touching_cells = get_touching_cells(seurat@meta.data,label='final_anno')
touching_cells$is_epiT = 0
touching_cells$nb_tnk = str_split(touching_cells$neighbors_label,',') %>% str_count(.,'T\\/NK')
touching_cells = touching_cells %>% mutate(is_epiT = ifelse(final_anno == 'Epi' & nb_tnk > 0,1,0))
touching_cells$is_epionly = 0
touching_cells$nb_epi = str_split(touching_cells$neighbors_label,',') %>% str_count(.,'Epi')
touching_cells = touching_cells %>% mutate(is_epionly = ifelse(nb_epi == n_neighbors & n_neighbors > 0 & final_anno == 'Epi',1,0))
touching_cells = touching_cells %>% mutate(epiType = ifelse(
							    is_epiT == 1,'epiT',
							    ifelse(
								   is_epionly == 1,'epiOnly',
								   ifelse(
									  final_anno != 'Epi','StromaImmune','epiOther'
									  )
								   )
							    )
)
touching_cells$nb_noepi = str_count(touching_cells$neighbors_label,paste0(all_celltypes[all_celltypes != 'Epi'],collapse='|'))

final_meta = snakemake@output[["final_meta"]]
saveRDS(touching_cells,final_meta)
