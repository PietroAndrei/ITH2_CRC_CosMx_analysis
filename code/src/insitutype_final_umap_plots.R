library(Seurat)
library(sf)
library(dplyr)
library(stringr)
library(ggplot2)
utils = snakemake@params[["utils"]]
source(utils)
sample = snakemake@wildcards[["sample"]]
seurat_path = 'Seurat/'
meta_path = snakemake@params[["meta_path"]]
#Load sample metadata
meta = list.files(meta_path)
meta = meta[str_detect(meta,sample)]
meta = readRDS(paste0(meta_path,meta))
#Define cell types and assign a color to each of them
celltypes = c('Epi','FibroEndoMuscle','Myeloid','Mast','Plasma/B','T/NK','T_Other')
cols = c('#0000FF','#E52628','#43CD80','#A2B5CD','#00FFFF','#FFFF00','#FFE4B5')
names(cols) = celltypes
ncells = table(meta$final_anno) %>% as.data.frame()
rownames(ncells) = ncells$Var1
colnames(ncells) = c('celltypes','n')

#Load seurat file
seurat = list.files(seurat_path)
seurat = seurat[str_detect(seurat,sample)]
seurat = readRDS(paste0(seurat_path,seurat))
#update seurat metadata
seurat = seurat[,WhichCells(seurat,cells=rownames(meta))]
seurat@meta.data = as.data.frame(meta)
Idents(seurat) = 'cell_ID'
#count normalization+dimensionality reduction
seurat = NormalizeData(seurat)
seurat = seurat %>%
FindVariableFeatures() %>%
ScaleData() %>%
RunPCA(npcs = 30) %>% 
RunUMAP(reduction='pca',dims=1:25)
Idents(seurat) = 'final_anno'
#Save PCA and UMAP embeddings
pca = Embeddings(seurat,reduction='pca')
umap = Embeddings(seurat,reduction='umap')
#Plot umap
p = DimPlot(seurat,reduction='umap',pt.size=1.5,cols=cols,raster.dpi=c(1270,1270))+
theme(axis.text = element_blank(),axis.ticks=element_blank(),legend.title=element_blank())+
scale_color_discrete(labels=paste0(levels(meta$final_anno)," (n = ",ncells[levels(meta$final_anno),]$n,")"),
		     type=cols[levels(meta$final_anno)])+
ggtitle(paste0(sample," (n = ",nrow(meta),")"))

rm(seurat)
pdf(snakemake@output[["plot"]],width=13,height=9)
p
dev.off()
rm(p)

saveRDS(list(pca=pca,umap=umap),snakemake@output[["dimreds"]])
