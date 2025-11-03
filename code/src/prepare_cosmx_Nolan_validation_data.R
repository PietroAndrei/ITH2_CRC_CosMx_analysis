library(Seurat)
library(dplyr)

path = snakemake@params[["data_path"]]
seurat = list.files(path,pattern='rds',full.names=T)
meta = list.files(path,pattern='csv',full.names=T)

seurat = readRDS(seurat)
seurat@misc[names(seurat@misc)] = NULL
seurat@graphs[names(seurat@graphs)] = NULL
DefaultAssay(seurat) = 'RNA'
seurat[['negprobes']] = NULL
seurat[['SCT']] = NULL
meta = read.delim(meta,sep=',')
all_meta = merge(seurat@meta.data,meta,by.x = 'fov',by.y='FOV_TAP')
rownames(all_meta) = all_meta$cell_ID
all_meta = all_meta %>% filter(indication == 'CRC')
seurat = seurat[,WhichCells(seurat,cells=rownames(all_meta))]
seurat@meta.data = all_meta
saveRDS(seurat,snakemake@output[["cosmx_clean"]])
