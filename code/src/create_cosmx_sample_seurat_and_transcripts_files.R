library(Seurat)
library(arrow)
library(sf)
library(sfarrow)
library(dplyr)
utils = snakemake@params[["utils"]]
source(utils)
options(Seurat.object.assay.version = 'v5')
#Define paths for sample-specific metadata/polygon files
meta_path = snakemake@params[["meta_path"]]
meta = list.files(meta_path)

#Load merged Seurat object (one of the main output downloaded from AtoMx)
seurat = snakemake@input[["seurat"]]
seurat = readRDS(seurat)

#Extract a Seurat object for each sample
samples = lapply(meta,function(x){
       file = paste0(meta_path,x)
       file = st_read_parquet(file)
       rownames(file) = file$cell_ID
       DefaultAssay(seurat) = 'RNA'
       seurat[['negprobes']] = NULL
       seurat[['falsecode']] = NULL
       seurat = seurat[,rownames(file)]
       seurat@meta.data = as.data.frame(file)
       return(seurat)
})
rm(seurat)
#Save sample Seurat objects
seurat_path = snakemake@params[["seurat_path"]]
if(!file.exists(seurat_path)){
system(paste0('mkdir ',seurat_path))
}
for(s in 1:length(samples)){
	sample_id = unique(samples[[s]]$sample)
	slide_id = paste0('Slide',unique(samples[[s]]$slide_ID_numeric))
	out_file = paste0('Cosmx_',slide_id,"_",sample_id,'_SeuratBase.rds')
	saveRDS(object=samples[[s]],file=paste0(seurat_path,out_file))
}
rm(samples)
for(i in 1:10){gc()}
