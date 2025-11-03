library(Seurat)
library(InSituType)
library(sf)
library(dplyr)
library(stringr)
library(ggplot2)
utils = snakemake@params[["utils"]]
source(utils)

sample = snakemake@wildcards[["sample"]]
meta_path = snakemake@params[["meta_path"]]
meta = list.files(meta_path)
meta = meta[str_detect(meta,sample)]

meta = readRDS(paste0(meta_path,meta))

#Select only T/NK cells
meta = meta %>% filter(final_anno == 'Myeloid')
seurat_path = 'Seurat/'
seurat = list.files(seurat_path)
seurat = seurat[str_detect(seurat,sample)]
seurat = readRDS(paste0(seurat_path,seurat))

seurat= seurat[,WhichCells(seurat,cells = rownames(meta))]

#Extract negative probes quantification for insitutype normalization and immunofluorescence values for cohorting
neg = seurat$nCount_negprobes/10
meta = as.data.frame(meta)
immunof = meta[,str_detect(colnames(meta),'Mean\\.')]

cohort <- fastCohorting(immunof,
			gaussian_transform = TRUE) 

#Extract raw count matrix 
mat = Matrix::t(LayerData(seurat,layer='counts'))
rm(seurat)
#Load reference pelka matrix and keep only relevant cell types
ref_mat = readRDS(snakemake@params[["ref_mat"]])
ref_mat = ref_mat[,c('Macro','Mono','DC','Granulo','Mast')]

#Run InSituType
sup = insitutypeML(x = mat,
		   neg = neg,
		   cohort = cohort,
		   reference_profiles = ref_mat)
saveRDS(sup,snakemake@output[["res"]])


fp_layout = flightpath_layout(logliks = sup$logliks, profiles = sup$profiles)
cols = c12[seq_along(unique(sup$clust))]
names(cols) = unique(sup$clust) 

#Plot clustering results with a flightpath plot ( it shows the average probability for each cluster)
fp = flightpath_plot(flightpath_result = fp_layout, insitutype_result = NULL, col = cols[sup$clust])
pdf(snakemake@output[["plot"]],width=13,height=8)
fp
dev.off()


