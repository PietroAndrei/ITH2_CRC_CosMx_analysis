library(Seurat)
library(sf)
library(MAST)
library(dplyr)
library(stringr)
utils = snakemake@params[["utils"]]
source(utils)
options(Seurat.object.assay.version = 'v5')

meta_path = snakemake@params[["meta_path"]]


samples = c('CR36','UH20','CR48','CR21')

all_meta = list.files(meta_path)
all_meta = all_meta[str_detect(all_meta,paste0(samples,collapse='|'))]

all = lapply(samples,function(x){
	     meta = readRDS(paste0(
				   meta_path,
				   all_meta[str_detect(all_meta,x)]
				   )
	     )
	     meta = meta %>% filter(final_anno == 'Myeloid')
	     meta$myelType = ifelse(meta$nb_tnk > 0,'myelT','myelNoT')
	     message('Comparing ',nrow(meta %>% filter(myelType == 'myelT')),
		     ' myelT vs ',nrow(meta %>% filter(myelType == 'myelNoT')),' myelNoT cells')
	     rownames(meta) = meta$cell_ID
	     return(meta)
})
all = purrr::reduce(all,rbind)

seurat = readRDS(snakemake@input[["seurat"]])
DefaultAssay(seurat) = 'RNA'
seurat[['negprobes']] = NULL
seurat[['falsecode']] = NULL
seurat = seurat[,WhichCells(seurat,cells=rownames(all))]
seurat@meta.data = as.data.frame(all)
rm(all)

Idents(seurat) = 'cell_ID'
LayerData(seurat,'data')=LayerData(seurat,'counts')
seurat = NormalizeData(seurat)
Idents(seurat) = 'myelType'

#assignInNamespace('MASTDETest', MASTDE_RandTest, asNamespace("Seurat"))
fc.func = function(x=NULL,pscount=1,base=2){return(log(x = rowMeans(x = expm1(x = x)) + pscount, base = base))}
markers = FindMarkers(seurat,ident.1='myelT',
		      only.pos=F,logfc.threshold = -Inf,
		      mean.fxn=fc.func,
		      test.use='MAST',latent.vars='sample')#,re.var='sample',ebayes=F)

de_res = snakemake@output[["de_res"]]
saveRDS(markers,de_res)
