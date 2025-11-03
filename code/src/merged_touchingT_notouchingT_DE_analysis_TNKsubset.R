library(Seurat)
library(sf)
library(MAST)
library(dplyr)
library(stringr)
utils = snakemake@params[["utils"]]
source(utils)

HM = snakemake@params[["HM"]]
meta_path = snakemake@params[["meta_path"]]

if(HM == 'HM'){
	samples = c('CR36','UH20')
}
if(HM == 'nHM'){
	samples = c('CR48','CR21')
}
all_meta = list.files(meta_path)
all_meta = all_meta[str_detect(all_meta,paste0(samples,collapse='|'))]

all = lapply(samples,function(x){
	     meta = readRDS(paste0(
				   meta_path,
				   all_meta[str_detect(all_meta,x)]
				   )
	     )
	     t = meta %>% filter(final_anno == 'CD8+/NK' | final_anno == 'T/NK')
	     t_touching = t %>% filter(nb_epi >=3) %>% mutate(tType = "t_touch")
	     t_notouching = t %>% filter(nb_epi == 0 ) %>% mutate(tType = "t_notouch")
	     rm(meta)
	     meta = rbind(t_touching,t_notouching)
	     message('Comparing T cells in touch (',nrow(meta %>% filter(tType == 't_touch')),
		     ') vs not in touch (',nrow(meta %>% filter( tType == 't_notouch')),') with Epithelial cells')
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
LayerData(seurat,'data') = LayerData(seurat,'counts')
seurat = NormalizeData(seurat)
Idents(seurat) = 'tType'

#Correct MAST pvalue by including Mean.PanCK value as a confounding factor
fc.func = function(x=NULL,pscount=1,base=2){return(log(x = rowMeans(x = expm1(x = x)) + pscount, base = base))}
markers = FindMarkers(seurat,ident.1='t_touch',latent.vars=c('sample','Mean.PanCK'),only.pos=F,logfc.threshold = -Inf,mean.fxn=fc.func,test.use='MAST')
de_res = snakemake@output[["de_res"]]
saveRDS(markers,de_res)
