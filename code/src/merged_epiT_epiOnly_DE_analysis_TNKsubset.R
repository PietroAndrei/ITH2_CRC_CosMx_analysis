library(Seurat)
library(sf)
library(MAST)
library(dplyr)
library(stringr)
utils = snakemake@params[["utils"]]
source(utils)
options(Seurat.object.assay.version = 'v5')

HM = snakemake@params[["HM"]]
meta_path = snakemake@params[["meta_path"]]
n_neigh = snakemake@params[["n_neigh"]]
n_epit = snakemake@params[["n_epit"]]

if(HM == 'HM'){
	samples = c('CR36','UH20')
}
if(HM == 'nHM'){
	samples = c('CR48','CR21')
}
if(HM == 'ALL'){
	samples = c('CR48','CR21','CR36','UH20')
}
all_meta = list.files(meta_path)
all_meta = all_meta[str_detect(all_meta,paste0(samples,collapse='|'))]



all = lapply(samples,function(x){
	     meta = readRDS(paste0(
				   meta_path,
				   all_meta[str_detect(all_meta,x)]
				   )
	     )
	     #Define all touching cell-cell pairs and identify epi-T cells
	     #Consider only T/NK cells in touch with >=4 epithelial cells (higher chance to be real intra-epithelial t/nk cells)
	     nb = st_intersects(meta)
	     nb = lapply(seq_len(length(nb)),function(x) nb[[x]][nb[[x]] != x])
	     attr(nb,'class') = c('sgbp','list')
	     nb = as.data.frame(nb)
	     nb$row.type = as.character(meta[nb$row.id,]$final_anno)
	     nb$col.type = as.character(meta[nb$col.id,]$final_anno)
	     nb$row.cellID = meta[nb$row.id,]$cell_ID
	     nb$col.cellID = meta[nb$col.id,]$cell_ID
	     tcells = meta %>% filter(final_anno == 'T/NK',nb_epi >= n_epit)
	     #Keep as epi-T only epithelial cells in touch with the previously defined intra-epithelial t/nk cells
	     nb_epit = nb %>% filter(col.type == 'Epi',row.cellID %in% tcells$cell_ID) 
	     epit = meta %>% filter(cell_ID %in% unique(nb_epit$col.cellID))
	     #epiOnly: in touch with >=4 other epithelial cells (higher chance to be immune-excluded epithelial cells)
	     meta = meta %>% 
	     filter(
		    (epiType == 'epiT' & cell_ID %in% epit$cell_ID) | 
		    (epiType == 'epiOnly' & n_neighbors >= n_neigh)
		    )
	     message('Comparing ',nrow(meta %>% filter(epiType == 'epiT')),
		     ' epiT vs ',nrow(meta %>% filter(epiType == 'epiOnly')),' epiOnly cells')
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
Idents(seurat) = 'epiType'

#assignInNamespace('MASTDETest', MASTDE_RandTest, asNamespace("Seurat"))
fc.func = function(x=NULL,pscount=1,base=2){return(log(x = rowMeans(x = expm1(x = x)) + pscount, base = base))}
markers = FindMarkers(seurat,ident.1='epiT',
		      only.pos=F,logfc.threshold = -Inf,
		      mean.fxn=fc.func,
		      test.use='MAST',latent.vars='sample')#,re.var='sample',ebayes=F)

de_res = snakemake@output[["de_res"]]
saveRDS(markers,de_res)
