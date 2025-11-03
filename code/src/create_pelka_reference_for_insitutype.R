library(Seurat)
library(dplyr)
library(stringr)
library(Matrix)
utils = snakemake@params[["utils"]]
source(utils)
cosmx_genes = readRDS(snakemake@params[["cosmx_genes"]])
#Set the directory containing Pelka related files
pelka_path = snakemake@params[["pelka_path"]]
pelka = list.files(pelka_path,full.names=T)
#Set the directory containing Pelka related files
h5 = pelka[str_detect(pelka,'h5')]
#h5 = paste0(pelka_path,h5)
#Load Pelka metadata and clusters
meta = pelka[str_detect(pelka,'metatables')]
meta = read.delim(meta,sep=',')

obs = pelka[str_detect(pelka,'cluster')]
obs = read.delim(obs,sep=',')

#Keep only immune/stromal cells (including those derived from normal adjacent tissue)
all_cells = obs[(obs$clTopLevel != 'Epi'),]$sampleID
message("Considering ",length(all_cells)," immune/stromal cells from Pelka dataset")
rownames(meta) = meta$cellID
rownames(obs) = obs$sampleID
meta = meta[all_cells,]
obs = obs[all_cells,]
meta$origCluster = obs$clMidwayPr
rm(obs)

#Load Pelka dataset
data = Read10X_h5(filename = h5, use.names = TRUE, unique.features = TRUE)
bcs = colnames(data)

#Filter Pelka dataset
common_bcs = intersect(rownames(meta),bcs)
message(length(all_cells),"/",length(common_bcs)," single cells matching with metadata in Pelka dataset")
data = data[,common_bcs]
meta = meta[common_bcs,]
genes = rownames(data)



#Load NCG7 gene list that will be used to filter the count matrix and update the gene symbols
ncg = snakemake@params[['ncg']]
ncg = read.delim(ncg)
#Update gene symbols in the dataset when necessary
new_symbols = sapply(genes,function(x){
	is_new = x %in% unique(ncg$NCG_symbol)
	is_old = !(x %in% unique(ncg$NCG_symbol)) & (x %in% unique(ncg$symbol))
	if(is_new == F & is_old == T){
		x = unique(ncg[which(ncg$symbol == x),]$NCG_symbol)
	} else {
		x = x
	}
	return(x)
}
)

#Filter Lee matrix to keep only the selected NCG genes
data = data[names(new_symbols),]
#Update gene symbols
rownames(data) = new_symbols
message("Keeping ",nrow(data)," genes included in the NCG7 catalogue")

#Adapt gene names to cosmx 1k plex format and keep only the overlapping genes between pelka and cosmx
data = gex2cosmx(data)
data = data[intersect(rownames(data),cosmx_genes),]

#Save matrix as Seurat objects
seurat = CreateSeuratObject(counts = data,meta.data = meta)

#Sum all counts across cells of the same cell type (to create a reference profile)
ref = AggregateExpression(seurat,group.by='origCluster',return_seurat=F)
ref = ref[[1]]
#Divide the total counts by the number of cells available for each cell type
ncells = table(seurat$origCluster)
cts = names(ncells)
ncells = as.vector(ncells)
names(ncells) = cts
ncells = ncells[colnames(ref)]
ref = ref %*% Diagonal(x = 1/ncells)
colnames(ref) = names(ncells)
ref = as.matrix(ref)
saveRDS(ref,snakemake@output[["ref_mat"]])
