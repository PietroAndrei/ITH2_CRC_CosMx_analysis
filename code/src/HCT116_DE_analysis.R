library(DESeq2)
library(dplyr)
library(stringr)
#Load count matrix and keep only NCG7 genes
ncg7 = snakemake@params[["ncg"]]
mat = read.delim(snakemake@input[["mat"]])
mat$gene_name = NULL
ncg7 = read.delim(ncg7)
colnames(mat)[colnames(mat) == "gene_id"] = "symbol"
ncg7_mat = merge(mat, ncg7, by = "symbol")
ncg7_mat = select(ncg7_mat, -NCG_symbol, -entrez)
rownames(ncg7_mat) = ncg7_mat$symbol
ncg7_mat$symbol = NULL

message('Number of NCG7 genes: ',nrow(ncg7_mat))

meta = data.frame(sample = colnames(ncg7_mat),condition = c(rep('EpiT',3),rep('EpiOnly',3)))
meta$condition = factor(meta$condition,levels=c('EpiOnly','EpiT'))
rownames(meta) = meta$sample
ncg7_mat = ncg7_mat %>% mutate_all(as.integer)
dds = DESeqDataSetFromMatrix(ncg7_mat,colData=meta,~condition)

##Keep only genes expressed in at least one sample (>10 counts) and perform DE test
keep = which(rowSums(assay(dds) > 10) > 0)
dds = dds[keep,]
message('Number of expressed genes considered for following analyses: ',length(keep))
dds = DESeq(dds)

#reorder DE results and save
res = results(dds) %>% as.data.frame()
res$gene = rownames(res) 
res = res %>% arrange(desc(stat))
saveRDS(dds,snakemake@output[["dds"]])
saveRDS(res,snakemake@output[["de_res"]])
