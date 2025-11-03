library(arrow)
library(sf)
library(sfarrow)
library(dplyr)
library(stringr)
library(BiocParallel)
source('../local/src/df2sf.R')
#Define samples for which metadata must be prepared
slide = snakemake@params[["slide"]]
if(slide == 'S0'){
	samples = c('CR48','UH20')
}else{
	samples = c('CR36','CR21')
}
#Define input files path
meta = list.files(snakemake@params[["input_path"]],full.names=T,pattern='metadata')
poly = list.files(snakemake@params[["input_path"]],full.names=T,pattern='polygon')

#merge metadata and polygon files for each sample, after having converted polygon coordinates into an sf object
all = lapply(samples,function(x){
	     meta = read_csv_arrow(meta[str_detect(meta,x)])
	     poly = read_csv_arrow(poly[str_detect(poly,x)])
	     poly$ID = poly$cell
	     meta$sample = x
	     cells = unique(meta$cell_id)
	     poly = poly %>% filter(ID %in% cells)
	     poly = df2sf(poly,spatialCoordsNames = c('x_global_px','y_global_px'))
	     poly$cell_id = poly$ID
	     merged = merge(meta,poly,by='cell_id')
	     global_centr = st_centroid(merged$geometry) %>% st_coordinates()
	     merged$x_global_px = global_centr[,1]
	     merged$y_global_px = global_centr[,2]
	     merged = st_as_sf(merged)
	     return(merged)
})
#Save final table as a sf data frame in parquet format (with sfarrow package)
st_write_parquet(all[[1]],snakemake@output[["merged_data1"]])
st_write_parquet(all[[2]],snakemake@output[["merged_data2"]])
rm(all)
