library(sf)
library(dplyr)
library(stringr)
library(ggplot2)
utils = snakemake@params[["utils"]]
source(utils)

#Load UH20 metadata and select a prototype region to highlight cell types of interest
cosmx_uh20_meta = readRDS('Seurat_meta/InSituType/V2/Cosmx_Slide1_UH20_InSitu_celltype_plus_touching_cells_V2_meta.rds')
cosmx_proto = cosmx_uh20_meta %>% filter(fov %in% c(7,8,9,10,11,12,51,176,178))

#Define the list of epithelial/T/macrophage cell types on the basis of their cell contact
nb = st_intersects(cosmx_proto)
nb = lapply(seq_len(length(nb)),function(x) nb[[x]][nb[[x]] != x]) 
attr(nb,'class') = c('sgbp','list')
nb = as.data.frame(nb)
nb$row.type = as.character(cosmx_proto[nb$row.id,]$final_anno)
nb$col.type = as.character(cosmx_proto[nb$col.id,]$final_anno)
nb$row.cellID = cosmx_proto[nb$row.id,]$cell_ID
nb$col.cellID = cosmx_proto[nb$col.id,]$cell_ID
tcells = cosmx_proto %>% filter(final_anno == 'T/NK')
ielt = cosmx_proto %>% filter(final_anno == 'T/NK',nb_epi >= 3)
stromalt = cosmx_proto %>% filter(final_anno == 'T/NK',nb_epi == 0)
#Keep as epi-T only epithelial cells in touch with the previously defined intra-epithelial t/nk cells
nb_epit = nb %>% filter(col.type == 'Epi',row.cellID %in% ielt$cell_ID) 
epit = cosmx_proto %>% filter(cell_ID %in% unique(nb_epit$col.cellID))
epionly = cosmx_proto %>% filter(epiType == 'epiOnly' & n_neighbors >= 4)
myelt = cosmx_proto %>% filter(final_anno == 'Myeloid',nb_tnk >= 1)
myelother = cosmx_proto %>% filter(final_anno == 'Myeloid',nb_tnk == 0)

#Create annotation to highlight TAM and T cell contacts
cosmx_proto = cosmx_proto %>% mutate(new_anno = ifelse( cell_ID %in% myelt$cell_ID,'TAM-T',ifelse(cell_ID %in% myelother$cell_ID,'TAM-Other',ifelse(final_anno == 'T/NK','T','Other'))))
#Create annotation to highlight Epi ant T cell contacts
cosmx_proto = cosmx_proto %>% mutate(new_anno2 = ifelse(cell_ID %in% ielt$cell_ID,'IE-T',ifelse(cell_ID %in% stromalt$cell_ID ,'Ts-T',ifelse(cell_ID %in% epit$cell_ID,'Epi-T',ifelse(cell_ID %in% epionly$cell_ID,'Epi-Only','Other')))))

cosmx_proto$new_anno = factor(cosmx_proto$new_anno,levels=c('TAM-T','TAM-Other','T','Other'))
cosmx_proto$new_anno2 = factor(cosmx_proto$new_anno2,levels=c('IE-T','Ts-T','Epi-T','Epi-Only','Other'))

#Plot a specfic region where all the contact events of interest are represented (after visual inspection of CosMx FOVs)
pdf(snakemake@output[["UH20_FOV11_TAMs_T"]],width=2.5,height=2.5)
cluster_fov_plot(cosmx_proto %>% filter(fov == 11)%>% filter(x_global_px >= 43830,x_global_px <= 45830,y_global_px >= 120700),lwd=0.22,class='new_anno',colors = c('#83C65D','#3B4E36','#FF1AFF','#A7B0BD'))+
theme(axis.text = element_blank(),axis.ticks = element_blank(),plot.title = element_blank(),
      legend.key.size = unit(0.2, 'cm'),
      legend.key.height = unit(0.2, 'cm'),
      legend.key.width = unit(0.2, 'cm'),
      legend.title = element_text(size=6),
      legend.text = element_text(size=5))
dev.off()

pdf(snakemake@output[["UH20_FOV11_Epi_T"]],width=2.5,height=2.5)
cluster_fov_plot(cosmx_proto %>% filter(fov == 11)%>% filter(x_global_px >= 43830,x_global_px <= 45830,y_global_px >= 120700),lwd=0.22,class='new_anno2',colors = c('#FFA8FF','#FF1AFF','cyan','blue','#A7B0BD'))+
theme(axis.text = element_blank(),axis.ticks = element_blank(),plot.title = element_blank(),
      legend.key.size = unit(0.2, 'cm'),
      legend.key.height = unit(0.2, 'cm'),
      legend.key.width = unit(0.2, 'cm'),
      legend.title = element_text(size=6),
      legend.text = element_text(size=5))
dev.off()
