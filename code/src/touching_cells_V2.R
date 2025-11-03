library(sf)
library(dplyr)
library(stringr)
utils = snakemake@params[["utils"]]
source(utils)

#Try to define touching cells with st_intersects instead of st_relate (it should be more precise...)
all_celltypes = c('Epi','FibroEndoMuscle','Myeloid','T/NK','T_Other','Plasma/B','Mast')
meta = readRDS(snakemake@input[["meta"]])
touching_cells = get_touching_cells2(meta,label='final_anno')
touching_cells$is_epiT = 0
touching_cells$nb_tnk = str_split(touching_cells$neighbors_label,',') %>% str_count(.,'T\\/NK')
touching_cells = touching_cells %>% mutate(is_epiT = ifelse(final_anno == 'Epi' & nb_tnk > 0,1,0))
touching_cells$is_epionly = 0
touching_cells$nb_epi = str_split(touching_cells$neighbors_label,',') %>% str_count(.,'Epi')
touching_cells = touching_cells %>% mutate(is_epionly = ifelse(nb_epi == n_neighbors & n_neighbors > 0 & final_anno == 'Epi',1,0))
touching_cells = touching_cells %>% mutate(epiType = ifelse(
							    is_epiT == 1,'epiT',
							    ifelse(
								   is_epionly == 1,'epiOnly',
								   ifelse(
									  final_anno != 'Epi','StromaImmune','epiOther'
									  )
								   )
							    )
)
touching_cells$nb_noepi = str_count(touching_cells$neighbors_label,paste0(all_celltypes[all_celltypes != 'Epi'],collapse='|'))

final_meta = snakemake@output[["final_meta"]]
saveRDS(touching_cells,final_meta)
