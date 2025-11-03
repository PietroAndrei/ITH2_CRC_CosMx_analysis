#color palettes for categorical data
c24 <- c(
	 'dodgerblue2', '#E31A1C', # red
	 'green4',
	 '#6A3D9A', # purple
	 '#FF7F00', # orange
	 'gold1',
	 'skyblue2', '#FB9A99', # lt pink
	 'palegreen2',
	 '#CAB2D6', # lt purple
	 '#FDBF6F', # lt orange
	 'gray70', 'khaki2',
	 'maroon', 'orchid1', 'deeppink1', 'blue1', 'steelblue4',
	 'darkturquoise', 'green1', 'yellow4', 'yellow3',
	 'darkorange4', 'brown'
	 )
c12 <- c('#88CCEE','#CC6677','#DDCC77','#117733',
	 '#332288','#AA4499','#44AA99','#999933',
	 '#882255','#661100','#6699CC','#888888')


#Count normalization based on the total cell area instead of total counts
AreaNorm = function(seurat = NULL,log=T,assay='RNA'){
	if(is(seurat,'Seurat') == F){
		stop('Input data must be a Seurat object')
	}
	area = seurat$Area
	mat = LayerData(seurat,layer='counts')
	mat = (mat %*% Matrix::Diagonal(x=1/area))
	if(log == T){
		mat = log1p(mat)
	}
	colnames(mat) = colnames(seurat)
	LayerData(seurat,layer='data') = mat
	rm(mat)
	return(seurat)
}


#sfarrow function to binarize sfc/geometry columns data in a sf/data.frame object
#To go back to a sfc/geometry format, use sf::st_as_sfc function on the column(s) of interest
encode_wkb = function(df=NULL){
	  geom_cols = lapply(df, function(i) inherits(i, "sfc"))
  geom_cols = names(which(geom_cols==TRUE))

    df = as.data.frame(df)

    for(col in geom_cols){
	        obj_geo = sf::st_as_binary(df[[col]])
        attr(obj_geo, "class") = c("arrow_binary", "vctrs_vctr", attr(obj_geo, "class"), "list")
	    df[[col]] = obj_geo
	  }
      return(df)
}


##Convert an sf object to an arrow table
sf2arrow = function(df = NULL){
	df = encode_wkb(df)
	df = arrow::as_arrow_table(df)
	return(df)
}
#Revert sf2arrow transformation
arrow2sf = function(df=NULL){
	geom_cols = lapply(dplyr::collect(df),function(i) inherits(i,"arrow_binary"))
	geom_cols = names(which(geom_cols==TRUE))
	df = dplyr::collect(df)

	for(col in geom_cols){
		df[[col]] = sf::st_as_sfc(df[[col]])
	}
	attr(df,'class') = c('sf','data.frame')
	return(df)
}



#Load transcript metadata from Nanostring CosMx raw files
#Specify the raw data folder and the FOV IDs for which you want to collect transcripts information
#The default output would be an arrow table. If return_tibble is TRUE, return a tibble instead (larger object size)
#params path is a string corresponding to the main CosMx data folder ('/path/to/CosMx/data/')
#params slide is a string corresponding to the slide id in the form of S0/S1 (for Slide1/2)
#params fovs is a vector of numeric fov ids 
load_transcripts = function(path = NULL,slide=NULL,fovs=NULL,return_tibble=F){
	path = paste0(path,slide)
        folders = list.files(path)
        subf = folders[stringr::str_detect(folders,'Logs|polygons',negate=T)]
        full = paste0(path,'/',subf,'/AnalysisResults/')
	full_path = paste0(full,list.files(full))
	my_fovs = data.frame(id = fovs,digits=nchar(fovs))
	my_fovs = my_fovs %>% mutate(folder = ifelse(digits == 1,paste0('FOV00',id),ifelse(digits == 2,paste0('FOV0',id),paste0('FOV',id))))
	err_fovs = setdiff(my_fovs$folder,list.files(full_path))
	if(length(err_fovs) > 0){
		stop('You are lookig for FOVs not reported in the main path. Please check your FOVs list')
	}else{
		tot = lapply(my_fovs$folder,function(x){
			     files = list.files(paste0(full_path,'/',x))
			     files = files[stringr::str_detect(files,'complete_code_cell')]
			     transcr = read_csv_arrow(paste0(full_path,'/',x,'/',files),as_data_frame=F) %>%
			     select(CellId,fov,target,target_idx,codeclass,CellComp,seed_x,seed_y,x,y,z) %>%
			     as_arrow_table()
			     return(transcr)})
		tot = purrr::reduce(tot,concat_tables)
		if(return_tibble == T){
			return(collect(tot))
		}else{
			return(tot)
		}
	}
}

#Inital filtering of single cells based on total counts and genes detected, and negProbes/transcripts count ratio
filter_cosmx_cells=function(data = NULL,verbose=T,minCounts=20,minGenes=15,NegPosRatio=0.05){
	if(!is.data.frame(data) & is(data,'Seurat') ==F){
		stop('Please provide a data frame or a Seurat Object as an input file')
	}
	if(is(data,'Seurat')){
		newdata = data@meta.data
	}else{
		newdata = data
		rm(data)
	}
	all = nrow(newdata)
	if(verbose==T){
		message('Initial number of unfiltered cells: ',all)
	}
	newdata = newdata %>% filter(nCount_RNA >= minCounts)
	no_transcr = all - nrow(newdata)
	if(verbose==T){
		message('Removing ',no_transcr,' cells with less than ',minCounts,' transcripts')
	}
	newdata = newdata %>% filter(nFeature_RNA >= minGenes)
	no_genes = all - no_transcr - nrow(newdata)
        if(verbose==T){
		message('Removing ',no_genes,' cells with less than ',minGenes,' total expressed genes')
	}
	newdata = newdata %>% filter(nCount_negprobes/nCount_RNA < NegPosRatio)
	negprob = all - no_transcr - no_genes - nrow(newdata)
	if(verbose==T){
		message('Removing ',negprob,' cells with negative probes accounting for more than ',NegPosRatio*100,'% of total transcripts')
	}
	if(verbose==T){
		message(nrow(newdata),' cells passed the filtering. ',all-nrow(newdata),' cells have been discarded')
	}
	rm(all)
	rm(no_transcr)
	rm(no_genes)
	rm(negprob)
	if(is(data,'Seurat')){
		data = data[,WhichCells(data,cells=rownames(newdata))]
		rm(newdata)
		return(data)
	}else{
		return(newdata)
	}
}

#Convert single cell dataframe to sf data.frame object by crreating geometry based on single cell centroids X/Y coordinates
make_geometry_df = function(df,xcoord,ycoord){st_geometry(df) <- st_geometry(st_as_sf(df,coords = c(xcoord, ycoord))); return(df)}

#Expand centroids for a certain distance/radius with st_buffer(), until they come in contact with other centroids/polygons boundaries
#Function taken from https://github.com/r-spatial/sf/issues/824
#Useful shortcut when no real cell boundaries are available
st_buffer_without_overlap = function(centroids, dist){
	# Voronoi tesselation
	voronoi = centroids %>% 
	st_geometry() %>%
	st_union() %>%
	st_voronoi() %>%
	st_collection_extract()	         
	# Put them back in their original order
	voronoi =
		voronoi[unlist(st_intersects(centroids,voronoi))]
	# Keep the attributes
	result = centroids

	# Intersect voronoi zones with buffer zones
	st_geometry(result) =
		mapply(function(x,y) st_intersection(x,y),
		       st_buffer(st_geometry(centroids),dist), 
		       voronoi,
		       SIMPLIFY=FALSE) %>%
	st_sfc(crs=st_crs(centroids))

	result
}

##Identify,for each single cell, the whole list of touching cells. 
##This function will return the input data object updated with 2 or 3 new columns (depending on the 'label' parameters)
##params data should be a data frame and/or a sf object with a geometry 'POLYGON' column
##params precision indicate the precision value to apply for identifying touching cells. Default is 1
##params label string indicating data colname that will be used to identify the touching cells for each single cells.
##If label == NULL, 2 new columns will be added to the initial df: n_neighbors/neighbors (total number/cell IDs as data rownames)
##If lavel != NULL,a third column 'neighbors_label' will be added to data df with a string of cell IDs defined on the basis of the 
## 'label' column (e.g. label = NULL --> neighbors = '1,10,25'; label = 'cell_type' --> neighbors_label = 'CD8+,Macro,NK')
get_touching_cells = function(data=NULL,precision = 1,label = NULL){
	if(any(str_detect(colnames(data),'geometry')) ==F){
		      stop('Geometry not found. Please use a data.frame containing a sf column')
	}
	if(is.data.frame(data) & is(data,'sf') == F){
	      data = st_as_sf(data)
	}
	data = data %>% st_set_precision(.,precision)
	#use st_relate() instead of st_intersects() to identify groups of touching cells.
	#By specifying the right pattern in st_relate(), 
	#we can already filter out single cell self-calls (cell 'A' will be always found 'in contact' with itself)
	neighbors = st_relate(data,pattern = "****0****")
	data = data %>% mutate(n_neighbors = (neighbors %>% lapply(.,length) %>% unlist()),
			       neighbors = (lapply(neighbors,paste0,collapse=',') %>% unlist()))
	rm(neighbors)
	if(is.null(label) == F){
		neigh = strsplit(data$neighbors,',')
		lab = lapply(neigh,function(x){
			     newlab = as_tibble(data)[x,label] %>% unlist()
			     return(newlab)
			       })
		rm(neigh)
		lab = lapply(lab,paste0,collapse=',') %>% unlist()
		data$neighbors_label = lab
		rm(lab)
	}
	return(data)
}

##Variation of the previous function using st_intersects() instead of st_relate(). It should work better
get_touching_cells2= function(data=NULL,precision = 1,label = NULL){
	if(any(str_detect(colnames(data),'geometry')) ==F){
		stop('Geometry not found. Please use a data.frame containing a sf column')
	}
        if(is.data.frame(data) & is(data,'sf') == F){
		data = st_as_sf(data)
	}
	data = data %>% st_set_precision(.,precision)
	#use st_intersects() instead of st_relate() to identify groups of touching cells.
	neighbors = st_intersects(data)
	#we need to filter out single cell self-calls (cell 'A' will be always found 'in contact' with itself)
	neighbors = lapply(seq_len(length(neighbors)),function(x) neighbors[[x]][neighbors[[x]] != x])
	attr(neighbors,'class') = c('sgbp','list')
	data$n_neighbors = neighbors %>% lapply(.,length) %>% unlist()
	data$neighbors = lapply(neighbors,paste0,collapse=',') %>% unlist()
	rm(neighbors)
	if(is.null(label) == F){
		neigh = strsplit(data$neighbors,',')
		lab = lapply(neigh,function(x){
			     newlab = as_tibble(data)[x,label] %>% unlist()
			     return(newlab)
			      })
		rm(neigh)
		lab = lapply(lab,paste0,collapse=',') %>% unlist()
		data$neighbors_label = lab
		rm(lab)
	}
	return(data)
}


##Identify,for each single cell, the top k nearest cells
##params data, data frame and/or a sf object with a geometry column
##params id string, data colname that will be used to identify cell identifier for each single cell
##params label string, data colname that will be used to identify the knn cells for each single cell
##params knn numeric, number of nearest neighbors to search for each single cell
##return updated data with an additional 'knn_neighbors_label' column reporting the label identity of the top knn cells for each single cell.
get_nearest_cells = function(data, id = NULL, label=NULL, knn=5){
	if(any(str_detect(colnames(data),'geometry')) ==F){
		stop('Geometry not found. Please use a data.frame containing a sf column')
	}
	if(is.data.frame(data) & is(data,'sf') == F){
		data = st_as_sf(data)
	}
	require(spdep)
	require(spatialreg)
	#Find top k nearest neighbors for each cell
	knn_obj  = knearneigh(st_centroid(st_geometry(data)),k=knn)
	colnames(data)[colnames(data) == label] = 'label'
	colnames(data)[colnames(data) == id] = 'cell_ID' 
	#knn --> nb --> listw format
	nb = knn2nb(knn_obj,row.names=data$cell_ID)
	listw = nb2listw(nb,style='B',zero.policy=T)
	mat  = as(listw,'CsparseMatrix')
	# Convert matrix to long df format reporting only non-zero entries
	neigh = which(mat!=0,arr.ind=TRUE)  %>% as.data.frame()
	#Convert data into data.frame
	data= as.data.frame(data)
	#Add cell ids and neighboring cell types to neigh object
	neigh$cell_ID = data[neigh$row,]$cell_ID
	neigh$col.type = as.character(data[neigh$col,]$label)
	neigh = neigh %>% group_by(row) %>% 
	summarise(cell_ID = unique(cell_ID),neighbors_label = paste0(col.type,collapse=',')) %>% 
	ungroup() %>% 
	select(cell_ID,neighbors_label)
	colnames(neigh)[colnames(neigh) == 'neighbors_label'] = paste0('knn',knn,'_label')
	data = merge(data, neigh, by = 'cell_ID')
	colnames(data)[colnames(data) == 'label'] = label
	colnames(data)[colnames(data) == 'cell_ID'] = id 
	return(data)
}

##Create voronoi diagram around single cell centroids. Use convex hull to limit voronoi diagram extension to the actual tissue/slide size
make_sc_voronoi = function(data){
	if(any(str_detect(colnames(data),'geometry')) ==F){
		stop('Geometry not found. Please use a data.frame containing a sf column')
	}
	if(is.data.frame(data) & is(data,'sf') == F){
		data = st_as_sf(data)
	}
	voronoi = data %>% 
	st_geometry() %>% 
	st_union() %>% 
	st_voronoi() %>% 
	st_collection_extract()
	#Keep original data point order in voronoi geometry
	voronoi = voronoi[unlist(st_intersects(data,voronoi))]
	#Create convex hull polygon wrapping data points
	hull = st_convex_hull(st_combine(data))
	#Intersect Voronoi diagram and convex hull
	vor_hull = st_intersection(st_cast(voronoi), hull)
	return(vor_hull)
}



#From a metadata file, given a list of cells with the corrisponding annotated cell type, build a neighbor matrix reporting the number of occurrences for each cell type
#Among the top k nearest neighbors for each cell

KnnMat = function(meta=NULL,knn=100,label='final_anno'){
	if(any(str_detect(colnames(meta),'geometry')) ==F){
		stop('Geometry not found. Please use a data.frame containing a sf column')
	}
	if(is.data.frame(meta) & is(meta,'sf') == F){
		data = st_as_sf(meta)
	}
	require(spdep)
	require(spatialreg)
	require(dplyr)
	#Find top k nearest neighbors for each cell
	knn  = knearneigh(st_centroid(st_geometry(meta)),k=knn)
	colnames(meta)[colnames(meta) == label] = 'label'
	#knn --> nb --> listw format
	nb = knn2nb(knn,row.names=meta$cell_ID)
	listw = nb2listw(nb,style='B',zero.policy=T)
	mat  = as(listw,'CsparseMatrix')
	# Convert matrix to long df format reporting only non-zero entries
	neigh = which(mat!=0,arr.ind=TRUE)  %>% as.data.frame()
	#Convert meta data into data.frame
	meta= as.data.frame(meta)
	#Add cell ids and neighboring cell types to neigh object
	neigh$cell_ID = meta[neigh$row,]$cell_ID
	neigh$col.type = as.character(meta[neigh$col,]$label)
	#Count the occurrences for each cell type among the neighbors of each cell
	neigh_mat = neigh %>% group_by(cell_ID,col.type) %>% mutate(n = row_number()) %>% ungroup() %>% select(c(cell_ID,col.type,n))
	#From tabular to matrix format
	neigh_mat = reshape2::dcast(neigh_mat,cell_ID~col.type,value.var='n')
	rownames(neigh_mat) = neigh_mat$cell_ID
	neigh_mat$cell_ID = NULL
	#Be sure that all column in the final matrix are in numeric format
	neigh_mat = neigh_mat %>% mutate_all(as.numeric)
	#Convert any possible NA value to 0
	neigh_mat[is.na(neigh_mat)] = 0
	#Convert to sparse format
	neigh_mat = as(as.matrix(neigh_mat),'CsparseMatrix')
	return(neigh_mat)
}

#Like the KnnMat, but based on dnearneigh() function from spdep. The total number of neighbors will be different across cells
#NB maximum distance must be in line with the unit system of the geometry coordinates (um, pixel...) 
dnnMat = function(data,d.max=20,label='final_anno',ids='cell_ID',sparse=F){
	if(any(str_detect(colnames(data),'geometry')) ==F){
		stop('Geometry not found. Please use a data.frame containing a sf column')
	}
	if(is.data.frame(data) & is(data,'sf') == F){
		data = st_as_sf(data)
	}
	#For each cell, find all neighbors under a certain distance d.max
	dnn = dnearneigh(data,d1=0,d2=d.max)
	listw = nb2listw(dnn,style='B',zero.policy=T)
	mat  = as(listw,'CsparseMatrix')
	# Convert matrix to long df format reporting only non-zero entries
	neigh = which(mat!=0,arr.ind=TRUE)  %>% as.data.frame()
	#Convert meta data into data.frame
	data= as.data.frame(data)
	#Add cell ids and neighboring cell types to neigh object
	neigh$cell_ID = data[neigh$row,colnames(data)[colnames(data) == ids]]
	neigh$col.type = as.character(data[neigh$col,colnames(data)[colnames(data) == label]])
	#Count the occurrences for each cell type among the neighbors of each cell
	neigh_mat = neigh %>% 
	group_by(cell_ID,col.type) %>% 
	mutate(n = row_number()) %>% 
	ungroup() %>% 
	select(c(cell_ID,col.type,n))
	#From tabular to matrix format
	neigh_mat = reshape2::dcast(neigh_mat,cell_ID~col.type,value.var='n',fun.aggregate=length)
	rownames(neigh_mat) = neigh_mat$cell_ID
	neigh_mat$cell_ID = NULL
	#Be sure that all column in the final matrix are in numeric format
	neigh_mat = neigh_mat %>% mutate_all(as.numeric)
	#Convert any possible NA value to 0
	neigh_mat[is.na(neigh_mat)] = 0
	#Since some cells may not have neighbors in a distance <= d.max,they must be added to the final matrix
	if(nrow(neigh_mat)!= nrow(data)){
		missing = setdiff(data$cell_ID,rownames(neigh_mat))
		missing_df = data.frame(matrix(0,nrow=length(missing),ncol=ncol(neigh_mat)))
		colnames(missing_df) = colnames(neigh_mat)
		rownames(missing_df) = missing
		neigh_mat = rbind(neigh_mat,missing_df)
	}
	if(sparse==T){
		#Convert to sparse format
		neigh_mat = as(as.matrix(neigh_mat),'CsparseMatrix')
	}
	return(neigh_mat)
}
  



#Plot FOVs as single cell polygons. You can additionally plot cell centroids
# The main input must be a data frame with global coordinates for cell boundaries or 
# a sf object with a geometry column where the polygon objects have been already specified
fov_plot = function(poly = NULL,fovID = c(),
		    cell_col = 'steelblue1',lwd=1,
		    centroid=F,centroid_col='indianred4',centroid_size=1){
	if(is.data.frame(poly) ==F){
		stop('polygon object must be a data frame')
	}
	if(length(fovID > 0)){
		poly = poly %>% filter(fov %in% fovID)
	}else{
		stop('Please specify the field(s) of view you want to plot')
	}
	if(length(setdiff(fovID,poly$fov)) > 0){
		stop('You are lookig for FOVs not reported in your data frame. Please check your FOVs list')
	}
	if(is(poly,'sf') == F){
		message('Creating polygon geomtries from initial list of cell boundaries...')
		poly = poly %>% tidyr::unite('ID',fov:cellID)
		poly = df2sf(poly,spatialCoordsNames = c('x_global_px','y_global_px'),geometryType = 'POLYGON')
		message('...Done!')
	}
	g = ggplot()+geom_sf(data=poly,fill=cell_col,lwd=lwd)+
	theme(
	      panel.background = element_rect(fill = "black",
					      colour = "black",
					      size  = 0.5, 
					      linetype = "solid"),
	      panel.grid.major = element_line(linewidth = 0.5, 
					      linetype = 'solid',
					      colour = "black"), 
	      panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
					      colour = "black")
	      )
	if(centroid == T){
		g = g + geom_sf(data=st_centroid(poly$geometry),size=1,color=centroid_col)
	}
	g
}


##Plot single cells across FOVs with colour gradient based on a numeric feature
#data must be a data.frame or a sf object with global single cell coordinates and the quantitative variable of interest
#max.feat is the percentile value that will be used as a limit for the color scale
fov_heatmap = function(data = NULL,
		       feature = NULL,
		       max.feat = 0.95,
		       min.feat = 0.05,
		       polygon =FALSE,lwd=0,
		       fovs = c(),
		       pt.size=0.3,
		       pt.shape=16,
		       cols = c(
			'midnightblue',
			'#1b129a',
			'#3429C4',
			'#FF6060',
			'#EEAA42'#'#FFFEA8'
			)
		       ){
if(is(data,'Seurat')){
	data = FetchData(data,c(feature,colnames(data@meta.data)))
  }
if(length(fovs)>0){
	data = data %>% filter(fov %in% fovs)
}
colnames(data)[colnames(data) == feature] = 'feat'
max.val = quantile(data$feat,max.feat)
min.val = quantile(data$feat,min.feat)
if(polygon==T){
	if(is(data,'sf') ==F){
		data = sf::st_as_sf(data)
	}
	g = data %>% ggplot()+geom_sf(aes(fill=feat),lwd=lwd)+
	scale_fill_gradientn(limits=c(min.val,max.val),
			     breaks = c(min.val,max.val),
			     labels = c('Min','Max'),
			     na.value=cols[length(cols)],colors=cols)+
	labs(fill = feature)
}else{
g = data %>% ggplot(.,aes(x=x_global_px,y=y_global_px))+
    theme_classic()
if(pt.shape %in% c(21:25)){
   g = g + geom_point(aes(fill=feat),shape=pt.shape,size=pt.size) + 
   scale_fill_gradientn(limits=c(0,max.val),na.value=cols[length(cols)],colors=cols)+
   labs(fill = feature)
}else{
	g = g+geom_point(aes(color=feat),shape=pt.shape,size=pt.size)+
	scale_color_gradientn(limits=c(0,max.val),na.value=cols[length(cols)],colors=cols)+
	labs(color=feature)
}
}
g +theme(
	 panel.background = element_rect(fill = "black",
					 colour = "black",
					 size  = 0.5, 
					 linetype = "solid"),
	 panel.grid.major = element_line(linewidth = 0.5, 
					 linetype = "solid",
					 colour = "black"), 
	 panel.grid.minor = element_line(linewidth = 0.25, linetype = "solid",
					 colour = "black"),
	 plot.background = element_rect(fill = 'black',
					colour = 'black'),
	 legend.background = element_rect(fill='black'),
	 legend.text = element_text(colour='white'),
	 legend.title=element_text(colour='white'),
	 legend.position = 'top'
	 )
}


#Plot single cells across FOVs with a color scheme based on a categorical variable
cluster_fov_plot = function(meta=NULL,class = NULL,colors=c12,pt.size=0.05,polygon=T,lwd=0.1){
	sample = unique(meta$sample)
	meta = meta %>% rename(.,'class'=class)
	if(polygon ==T){
		g = meta %>% sf::st_as_sf() %>%
		ggplot()+geom_sf(aes(fill=class),lwd=lwd)+
		scale_fill_manual(values=colors)+
		labs(fill='Cell Type')
        }else{
		g = meta %>% select(x_global_px,y_global_px,class) %>%
		ggplot(.,aes(x=x_global_px,y=y_global_px))+theme_classic()+
		geom_point(aes(color=class),size=pt.size) +
		scale_color_manual(values=colors)+
		labs(color='Cell Type')
	}
	g + theme(
	      panel.background = element_rect(fill = "black",
					      colour = "black",
					      linewidth  = 0.5, linetype = "solid"),
	      panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
					      colour = "black"), 
	      panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
					      colour = "black")
	      )+
	ggtitle(sample)+theme(plot.title = element_text(color='black',hjust=0.5))+
	guides(color = guide_legend(override.aes = list(size = 7.5)))+
	theme(legend.title = element_text(size=12,face='bold',color='black'),
	      legend.text = element_text(face='bold',color='black'))
}




cosmx_featureplot = function(data=NULL,feature = NULL,
			     label=NULL,grouping = NULL,
			     geom='boxplot',geom_args=list(),cols=c12,
			     my_comps = NULL){
	if(!is.data.frame(data) & !is(data,'Seurat')){
		stop('data must be a data frame or a Seurat Object')
	}
	if(is(data,'Seurat')){
		columns = c(feature,label,grouping)
		data = data %>% FetchData(.,columns)
	}
	colnames(data)[colnames(data) == feature] = 'gene'
	colnames(data)[colnames(data) == label] = 'label'
	if(!is.null(grouping)){
		colnames(data)[colnames(data) == grouping] = 'group'
	}
	if(!is.null(grouping)){
		g = ggplot(data,aes(x=group,y=gene,fill=label))+
		theme_classic()+
		do.call(what=paste0('geom_',geom),
			args=geom_args
			)+
		scale_fill_manual(values=cols)+
		facet_wrap(vars(group))
	}else{
		g = ggplot(data,aes(x=label,y=gene,fill=label))+
		theme_classic()+
		do.call(what=paste0('geom_',geom),
			args=geom_args
			)+
		scale_fill_manual(values=cols)
	}
	g+
	#stat_compare_means(comparisons = my_comps)+ install ggpubr for this 
	stat_summary(fun='median',fill='white',geom='label',aes(label=round(after_stat(y),2)))+
	ylab(feature)+xlab(label)+labs(fill=label)+
	theme(axis.text.x = element_text(size=10,color='black'))
}



###Volcano plot to visualize differential expression results. Works with Seurat DE output format (colnames:p_val/avg_log2FC/pct.1/pct.2/p_val_adj)
#'
#' @param seurat_de DE results in a Seurat format (data.frame,colnames:p_val/avg_log2FC/pct.1/pct.2/p_val_adj)
#' @param clust character/numeric cluster identifier, when DE results are reported for multiple clusters/groups (e.g. FindAllMarkers() output)
#' @param padj_threshold numeric, minimum adjusted pvalue required for a gene to be considered differentially expressed
#' @param min.pct numeric, minimum % expressing cells in either one of the two clusters compared (or both) required for a gene to be considered differentially expressed
#' @param logfc_thresh_up numeric, minimum log2(FC) value required for a gene to be considered significantly UP-regulated
#' @param logfc_thresh_down numeric, maximum log2(FC) value required for a gene to be considered significantly DOWN-regulated
#' @param label.default boolean, whether or not using the default gene labeling settings (i.e., top 25 UP + 25 DOWN genes)
#' @param label.genes character vector, gene names to highlight in the plot when label.default=FALSE
#' @param label.size numeric, gene name labels fontsize
#' @param pt.size numeric, size of the volcano dots
#' @param col_up character string, color that will be used to highlight UP-regulated genes
#' @param col_down character string, color that will be used to highlight DOWN-regulated genes
#' @param outType character string, one between 'seurat' or 'deseq'. If outType = 'deseq', absolute value of Wald statistics will be plotted on the y axis instead of FDR
#' @return ggplot2 volcano plot reporting all genes tested for differential expression with a colorscheme based on the FDR and Fold Change
#' @examples
#'de_volcano_new(de_res,logfc_thresh_up = 0.5,logfc_thresh_down = -0.5,outType = 'seurat',label.default = T,pt.size=0.5,label.size=1.5,col_up = 'cyan',col_down = 'blue')
de_volcano = function(seurat_de = NULL,clust=0,
		      padj_thresh = 0.01,min.pct=0,
		      logfc_thresh_up = log2(1.1),logfc_thresh_down = log2(0.9),
		      label.default=T,label.genes = c(),label.size=3,pt.size=2,
		      col_up = 'goldenrod2',col_down = 'steelblue3',
		      outType = c('seurat','deseq')){
	if(length(intersect(colnames(seurat_de),"cluster")) > 0){
		seurat_de = seurat_de %>% filter(cluster == clust)
	}
	if(length(intersect(colnames(seurat_de),"gene")) == 0){
		seurat_de$gene = rownames(seurat_de)
	}
	if(length(colnames(seurat_de)[str_detect(colnames(seurat_de),'pct.1|pct.2')]) <= 1){
		seurat_de$pct.1 = 1
		seurat_de$pct.2 = 1
	}
	colnames(seurat_de)[str_detect(colnames(seurat_de),'avg_log2FC|log2FoldChange')] = 'avg_log2FC'
	colnames(seurat_de)[str_detect(colnames(seurat_de),'padj|p_val_adj')] = 'p_val_adj'
	seurat_de$FC = 2^seurat_de$avg_log2FC
	seurat_de = seurat_de %>% arrange(desc(FC)) 
	seurat_de$label = ''
	seurat_de= seurat_de %>% filter(!is.na(p_val_adj))
	if(label.default ==T){
		seurat_de[1:25,]$label = seurat_de[1:25,]$gene
		n = nrow(seurat_de)
		seurat_de[(n-24):n,]$label = seurat_de[(n-24):n,]$gene
		seurat_de[(seurat_de$avg_log2FC < logfc_thresh_up & seurat_de$avg_log2FC > logfc_thresh_down) | 
			  (seurat_de$p_val_adj > padj_thresh) | (seurat_de$pct.1 <= min.pct & seurat_de$pct.2 <= min.pct),]$label = ''
	}else if(label.default ==F & length(label.genes) > 0){
		seurat_de[seurat_de$gene %in% label.genes,]$label = seurat_de[seurat_de$gene %in% label.genes,]$gene
		seurat_de[(seurat_de$avg_log2FC < logfc_thresh_up & seurat_de$avg_log2FC > logfc_thresh_down) | 
			  (seurat_de$p_val_adj > padj_thresh)| (seurat_de$pct.1 <= min.pct & seurat_de$pct.2 <= min.pct),]$label = ''
	}
	seurat_de$color = ifelse(
				 seurat_de$avg_log2FC > logfc_thresh_up & seurat_de$p_val_adj < padj_thresh & (seurat_de$pct.1 > min.pct | seurat_de$pct.2 > min.pct),'UP',ifelse(
																						  seurat_de$avg_log2FC < logfc_thresh_down & seurat_de$p_val_adj < padj_thresh & (seurat_de$pct.1 > min.pct | seurat_de$pct.2 > min.pct),'DOWN','N.S.')
				 )
	seurat_de$log10padj = -log10(seurat_de$p_val_adj)
	max_padj = max(seurat_de$log10padj)
	if(max_padj == Inf){
		max_padj_new = max(seurat_de[seurat_de$log10padj != Inf,]$log10padj)+10
		seurat_de[seurat_de$log10padj == Inf,]$log10padj = round(max_padj_new)
	}else{
		max_padj_new= max_padj
	}
	if(is.null(outType)){
		outType = 'seurat'
	}
	if(outType == 'seurat'){
		g = ggplot(seurat_de,aes(x=FC,y=log10padj,label=label,colour=color))+
		theme_classic()+
		geom_point(stat='identity',size=pt.size)+
		scale_color_manual(values=c('UP' = col_up,
					    'N.S.' = 'gray90',
					    'DOWN' = col_down))+
		geom_hline(yintercept = -log10(padj_thresh),linetype='dashed',color ='black')+
		geom_vline(xintercept = 2^(logfc_thresh_up),linetype='dashed',color ='black')+
		geom_vline(xintercept = 2^(logfc_thresh_down),linetype='dashed',color ='black')+
		ylab('Adjusted P-value')+labs(fill = 'Sign')+xlab('FC')+
		geom_text_repel(max.overlaps = Inf,seed=1,colour='black',force=2,box.padding=0.5,size=label.size,segment.size=0.1)+
		coord_trans(x = "log2")+
		theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
		#Keep log scale but shows real padj on the y axis
		scale_y_continuous(breaks = c(round(max_padj_new/seq_len(4)),2),
				   labels=function(x) ifelse(x > 0 & x < max_padj_new ,label_parsed(paste("10^",-1*x)),
							     ifelse(x >= max_padj_new & max_padj == Inf,'0',
								    ifelse(x >= max_padj_new & max_padj != Inf,label_parsed(paste("10^",-1*x)),'1')
								    )
							     )
				   )+
		scale_x_continuous(breaks=c(round(max(seurat_de$FC)),3,2,1.5,1.25,1,0.75,0.5,0.25,0.1))+
		theme(axis.ticks = element_line(linewidth =0.5,colour='gray40'),axis.ticks.length=unit(3,'pt'))
	}
	if(outType == 'deseq'){
		max_stat = max(abs(seurat_de$stat))
		min_stat = min(abs((seurat_de %>% filter(abs(avg_log2FC) > 1,p_val_adj < 0.05))$stat))
		g = ggplot(seurat_de,aes(x=FC,y=abs(stat),label=label,color=color))+
		theme_classic()+
		#annotate('rect',xmin=2^(logfc_thresh_up),xmax=Inf,ymin=-log10(padj_thresh),ymax =Inf,alpha=0.2,fill=col_up)+
		#annotate('rect',xmin=-Inf,xmax=2^(logfc_thresh_down),ymin=-log10(padj_thresh),ymax =Inf,alpha=0.2,fill=col_down)+
		geom_point(size=pt.size)+
		scale_color_manual(values=c('UP' = col_up,
					    'N.S.' = 'gray90',
					    'DOWN' = col_down)
		)+
		geom_hline(yintercept = min_stat,linetype='dashed',color ='black')+
		geom_vline(xintercept = 2^(logfc_thresh_up),linetype='dashed',color ='black')+
		geom_vline(xintercept = 2^(logfc_thresh_down),linetype='dashed',color ='black')+
		ylab('|Wald Statistic|')+labs(fill = 'Sign')+xlab('FC')+
		geom_text_repel(max.overlaps = Inf,seed=1,colour='black',size=label.size,segment.size=0.1)+
		coord_trans(x = "log10")+
		theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
		#Keep log scale but shows real padj on the y axis
		scale_y_continuous(breaks = c(round(max_stat/seq_len(4)),round(min_stat,1)))+
		scale_x_continuous(breaks=c(round(max(seurat_de$FC)/c(1,10,100,1000)),1,0.5,0.25,0.1))+
		theme(axis.ticks = element_line(linewidth =0.5,colour='gray40'),axis.ticks.length=unit(3,'pt'))
	}
	g
}



##Check https://github.com/satijalab/seurat/issues/3712
# We want a different behavior from MASTDETest, where the latent variable is
# used as a random effect in the model
MASTDE_RandTest <- function(
			   data.use,
			   cells.1,
			   cells.2,
			   latent.vars = NULL,
			   verbose = TRUE,
			   # New option - random effect variable (should be included in latent.vars)
			   re.var = NULL,
			   ...
			   ) {
	# Check for MAST
	if (!PackageCheck('MAST', error = FALSE)) {
		stop("Please install MAST - learn more at https://github.com/RGLab/MAST")
	}
	group.info <- data.frame(row.names = c(cells.1, cells.2))
	latent.vars <- latent.vars %||% group.info
	group.info[cells.1, "group"] <- "Group1"
	group.info[cells.2, "group"] <- "Group2"
	group.info[, "group"] <- factor(x = group.info[, "group"])
	latent.vars.names <- c("condition", colnames(x = latent.vars))
        latent.vars <- cbind(latent.vars, group.info)
	latent.vars$wellKey <- rownames(x = latent.vars)
	fdat <- data.frame(rownames(x = data.use))
	colnames(x = fdat)[1] <- "primerid"
	rownames(x = fdat) <- fdat[, 1]
	sca <- MAST::FromMatrix(
				exprsArray = as.matrix(x = data.use),
				check_sanity = FALSE,
				cData = latent.vars,
				fData = fdat
				)
	cond <- factor(x = SummarizedExperiment::colData(sca)$group)
	cond <- relevel(x = cond, ref = "Group1")
	SummarizedExperiment::colData(sca)$condition <- cond
	# This is the main change in the code - we want ~ ( 1 | re.var) in the formula:
	# ~ condition + lat.vars + (1 | re.var)
	if (!is.null(re.var)) {
		if (!re.var %in% latent.vars.names) {
			stop("Random effect variable should be included in latent variables!")
		}
		latent.vars.names <- latent.vars.names[!latent.vars.names %in% re.var]
		fmla <- as.formula(
				   object = paste0(
						   " ~ ", paste(latent.vars.names, collapse = "+"), glue::glue("+ (1|{re.var})")
						   )
				   )
		# print(fmla)
		# We need glmer to make this work
		method <-  "glmer" # trying to troubleshoot this - it can clash with the already existing method var in the function call    
		zlmCond <- MAST::zlm(formula = fmla, sca = sca, method = "glmer", ...)
	} else {
		# Original code
		fmla <- as.formula(
				   object = paste0(" ~ ", paste(latent.vars.names, collapse = "+"))
				   )
		zlmCond <- MAST::zlm(formula = fmla, sca = sca, ...)
	}
	summaryCond <- MAST::summary(object = zlmCond, doLRT = 'conditionGroup2')
	summaryDt <- summaryCond$datatable	    
	# The output format is slightly different, so we need adapt the code
	if(!is.null(re.var)) {
		p_val <- summaryDt[summaryDt$"component" == "H", 4]$`Pr(>Chisq)`
		genes.return <- summaryDt[summaryDt$"component" == "H", 1]$primerid
	} else {
		p_val <- summaryDt[summaryDt[, "component"] == "H", 4]
		genes.return <- summaryDt[summaryDt[, "component"] == "H", 1]
	}	      
	to.return <- data.frame(p_val, row.names = genes.return)
	return(to.return)
}

## We will replace the original function with ours inside the Seurat namespace
#assignInNamespace('MASTDETest', MASTDE_RandTest, asNamespace("Seurat"))



##Adjust any vector of gene names (e.g. MSigDB gene sets) to Cosmx scheme
#Input x must be a vector of gene names
hallmark2cosmx = function(x,collapse=T){
	x[x %in% c('CCL3','CCL3L1','CCL3L3')] = 'CCL3/L1/L3'
	x[x %in% c('CCL4','CCL4L1','CCL4L2')] = 'CCL4/L1/L2'
	x[x %in% c('CXCL1','CXCL2','CXCL3')] = 'CXCL1/2/3'
	x[x %in% c('EIF5A','EIF5AL1')] = 'EIF5A/L1'
	x[x %in% c('FCGR3A','FCGR3B')] = 'FCGR3A/B'
	x[x %in% c('HBA1','HBA2')] = 'HBA1/2'
	x[x %in% c('HCAR2','HCAR3')] = 'HCAR2/3'
	x[x %in% c('HLA-DQB1','HLA-DQB2')] = 'HLA-DQB1/2'
	x[x %in% c('HSPA1A','HSPA1B')] = 'HSPA1A/B'
	x[x %in% c('IFNA1','IFNA13')] = 'IFNA1/13'
	x[x %in% c('IFNL2','IFNL3')] = 'IFNL2/3'
	x[x %in% c('KRT6A','KRT6B','KRT6C')] = 'KRT6A/B/C'
	x[x %in% c('MAP1LC3B','MAP1LC3B2')]  = 'MAP1LC3B/2'
	x[x %in% c('MZT2A','MZT2B')] = 'MZT2A/B'
	x[x %in% c('PF4','PF4V1')] = 'PF4/V1'
	x[x %in% c('SAA1','SAA2')] = 'SAA1/2'
	x[x %in% c('TNXA','TNXB')] = 'TNXA/B'
	x[x %in% c('TPSAB1','TPSB2')] = 'TPSAB1/B2'
	x[x %in% c('XCL1','XCL2')] = 'XCL1/2'
	x[x %in% c('HLA-A','HLA-B','HLA-C')] = 'MHC I'
	if(collapse == T){
		x = unique(x)
	}
	return(x)
}

#Adapt CellChatDB.human interaction dataframe (1st object of CellChat.human list) to work with CosMx 1000-plex panel
cellchat2cosmx = function(cc_df = NULL){
	res = ddply(cc_df,.(interaction_name),function(x){
		    receptors = unique(unlist(strsplit(x$receptor.symbol,', ')))
		    ligands = unique(unlist(strsplit(x$ligand.symbol,', ')))
		    all_combs = expand.grid(ligands,receptors)
		    colnames(all_combs) = c('ligand.symbol','receptor.symbol')
		    return(all_combs)})

	res$receptor.symbol = hallmark2cosmx(res$receptor.symbol,collapse=F)
	res$ligand.symbol = hallmark2cosmx(res$ligand.symbol,collapse=F)
	return(res)
}

#Adapt (sparse) gene expression matrix to work with CosMx 1000-plex panel names
gex2cosmx = function(mat){
	### Define genes that will be collapsed in the gene expression matrix 
	Genes = list(
		     'CCL3/L1/L3' = c('CCL3','CCL3L1','CCL3L3'),
		     'CCL4/L1/L2' = c('CCL4','CCL4L1','CCL4L2'),
		     'CXCL1/2/3' = c('CXCL1','CXCL2','CXCL3'),
		     'EIF5A/L1' = c('EIF5A','EIF5AL1'),
		     'FCGR3A/B' = c('FCGR3A','FCGR3B'),
		     'HBA1/2' = c('HBA1','HBA2'),
		     'HCAR2/3' = c('HCAR2','HCAR3'),
		     'HLA-DQB1/2' = c('HLA-DQB1','HLA-DQB2'),
		     'HSPA1A/B' = c('HSPA1A','HSPA1B'),
		     'IFNA1/13' = c('IFNA1','IFNA13'),
		     'IFNL2/3' = c('IFNL2','IFNL3'),
		     'KRT6A/B/C' = c('KRT6A','KRT6B','KRT6C'),
		     'MAP1LC3B/2' = c('MAP1LC3B','MAP1LC3B2'),
		     'MZT2A/B' = c('MZT2A','MZT2B'),
		     'PF4/V1' = c('PF4','PF4V1'),
		     'SAA1/2' = c('SAA1','SAA2'),
		     'TNXA/B' = c('TNXA','TNXB'),
		     'TPSAB1/B2' = c('TPSAB1','TPSB2'),
		     'XCL1/2' = c('XCL1','XCL2'),
		     'MHC I' = c('HLA-A','HLA-B','HLA-C')
		     )
	require(Matrix)
	### Create a new matrix with collapsed gene expression counts
	newlines = lapply(names(Genes),function(x){
			  genes = intersect(Genes[[x]],rownames(mat))
			  if(length(genes) > 1){
			  	gene = colSums(mat[genes,])
			  	gene = as(gene,'sparseVector')
				gene = t(as(gene,'sparseMatrix'))
			  }else if(length(genes) == 1){
				gene = mat[genes,]
			  	gene = as(gene,'sparseVector')
			  	gene = t(as(gene,'sparseMatrix'))
			  }else{
				gene = NULL
			  }
			  return(gene)
		     })
	names(newlines) = names(Genes)
	#If present, remove empty vectors
	newlines = newlines %>% purrr::compact()
	### Add the collapsed matrix to the original one, and remove rows corresponding to those genes which have been now merged
	newmat = purrr::reduce(newlines,rbind)
	rownames(newmat) = names(newlines)
	mat=rbind(mat, newmat)
	mat=mat[setdiff(rownames(mat),unique(unlist(Genes))),]
	return(mat)
}



#Compute Spatial(Moran) Eigenvectors from single cell connectivity matrix (binary matrix reporting cell-cell contacts)
#The resulting eigenvectors can be used as covariates to remove spatial correlation effect from downstream analysis
#' @param data must be either a sf dataframe or a binary weight matrix
#' @param n_eigen specify the number of eigenvectors/values to extract from the centered weight matrix
#' @param nb_type character stringr specifying how to define neighboring cells. It can be one of the following:
#' 'touching','voronoi','knn' (knn still need to be implemented)
#' @param k numeric, used when nb_type is 'knn'. K number of nearest neighbors to return for each cell (NB* based on centroids)
cosmx_meigen = function(data = NULL,nb_type = NULL,k = 10,n_eigen=200){#l_chunks = (nrow(data)/1000)){
	if(!is(data,'sf') & !is(data,'CsparseMatrix') & !is(data,'matrix')){
		stop('Data must be an sf object or a sparse/dense weight matrix')
	}
	if(is.null(nb_type)){
		nb_type = 'touching'
	}
	if(is(data,'sf')){
		require(sf)
		data = st_set_precision(data,1)
		if(nb_type == 'touching'){#Define cells as neighbors if they are physically in contact
			nb = st_intersects(data)
			#we need to filter out single cell self-calls (cell 'A' will be always found 'in contact' with itself)
			nb = lapply(seq_len(length(nb)),function(x) nb[[x]][nb[[x]] != x])
			attr(nb,'class') = c('sgbp','list')
		}
		if(nb_type == 'voronoi'){#Define neighbors as touching polygons resulting from Voronoi tessellation
			# fov-wise voronoi
			datas = split(data,f=data$fov)
			vors = lapply(datas,function(x){
				cents = st_centroid(st_geometry(x))
				vor =st_voronoi(do.call(c,cents))
				vor = st_collection_extract(st_sfc(vor))
				x$cents = cents
				st_geometry(x) = cents
				vor = vor[unlist(st_intersects(x,vor))]
				return(vor)
				   }
			)
			#merge back voronoi-derived polygons and find neighbors (with st_touches())
			datas = lapply(1:length(datas),function(x){datas[[x]]$voronoi = vors[[x]]; return(datas[[x]])})
			datas = purrr::reduce(datas,rbind) %>% arrange(cell_ID)
			st_geometry(datas) = 'voronoi'
			nb  = st_touches(datas)
			rm(datas)
			rm(vors)
		}
		if(nb_type == 'knn'){#use spdep knn implementation to define neighbors
			kns = spdep::knearneigh(st_centroid(data$geometry),k=k)
			nb = spdep::knn2nb(kns)
			rm(kns)
		}
		#Convert neighbors list to a nb object
		if(nb_type != 'knn'){
		#st_as_nb
		attrs = attributes(nb)
		nb = lapply(nb, function(i) { if(length(i) == 0L) 0L else i } )
		attributes(nb) = attrs
		class(nb) = "nb"
		}
		listw = spdep::nb2listw(nb,zero.policy=T,style='B')
		require(spatialreg)
		data = as(listw,'CsparseMatrix')
		rm(listw)
		rm(nb)
	}else{
		require(spatialreg)
		#if data input is a weight matrix, verify its binary format
		if(length(data@x[!(data@x %in% c(0,1))] > 0)){
			data@x[data@x != 0] = 1
			data = as(data,'CsparseMatrix')
		}
	}
	require(BPCells)
	#Compress matrices with BPCells, so you don't need to chunk them (and you also save memory)
	message('Compressing connectivity matrix with BPCells.')
	tmpdir_C = tempfile('mat')
	tmpdir_M = tempfile('mat')
	message(paste0('Writing matrix in ',tmpdir_C,' directory...'))
	M = Matrix::Diagonal(nrow(data))
	data = data %>% write_matrix_dir(.,dir=tmpdir_C)
	#Convert Diagonal matrix to a CSC format for BPCells compatibility
	M = as(M,'dgCMatrix') %>% write_matrix_dir(.,dir=tmpdir_M)
	#Prepare centering matrix
	M = M - 1/nrow(data)
	#MCM transformation to center the (symmetric) connectivity  matrix
	data = M %*% data %*% M
	#Use RSpectra function interface for BPCells compatibility
	f = function(x,args){args %*% x}
	if(nb_type != 'knn'){#KNN-based connectivity matrix is not symmetric
		eigs = RSpectra::eigs_sym(f,k=n_eigen,n=nrow(data),args=data,which='LA')
	}else{
		eigs = RSpectra::eigs(f,k=n_eigen,n=nrow(data),args=data,which='LM')
	}
	message(paste0('Removing ',tmpdir_C,' and ',tmpdir_M))
	system(paste0('rm -r ',tmpdir_C))
	system(paste0('rm -r ',tmpdir_M))
	return(eigs)
}

#'Calculate the shortest distance between cells belonging to two different groups 
#'(i.e., find the closest cell of type B for each cell of type A), then save the corresponding euclidean distance in a new column
#' If plot is TRUE, plot a spatial heatmap of cells (type A/B) with the colorscale based on the euclidean distance of cells A to the closest cell B
#'
#' @param data must be a data frame (ideally an sf data frame)
#' @param ct.1 cell type A label
#' @param ct.2 cell type B label
#' @param ct1.lwd linewidth to use in plot for cell type A
#' @param ct2.lwd linewidth to use in plot for cell type B
#' @param class column name where the cell type labels are stored
#' @param plot boolean, wheter or not include a spatial splot in addition to the main CT1-->CT2 distance calculation
#' @param fov list fov IDs to include in the analysis. If NULL, no fov-based filter will be applied to the input data
#' @param ct1.colthresh quantile threshold to apply to the distance colorscheme for cells belonging to cell type A in the spatial plot 
#' (e.g., ct1.colthresh = 0.99 --> any distance greater than the 99% percentile will not result in a more 'intense' (brighter/darker) cell color
#' @param ct2.col discrete color to use for plotting cells belonging to cell type A
#' @return if plot is FALSE, a sf data frame reporting all the original metadata plus the euclidean distance between each cell 'A' and its closest cell 'B'. 
#' If plot is TRUE, a list object containing the result data frame plus a spatial plot reporting cells 'A' colored according to their distance from cells 'B'
get_cosmx_ct_dists = function(data=NULL,ct.1=NULL,ct1.lwd=0,ct.2=NULL,ct2.lwd=0,class=NULL,plot=T,fovs=NULL,ct1.colthresh = 0.99,ct2.col = 'white'){
	if(length(fovs)>0){
		data = data %>% filter(fov %in% fovs)
	}
	colnames(data)[colnames(data) == class] = 'ct_var'
	ct1 = data %>% filter(ct_var %in% ct.1)
	ct2 = data %>% filter(ct_var %in% ct.2)
	nearests = st_nearest_feature(ct1$geometry,ct2$geometry)
	nearest_dist = st_distance(x=ct1$geometry,y=ct2[nearests,]$geometry,by_element = T)
	ct1$dist = nearest_dist
        ct2$dist = 0
	if(is(ct1,'sf') == F){
		ct1 = st_as_sf(ct1)
	}
	if(is(ct2,'sf') == F){
		ct2 = st_as_sf(ct2)
	}
	if(plot == T){
		g = ct2 %>% ggplot()+
		geom_sf(aes(fill =ct_var),lwd=ct2.lwd)+
		scale_fill_manual(values= ct2.col)+
		theme(legend.key = element_rect(fill = "white", colour = "black",linewidth=1))+
		labs(fill= '')+
		theme(legend.position = 'top')+
		guides(fill = guide_legend(title.position = "top"))+
		ggnewscale::new_scale_fill()+
		geom_sf(data= ct1,aes(fill=dist),lwd=ct1.lwd)+
		scale_fill_gradientn(limits=c(0,quantile(ct1$dist,ct1.colthresh)),breaks=c(0,quantile(ct1$dist,ct1.colthresh)),labels=c('0','Max'),
				     na.value = 'midnightblue',
				     colors=c('#EEAA42','#FF6060','#3429C4','#1b129a','midnightblue')
				     )+
		labs(fill=glue::glue('{ct.1} --> {ct.2} Dist'))+
            	theme(
			  panel.background = element_rect(fill = "black",
							  colour = "black",
							  size  = 0.5, 
							  linetype = "solid"),
		          panel.grid.major = element_line(linewidth = 0.5, 
							  linetype = 'solid',
							  colour = "black"), 
		          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
							  colour = "black"),
		          plot.background = element_rect(fill = 'black',
							 colour = 'black'),
		          legend.background = element_rect(fill='black'),
			  legend.text = element_text(colour='white'),
			  legend.title=element_text(colour='white'),
			  legend.position = 'top'
			  )+
		guides(fill = guide_colorbar(title.position = "top"))
		ct.2 = gsub("\\+","",ct.2)
          	colnames(ct1)[colnames(ct1) == 'dist'] = glue::glue('{ct.2}_dist')
	      	return(list(dist = ct1,plot=g))
	}else{
		colnames(ct1)[colnames(ct1) == 'dist'] = glue::glue('{ct.2}_dist')
	        return(ct1)
	}
}

