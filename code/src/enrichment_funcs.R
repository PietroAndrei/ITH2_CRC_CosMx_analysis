
#Wrapper function for ClusterProfiler::enricher() function. Assume to receive a data frame with gene symbols as rownames or a vector of gene symbols
#paths should be a data frame consiting of at least two columns: mapped pathway name and entrez id.
#A third paths column is needed if you want to use which.ids = 'symbol'
wrap_enrich = function(genes=NULL,which.ids = c('entrez','symbol'),paths=NULL){
	if(is.data.frame(genes)){
		message("Considering ",nrow(genes)," DE genes")
		genes = rownames(genes)
	}else if(is.vector(genes)){
		message("Considering ",length(genes)," DE genes")
	}else{
		stop("Please provide a list of genes as df rownames or as a vector of gene names")
	}
	if(which.ids %in% c('entrez','symbol')){
		if(which.ids == 'entrez'){
			sig.gene = bitr(genes,
					 fromType="SYMBOL",
					 toType = "ENTREZID",
				 	OrgDb = org.Hs.eg.db)
			sig.gene = sig.gene[!duplicated(sig.gene$SYMBOL),]
			rownames(sig.gene) = sig.gene$SYMBOL
			sig.gene = sig.gene$ENTREZID
			x = enricher(gene=sig.gene, pvalueCutoff = 0.05,TERM2GENE=paths[,c(1,2)])
			x = setReadable(x, 'org.Hs.eg.db', 'ENTREZID')
		}else if(which.ids == 'symbol'){
			x = enricher(gene=genes, pvalueCutoff = 0.05,TERM2GENE=paths[,c(1,3)])
		}
		message(nrow(x)," Reactome terms enriched")
		x = pairwise_termsim(x)
		return(x)
	}else{
		stop("which.ids parameter should be set as either 'entrez' or 'symbol'")
	}
}	


#Alternative to the enricher() function of clusterProfiler
#genes and background should be vectors of gene identifiers
#sigs should be a names list of signatures(each list element contains the gene identifiers of each signature)
my_enricher = function(genes = NULL,sigs = NULL,background=NULL,padj_thresh=0.05){
	hyper = lapply(1:length(sigs),function(x){
	n = intersect(sigs[[x]],background) %>% length()
	common = intersect(sigs[[x]],genes) %>% length()
	geneRatio = paste0(common,"/",length(genes))
	bgRatio = paste0(n,"/",length(background))
	common_genes = paste0(intersect(sigs[[x]],genes),collapse=',')
	pval = phyper(common-1,n,length(background)-n,length(genes),lower.tail=F)
	df=data.frame(Description = names(sigs)[x],geneRatio = geneRatio,bgRatio=bgRatio,overlap=common_genes,Count = common,pval=pval)
	return(df)})
	all = do.call('rbind',hyper)
	all$p.adjust = p.adjust(all$pval,method='BH')
	all = all[all$p.adjust < padj_thresh,]
	return(all)
}

#Dotplot of enriched Reactome/Hallmark/... terms. 
#The input should be a data.frame reporting the pathway enrichment results obtained through the clusterProfiler::enricher() function
enrich_dotplot = function(x=NULL,top_n=25,title=""){
	#Define scale function to log-scale adj pvalues (-log10) 
	neg_log10 = function(x){-log10(x)}
	elevate_10 = function(x){10^(-x)}
	tn = scales::trans_new("neg_log10",
			       function(x) -log10(x),
			       function(y) 10^(-y),
			       domain=c(Inf, 0))
	if(is.data.frame(x) == F){
		x = as.data.frame(x)
	}
	if(nrow(x) > 0){
		x = mutate(x,qscore = -log10(p.adjust),label=round(p.adjust,2)) %>% arrange(qscore)
		x = tail(x,n)
		x$Description = factor(x$Description,levels=unique(x$Description))
		#Define which adj pvalues breaks you want to show on the x axis
		labels.min = min(x$p.adjust)
		labels.max = max(x$p.adjust)
		#log-transform padj and calculate the range btw min and max padj, then split the range into 5 breaks
		labels.range = neg_log10(labels.max) - neg_log10(labels.min)
		labels.unit = labels.range/5
		#Store the list of padj breaks to show in the plot
		labels.final = (neg_log10(labels.min)+seq(labels.unit,labels.range,labels.unit)) %>% elevate_10()
		x$padj_2 = prettyNum(x$p.adjust, scientific = T, digits = 2)
	        g = ggplot(x,aes(x = p.adjust,y=Description))+
		theme_linedraw()+
		geom_point(stat='identity',aes(color=p.adjust,size=Count))+
		#xlab(expression(-log[10](p.adjust)))+
		xlab('p.adjust')+
		scale_colour_gradient(high="dodgerblue3",low="indianred3")+
		theme(axis.text.y = element_text(size=7,color='black',face='bold'))+
		theme(axis.text.x = element_text(size=12,color='black',angle = 45,vjust=0.5))+
		theme(axis.title = element_text(size=14,color='black'))+
		guides(size  = guide_legend(order = 1),
		       color = guide_colorbar(order = 2))+
		scale_x_continuous(trans=tn,
				   breaks=(c((labels.min),labels.final,(labels.max))),
				   labels=prettyNum(c((labels.min),labels.final,(labels.max)), scientific = T, digits = 2))+
		ggtitle(title)+theme(plot.title=element_text(hjust = 0.5))
		g
	}else{
		ggplot()+theme_linedraw()+
		ggtitle(title)+theme(plot.title=element_text(hjust = 0.5))
        }
}
	
#Function to plot multiple dotplot given a list of dataframes/vectors of DE genes
multi_enrich_dotplot = function(gene_lists = NULL,which.ids = c('entrez','symbol'),paths=NULL,top_n=25){
	all_enrich = lapply(gene_lists,function(x){
		if(is.data.frame(x)){
			if(nrow(x) > 0){
				x_enrich = try(wrap_enrich(x,which.ids = which.ids,paths=paths))
			}else{
				return(message("No DE genes found"))
			}
		}else if(is.vector(x)){
			if(length(x) > 0){
				x_enrich = try(wrap_enrich(x,which.ids = which.ids,paths=paths))
			}else{
				return(message("No DE genes found"))
			}
		}
		if(!is(x_enrich,'try-error')){
			x_dot = enrich_dotplot(x_enrich,top_n=top_n)
		}
		if(exists('x_dot') == T){
			x_dot
		}else{
			return(message("No enriched terms found among the chosen pathways"))
		}
	       }
	)
	return(all_enrich)	
}


#' **Volcano plot based on Seurat DE markers results**
#'
#' @param seurat_de data.frame reporting DE analysis results from Seurat FindMarkers/FindAllMarkers functions
#' @param clust numeric,cluster ID to select if FindAllMarkers has been used
#' @param padj_thresh numeric, p.adjust threshold to consider to plot and label significant DE genes
#' @param logfc_thresh numeric, log2FC threshold to consider to plot and label significant DE genes
#' @param label.default boolean, whether or not label significant DE genes with default settings (top 25 UP and top 25 DOWN genes)
#' @param label.genes vector of gene symbols to label in the volcano plot. It will be applied only if label.default is FALSE
#' @param plot.title string, title to add at the top of the plot. If NULL, no default title is added to the plot
#' @return ggplot2 volcano plot
#' @example de_volcano(seurat_de=tum_hm_de,clust=0,padj_thresh=0.01,logfc_thresh=0.5,label.default=T)
de_volcano = function(seurat_de = NULL,
		      clust=0,
		      padj_thresh = 0.01,logfc_thresh = 0.25,
		      label.default=T,label.genes = c(),
		      plot.title = NULL,col.up='goldenrod2',col.down='steelblue3'){
	if(length(intersect(colnames(seurat_de),"cluster")) > 0){
		seurat_de = seurat_de %>% filter(cluster == clust)
	}
	if(length(intersect(colnames(seurat_de),"gene")) == 0){
		seurat_de$gene = rownames(seurat_de)
	}
	seurat_de$FC = 2^seurat_de$avg_log2FC
	seurat_de = seurat_de %>% arrange(desc(FC))
	seurat_de$label = ''
	if(label.default ==T){
		seurat_de[1:25,]$label = seurat_de[1:25,]$gene
		n = nrow(seurat_de)
		seurat_de[(n-24):n,]$label = seurat_de[(n-24):n,]$gene
		seurat_de[abs(seurat_de$avg_log2FC) < logfc_thresh | 
			  seurat_de$p_val_adj > padj_thresh,]$label = ''
	}else if(label.default ==F & length(label.genes) > 0){
		seurat_de[seurat_de$gene %in% label.genes,]$label = seurat_de[seurat_de$gene %in% label.genes,]$gene
		seurat_de[abs(seurat_de$avg_log2FC) < logfc_thresh | 
			  seurat_de$p_val_adj > padj_thresh,]$label = ''
	}
	seurat_de$color = ifelse(
				 (seurat_de$avg_log2FC) > logfc_thresh & seurat_de$p_val_adj < padj_thresh,'UP',
				 ifelse(
					seurat_de$avg_log2FC < -1*(logfc_thresh) & seurat_de$p_val_adj < padj_thresh,'DOWN','N.S.')
				 )
        seurat_de$log10padj = -log10(seurat_de$p_val_adj)
	g = ggplot(seurat_de,aes(x=FC,y=log10padj,fill=color,label=label))+
	theme_minimal()+geom_point(shape=21,alpha=0.75,size=3)+
	scale_fill_manual(values=c('UP' = col.up,
				   'N.S.' = 'gray40',
				   'DOWN' = col.down)
	)+
	geom_hline(yintercept = -log10(padj_thresh),linetype='dashed',color ='black')+
        geom_vline(xintercept = 2^(logfc_thresh),linetype='dashed',color ='black')+
	geom_vline(xintercept = 2^(-1*(logfc_thresh)),linetype='dashed',color ='black')+
	annotate('text',x = 2,y = -log10(padj_thresh/1000),label= paste0('FDR = ',padj_thresh*100,'%',size=4))+
	ylab('-log10padj')+labs(fill = 'Sign')+xlab('FC')+
	ggrepel::geom_label_repel(max.overlaps = Inf,seed=1)
	if(!is.null(plot.title)){
		g = g + ggtitle(plot.title)+
		theme(plot.title = element_text(color='black',hjust=0.5))
	}
      g
}
#Same as de_volcano, but with pct.1 on the y axis
de_volcano2 = function(seurat_de = NULL,
		       clust=0,padj_thresh = 0.01,
		       logfc_thresh = 0.25,label.default=T,n_label=25,
		       label.genes = c(),col.up='goldenrod2',col.down='steelblue3'){
	if(length(intersect(colnames(seurat_de),"cluster")) > 0){
		seurat_de = seurat_de %>% filter(cluster == clust)
	}
	if(length(intersect(colnames(seurat_de),"gene")) == 0){
		seurat_de$gene = rownames(seurat_de)
	}
	seurat_de$FC = 2^seurat_de$avg_log2FC
        seurat_de = seurat_de %>% arrange(desc(FC)) 
	seurat_de$label = ''
	seurat_de= seurat_de %>% filter(!is.na(p_val_adj))
	if(label.default ==T){
		seurat_de[1:n_label,]$label = seurat_de[1:n_label,]$gene
		n = nrow(seurat_de)
		seurat_de[(n-(n_label-1)):n,]$label = seurat_de[(n-(n_label-1)):n,]$gene
		seurat_de[abs(seurat_de$avg_log2FC) < logfc_thresh | 
			  seurat_de$p_val_adj > padj_thresh,]$label = ''
	}else if(label.default ==F & length(label.genes) > 0){
		seurat_de[seurat_de$gene %in% label.genes,]$label = seurat_de[seurat_de$gene %in% label.genes,]$gene
		seurat_de[abs(seurat_de$avg_log2FC) < logfc_thresh | 
			  seurat_de$p_val_adj > padj_thresh,]$label = ''
	}
	seurat_de$color = ifelse(
		(seurat_de$avg_log2FC) > logfc_thresh & seurat_de$p_val_adj < padj_thresh,'UP',
		ifelse(
		       seurat_de$avg_log2FC < -1*(logfc_thresh) & seurat_de$p_val_adj < padj_thresh,'DOWN','N.S.')
		)
	g = ggplot(seurat_de,aes(x=FC,y=pct.1,fill=color,label=label))+
	theme_minimal()+geom_point(shape=21,alpha=0.75,size=3)+
	scale_fill_manual(values=c('UP' = col.up,
				   'N.S.' = 'gray40',
				   'DOWN' = col.down))+
	geom_vline(xintercept = 2^(logfc_thresh),linetype='dashed',color ='black')+
	geom_vline(xintercept = 2^(-1*(logfc_thresh)),linetype='dashed',color ='black')+
	ylab('% Expressing cells')+labs(fill = 'Sign')+xlab('FC')+
	scale_y_continuous(breaks=seq(0,1,0.1),labels=seq(0,1,0.1))+
	ggrepel::geom_label_repel(max.overlaps = Inf,seed=1)
	g
}


#' **Combined UP and Down dotplot for a specified pathway enrichment analysis
#'
#' @param x data.frame reporting DE analysis results from Seurat FindMarkers/FindAllMarkers functions. Default is NULL
#' @param n numeric, maximum number of significant pathways to report in the plot. Default is 50
#' @param title, character string to use as plot title. Default is an empty string
#' @param logfc.pos numeric, log2FC threshold to apply to DE UP genes
#' @param logfc.neg numeric, log2FC threshold to apply to DE DOWN genes
#' @param padj.thresh numeric,p-adjusted threshold to apply to DE genes. Default is 0.01
#' @param enrichment.thresh numeric,p-adjusted threshold to apply to enriched pathways. Default is 0.05
#' @param filter.paths boolean, whether or not discard from the final plot those paths which are enriched for both UP and DOWN genes. Default is FALSE
#' @param background, character vector indicating the list of background genes to use for the enrichment test. Default is NULL
#' @param sigs, list of vectors of gene symbol corresponding to the gene set to test for enrichment. Default is NULL
#' @param select.paths, character vector indicating a selected list of gene set names to be reported in the final plot. Default is NULL
#' @return a two dodged dotplots, each one reporting gene sets enriched for either DOWN- or UP- regulated genes from the tested comparison
#' @example single_dot(x=de_hm,n=50,title='HM',padj.thresh=0.01,enrichment.thresh=0.1,background=rownames(de_hm),sigs=hallmark)

single_dot = function(x = NULL,n=50,title = "", 
		      logfc.pos=log2(1.1),logfc.neg=log2(0.9),
		      padj.thresh=0.01,enrichment.thresh=0.05,
		      filter.paths=F,
		      background = NULL,sigs=NULL,
		      select.paths = NULL){
	#Define scale function to log-scale adj pvalues (-log10) 
	neg_log10 = function(x){-log10(x)}
	elevate_10 = function(x){10^(-x)}
	tn = scales::trans_new("neg_log10",
			       function(x) -log10(x),
			       function(y) 10^(-y),
			       domain=c(Inf, 0))
	x_neg = my_enricher(rownames(x[x$avg_log2FC < logfc.neg & x$p_val_adj < padj.thresh,]),background = background,sigs=sigs,padj_thresh = enrichment.thresh)
	x_pos = my_enricher(rownames(x[x$avg_log2FC > logfc.pos & x$p_val_adj < padj.thresh,]),background = background,sigs=sigs,padj_thresh = enrichment.thresh)
	
	if(nrow(x_pos) > 0){
		x_pos = mutate(x_pos,qscore = -log10(p.adjust),label=round(p.adjust,2)) %>% arrange(qscore)
		x_pos = tail(x_pos,n)
	}
	if(nrow(x_neg) > 0){
		x_neg = mutate(x_neg,qscore = -log10(p.adjust),label=round(p.adjust,2)) %>% arrange(qscore)
	        x_neg = tail(x_neg,n)
	}
	if((nrow(x_pos) != 0) | (nrow(x_neg) != 0)){
		if((nrow(x_pos) != 0)){
			x_pos = x_pos %>% arrange(qscore)
			x_pos$Sign = 'UP'
		}
		if((nrow(x_neg) != 0)){
			x_neg = x_neg %>% arrange(desc(qscore))
			x_neg$Sign = 'DOWN'
		}
		x_tot = rbind(x_pos,x_neg)
		x_tot$Description = factor(x_tot$Description,levels=unique(x_tot$Description))
		paths = rev(as.character(unique(x_tot$Description)))
		if(filter.paths==T){
			if((nrow(x_pos) != 0)){
				up_paths = (x_tot %>% count(Description,Sign) %>% filter(Sign == 'UP'))$Description
		          }
			if((nrow(x_neg) != 0)){
				down_paths = (x_tot %>% count(Description,Sign) %>% filter(Sign == 'DOWN'))$Description
			}
			common_paths = intersect(up_paths,down_paths)
			x_tot = x_tot[!x_tot$Description %in% common_paths,]
		}
		if(!is.null(select.paths)){
			x_tot = x_tot[!x_tot$Description %in% select.paths,]
		}        
		#Define which adj pvalues breaks you want to show on the x axis
		labels.min = min(x_tot$p.adjust)
		labels.max = max(x_tot$p.adjust)
		#log-transform padj and calculate the range btw min and max padj, then split the range into 5 breaks
		labels.range = neg_log10(labels.max) - neg_log10(labels.min)
		labels.unit = labels.range/5
		#Store the list of padj breaks to show in the plot
		labels.final = (neg_log10(labels.min)+seq(labels.unit,labels.range,labels.unit)) %>% elevate_10()
		#x_tot$Sign = factor(x_tot$Sign, levels = c('UP','DOWN'))
		g=ggplot(x_tot,aes(x=p.adjust,y=Description,fill=Sign))+
		theme_linedraw()+
		geom_errorbarh(stat='identity',aes(xmin = 1, xmax = p.adjust),color='black',
			       height=0)+#,position=position_dodge(width=0.7))+
		geom_point(stat='identity',aes(size=Count),shape=21)+#,position=position_dodge(width=0.7))+
		xlab('p.adjust')+
		scale_fill_manual(values = c('DOWN' = 'steelblue3','UP' = 'goldenrod2'))+
		#scale_colour_gradient(high="dodgerblue3",low="indianred3")+
		theme(axis.text.y = element_text(size=7,color='black',face='bold'))+
		theme(axis.text.x = element_text(size=12,color='black',angle = 45,vjust=0.5))+
		theme(axis.title = element_text(size=14,color='black'))+
		guides(size  = guide_legend(order = 1),
		       fill = 'none',
		       color = 'none'
		       )+
		scale_x_continuous(trans=tn,
				   breaks=(c((labels.min),
					     labels.final,
					     (labels.max))),
				   labels=prettyNum(c((labels.min),
						      labels.final,
						      (labels.max)), 
						    scientific = T, 
						    digits = 2
						    )
				   )+
		ggtitle(title)+theme(plot.title=element_text(hjust = 0.5))+
		facet_grid(cols=vars(Sign))
		return(list(paths = paths,plot=g))
	}else{
		g=ggplot()+theme_linedraw()+
		ggtitle(title)+theme(plot.title=element_text(hjust = 0.5))
		return(list(paths = paths,plot=g))
	}
}

#' **Run fgsea pre-ranked analysis on a list of DE genes (Seurat format expected)**
#'
#' @param seurat_de data.frame reporting DE analysis results from Seurat FindMarkers/FindAllMarkers functions. Default is NULL
#' @param clust numeric,cluster ID to select if FindAllMarkers has been used. Default is 0
#' @param sigs list, vectors of gene symbols (i.e. signatures) to test for differential activity. Default is NULL
#' @param logfc.filter boolean, whether or not apply an initial filter on the log2FC of the DE genes list. Default is FALSE
#' @param logfc.thresh numeric, log2FC threshold to apply to DE genes if logfc.filter parameter is TRUE. Default is log2(1.1)
#' @param padj.filter boolean, whether or not apply an initial filter on the adjusted p-values of the DE genes list. Default is FALSE 
#' @param padj.threshold numeric,p-adjusted threshold to apply to DE genes if padj.filter parameter is TRUE. Default is 0.01
#' @param minSize numeric, minimum number of genes of the tested gene set that must be included in the DE genes list. Default is 5
#' @param maxSize numeric, maximum number of genes of the tested gene set that must be included in the DE genes list. Default is 500
#' @param seed numeric, random seed to set before running fgsea function (for reproducibility). Default is 42
#' @param plot boolean, whether or not making a barplot of the top altered pathways based on fgsea NES and padj. Default is FALSE
#' @param plot.signif.thresh numeric, padj threshold for significant pathways from fgsea analysis if plot parameter is TRUE
#' @param plot.title string, title to add at the top of the plot when plot=TRUE. If NULL, no default title is added to the plot
#' @return a list containing an fgsea result data.frame and a ggplot2 object if plot is TRUE. if FALSE, only the data.frame is returned
#' @example fgsea_seurat(seurat_de=tum_hm_de,clust=0,sigs=hallmark,padj.filter=T,minSize=5,plot=T,plot.signif.thresh=0.1)
fgsea_seurat = function(seurat_de=NULL,
			 clust = 0,
			 sigs = NULL,
			 logfc.filter=F,logfc.thresh=log2(1.1),
			 padj.filter=F,padj.thresh=0.01,
			 minSize=5,maxSize=500,
			 seed = 42,
			 plot=F,
			 plot.signif.thresh=0.25,
			 plot.title = NULL){
	if(length(intersect(colnames(seurat_de),"cluster")) > 0){
		seurat_de = seurat_de %>% filter(cluster == clust)
	}
	if(length(intersect(colnames(seurat_de),"gene")) == 0){
		seurat_de$gene = rownames(seurat_de)
	}
	if(logfc.filter ==T){
		seurat_de = seurat_de %>% filter(abs(avg_log2FC) > logfc.thresh)
	}
	if(padj.filter == T){
		seurat_de = seurat_de %>% filter(p_val_adj < padj.thresh)
	}
	seurat_de = seurat_de %>% arrange(desc(avg_log2FC))
	stats = seurat_de$avg_log2FC
	names(stats) = seurat_de$gene
	set.seed(seed)
	res = fgsea::fgsea(pathways = sigs,stats=stats,eps=0.0,minSize=minSize,maxSize=maxSize)
	if(nrow(res[res$padj < plot.signif.thresh,]) > 0){
		    res = res[res$padj < plot.signif.thresh,]
	}
	if(plot == T){
		res = as.data.frame(res)
		res = res[order(res$padj),]
		res = res[!is.na(res$pval),]
		if(nrow(res) > 25){
			res = head(res,50)
		}
		res = res[order(res$NES,decreasing=F),]
		res$pathway = factor(res$pathway,levels=res$pathway)
		res$color = 'N.S'
		if(nrow(res[res$padj < plot.signif.thresh & res$NES > 0,]) > 0){
		res[res$padj < plot.signif.thresh & res$NES > 0,]$color = 'UP'
		}
		if(nrow(res[res$padj < plot.signif.thresh & res$NES < 0,]) > 0){
		res[res$padj < plot.signif.thresh & res$NES < 0,]$color = 'DOWN'
		}
		g=ggplot(res,aes(x=pathway,y=NES,fill=color))+
		theme_classic()+
		geom_bar(stat='identity')+
		scale_fill_manual(values=c('UP' = 'goldenrod2',
					   'N.S' = 'gray40',
					   'DOWN' = 'steelblue3')
		)+
		coord_flip()+
		labs(fill = 'Regulation')+
		xlab('Pathway')
		if(!is.null(plot.title)){
		   g = g + ggtitle(plot.title)+
		   theme(plot.title = element_text(color='black',hjust=0.5))
		}
		return(list(df=res,plot=g))
	}else{
		return(res)
	}
}

