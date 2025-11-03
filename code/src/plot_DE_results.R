library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
utils = snakemake@params[["utils"]]
enrich = snakemake@params[["enrich_utils"]]
source(enrich)
source(utils)

#Define enrichment plot function
single_dot2 = function(x = NULL,n=50,title = "", 
		       logfc.pos=log2(1.1),logfc.neg=log2(0.9),
		       padj.thresh=0.01,enrichment.thresh=0.05,
		       filter.paths=F,
		       background = NULL,sigs=NULL,
		       select.paths = NULL,col.up='goldenrod2',col.down='steelblue3'){
	#Define scale function to log-scale adj pvalues (-log10) 
	neg_log10 = function(x){-log10(x)}
	elevate_10 = function(x){10^(-x)}
	tn = scales::trans_new("neg_log10",
			       function(x) -log10(x),
			       function(y) 10^(-y),
			       domain=c(Inf, 0))
	x_neg = my_enricher(rownames(x[x$avg_log2FC < logfc.neg & x$p_val_adj < padj.thresh,]),background = background,sigs=sigs,padj_thresh = enrichment.thresh)
	x_neg = x_neg %>% mutate(weight = round(Count / as.numeric(sapply(bgRatio,function(x){strsplit(x,'/')[[1]][[1]]})),1))
	x_pos = my_enricher(rownames(x[x$avg_log2FC > logfc.pos & x$p_val_adj < padj.thresh,]),background = background,sigs=sigs,padj_thresh = enrichment.thresh)
        x_pos = x_pos %>% mutate(weight = round(Count / as.numeric(sapply(bgRatio,function(x){strsplit(x,'/')[[1]][[1]]})),1))
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
		geom_errorbarh(stat='identity',aes(xmin = enrichment.thresh, xmax = p.adjust),color='black',
			       height=0)+#,position=position_dodge(width=0.7))+
		geom_point(stat='identity',aes(size=weight),shape=21)+#,position=position_dodge(width=0.7))+
		scale_size(range=c(1,3),breaks=c(0.4,0.6,0.8))+
		xlab('p.adjust')+
		scale_fill_manual(values = c('DOWN' = col.down,'UP' = col.up))+
		theme(axis.text.y = element_text(size=4.5,color='black',face='bold'))+
		theme(axis.text.x = element_text(size=7,color='black',angle = 45,vjust=0.5))+
		theme(axis.title = element_text(size=9,color='black'))+
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

#Load Reactome pathways (levels 3-5)
react_df = read.delim(snakemake@params[["react"]])
react_df_list = split(react_df$gene_symbol,react_df$pathwayName)
react_df_list = lapply(react_df_list,hallmark2cosmx,collapse=T)
#Keep genes associated with antigen presentation/processing (AP) or Interferon (IFN) signaling (use keywords)
mhc_genes = (react_df %>% filter(str_detect(pathwayName,'(?i)antigen presentation|antigen processing|Interferon')))$gene_symbol %>% unique()
mhc_genes_cosmx = hallmark2cosmx(mhc_genes,collapse=T)

#Load Epi-T/Epi-Only DE results
all_de_epi = readRDS('intermediate_outputs/DE_res/epiT_vs_epiOnly/InSituType/Cosmx_ALL_InSitu_epiT_epiOnly_4Neigh_DE_res_TNKsubset.rds')
up_epi = rownames(all_de_epi %>% filter(p_val_adj < 0.01) %>% filter(avg_log2FC > log2(1.1)))
mhc_epi = intersect(mhc_genes_cosmx,up_epi)
#Load TAM-T/TAM-Other DE results
all_de_macro = readRDS('intermediate_outputs/DE_res/myelT_vs_myelNoT/Cosmx_ALL_myelT_myelNoT_DE_res.rds')
up_macro = rownames(all_de_macro %>% filter(p_val_adj < 0.01) %>% filter(avg_log2FC > log2(1.1)))
mhc_macro = intersect(mhc_genes_cosmx,up_macro)
#Load IE-T/Ts-T DE results
all_de_t = readRDS('intermediate_outputs/DE_res/Ttouch_vs_Tnotouch/InSituType/Cosmx_ALL_Ttouch_Tnotouch_DE_res_TNKsubset.rds')
up_t = rownames(all_de_t %>% filter(p_val_adj < 0.01) %>% filter(avg_log2FC > log2(1.1)))
#Load HCT116 coculture DE results
hct116_de = readRDS('intermediate_outputs/HCT116/res/HCT116_EpiT_EpiOnly_DE_res_NCG7.rds')
up_hct = rownames(hct116_de %>% filter(log2FoldChange > 1,padj < 0.05))
up_hct = data.frame(gene_name = up_hct,cosmx_name = hallmark2cosmx(up_hct,collapse=F))
mhc_hct = intersect(mhc_genes,up_hct$gene_name)

#Find common AP/IFN genes upregulated across the three comparisons
mhc_common = intersect(mhc_epi,mhc_macro)
mhc_common = intersect(mhc_common,up_hct$cosmx_name)

#Plot TAM-T/TAM-Other DE volcano plot
pdf(snakemake@output[["volcano_tam"]],width=4,height=1.8)
de_volcano(all_de_macro,label.default = F, logfc_thresh_up = log2(1.1),min.pct = 0,logfc_thresh_down = log2(0.9),label.genes=c('CXCL10','CXCL9',mhc_common),pt.size=0.5,label.size=1.3,col_up = '#83C65D',col_down='#3B4E36',outType='seurat')+
guides(color='none')+
ggtitle('TAMs-T (13927) vs TAMs-Other (59404)')+
annotate('text',size=2,x=0.68,y=193,label = paste0('AP-IFN genes/UP genes = ',length(mhc_macro),'/',length(up_macro)))
dev.off()
#Plot Epi-T/Epi-Only DE volcano plot
pdf(snakemake@output[["volcano_epi"]],width=4,height=1.8)
de_volcano(all_de_epi,label.default = F, logfc_thresh_up = log2(1.1),min.pct=0,logfc_thresh_down = log2(0.9),label.genes=c(mhc_common),pt.size=0.5,label.size=1.3,col_up = 'cyan',col_down='blue',outType='seurat')+
guides(color='none')+
ggtitle('Epi-T (37899) vs Epi-Only (186371)')+
annotate('text',size=2,x=0.4,y=250,label = paste0('AP-IFN genes/UP genes = ',length(mhc_epi),'/',length(up_epi)))
dev.off()
#Plot HCT116+TCD8/HCT116-Only DE volcano plot
pdf(snakemake@output[["volcano_hct"]],width=4,height=3)
de_volcano(hct116_de,logfc_thresh_up = 1,logfc_thresh_down = -1,min.pct=0,outType = 'deseq',label.default = F,label.genes = c(up_hct[up_hct$cosmx_name %in% mhc_common,]$gene_name),pt.size=0.5,label.size=1.5,col_up = 'cyan',col_down = 'blue')+
guides(color = 'none')+
ggtitle('HCT116+TCD8 vs HCT116-Only')+theme(axis.text = element_text(size=5))+
annotate('text',size=2,x=75,y=125,label = paste0('AP-IFN genes/UP genes = ',length(mhc_hct),'/',length(up_hct$gene_name)))
dev.off()

#List of Cytotoxic-23 genes + ITGAE/ENTPD1
cd8_genes = c('RORA','CTSW','RUNX3','DUSP2','KLRB1','KLRF1','SIGIRR','CD8B','CD8A','GZMB','GZMA','GNLY','GZMH','PDCD1','ITGAE','ENTPD1','PRF1','KLRK1','NKG7','LAG3','HAVCR2','CCL3/L1/L3')
#Plot IE-T/Ts-T DE volcano plot
pdf(snakemake@output[["volcano_t"]],width=4,height=1.8)
#Remove CXCL14 from the plot just to make it less 'squeezed'
de_volcano(all_de_t[rownames(all_de_t) != 'CXCL14',],label.default = F,min.pct=0,label.genes=intersect(cd8_genes,up_t),label.size = 1.3,pt.size=0.5,outType='seurat',col_up = '#FFA8FF',col_down='#FF1AFF')+
guides(color='none')+
ggtitle('IE-T (11099) vs Ts-T (13468)')
dev.off()



###TAM pathway enrichment boxplot
pdf(snakemake@output[["path_dot_tam"]],width=4,height=2)
single_dot2(x=all_de_macro,background=rownames(all_de_macro),sigs=react_df_list,enrichment.thresh=0.01,title = '',col.up='#83C65D',col.down='#3B4E36')$plot
dev.off()
