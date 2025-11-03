library(dplyr)
library(stringr)
library(gpglot2)

myeloid_path = snakemake@params[["myeloid_path"]]
tnk_path = snakemake@params[["tnk_path"]]

#Load Myeloid subcluster annotations and plot subtype proportions for each sample
cr21_sup = readRDS(list.files(myeloid_path,pattern='CR21',full.names=T))$clust
cr21_sup = data.frame(cell_ID = names(cr21_sup),cellType=cr21_sup,sample='CR21')
cr48_sup = readRDS(list.files(myeloid_path,pattern='CR48',full.names=T))$clust
cr48_sup = data.frame(cell_ID = names(cr48_sup),cellType=cr48_sup,sample='CR48')
cr36_sup = readRDS(list.files(myeloid_path,pattern='CR36',full.names=T))$clust
cr36_sup = data.frame(cell_ID = names(cr36_sup),cellType=cr36_sup,sample='CR36')
uh20_sup = readRDS(list.files(myeloid_path,pattern='UH20',full.names=T))$clust
uh20_sup = data.frame(cell_ID = names(uh20_sup),cellType=uh20_sup,sample='UH20')

pdf(snakemake@output[["myeloid_bar"]],width=5,height=4)
all_sup = purrr::reduce(list(cr21_sup,cr48_sup,cr36_sup,uh20_sup),rbind)
all_sup$cellType = factor(all_sup$cellType,levels=c('Mast','Granulo','DC','Mono','Macro'))
all_sup$sample = factor(all_sup$sample,levels=c('CR21','CR48','CR36','UH20'))
all_sup %>% 
dplyr::count(sample,cellType) %>% 
group_by(sample) %>% 
mutate(prop=n/sum(n)) %>% 
ungroup() %>% 
ggplot(.,aes(x=sample,y=100*prop,fill=cellType))+
theme_classic()+geom_bar(stat='identity')+
scale_fill_manual(values=c12)+scale_y_continuous(breaks=c(0,50,100),labels=c('0','','100'))+
ylab('Myeloid cells (%)')+xlab('')+labs(fill = 'Cell Type')
dev.off()

#Load T/NK subcluster annotations and plot subtype proportions for each sample
cr21_sup = readRDS(list.files(tnk_path,pattern='CR21',full.names=T))$clust
cr21_sup = data.frame(cell_ID = names(cr21_sup),cellType=cr21_sup,sample='CR21')
cr48_sup = readRDS(list.files(tnk_path,pattern='CR48',full.names=T))$clust
cr48_sup = data.frame(cell_ID = names(cr48_sup),cellType=cr48_sup,sample='CR48')
cr36_sup = readRDS(list.files(tnk_path,pattern='CR36',full.names=T))$clust
cr36_sup = data.frame(cell_ID = names(cr36_sup),cellType=cr36_sup,sample='CR36')
uh20_sup = readRDS(list.files(tnk_path,pattern='UH20',full.names=T))$clust
uh20_sup = data.frame(cell_ID = names(uh20_sup),cellType=uh20_sup,sample='UH20')

pdf(snakemake@output[["tnk_bar"]],width=5,height=4)
all_sup = purrr::reduce(list(cr21_sup,cr48_sup,cr36_sup,uh20_sup),rbind)
all_sup$cellType = factor(all_sup$cellType,levels=c('Tgd','NK','TCD4','TCD8'))
all_sup$sample = factor(all_sup$sample,levels=c('CR21','CR48','CR36','UH20'))
all_sup %>% 
dplyr::count(sample,cellType) %>% 
dplyr::group_by(sample) %>% 
dplyr::mutate(prop=n/sum(n)) %>% 
dplyr::ungroup() %>% 
ggplot(.,aes(x=sample,y=100*prop,fill=cellType))+
theme_classic()+geom_bar(stat='identity')+
scale_fill_manual(values=c12)+scale_y_continuous(breaks=c(0,50,100),labels=c('0','','100'))+
ylab('T/NK cells (%)')+xlab('')+labs(fill = 'Cell Type')
dev.off()


