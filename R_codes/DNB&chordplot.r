library(Seurat)
library(DNBr)
library(dplyr)

stages_BP_all_rna_filtered_harmony <- readRDS( "./all_3_stage_delbf_all_filtered_blood(2).rds")

data <- stages_BP_all_rna_filtered_harmony@assays$RNA@scale.data
meta <- stages_BP_all_rna_filtered_harmony@meta.data$group
names(meta) <- rownames(stages_BP_all_rna_filtered_harmony@meta.data)

a <- DNBcompute(
    data = data, 
    meta = meta, 
    # diffgenes = diffgenes.example
)

b <- DNBfilter(
  DNB_output = a, 
  ntop = 5
)

df_score

df_allscore <- resultAllExtract(
    object = b,
    group = "ef",
    slot = "result"
)
df_allscore

DNBplot(b, ranking = 1, group = "mf", show = TRUE, save_pdf = FALSE, meta_levels = c('ef','mf','lf'))

df_allscore[order(df_allscore$rank_all,decreasing = F),]

b <- readRDS('3_stage_DNB_b_8_24.rds')

df_allscore <- resultAllExtract(
    object = b,
    group = "mf",
    slot = "result"
)
#df_allscore

new_ann  <- read.csv('./annotation new_change.csv', header =  F)
ann <- data.frame(substring(new_ann$V1, 8, 17), new_ann$V6)
colnames(ann) <- c('gene','genename')
ann <- ann[!(ann$genename == ''),]
#ann$genename <- unlist(lapply(strsplit(ann$genename,'_'),function(x){return(x[[1]])}))
ann <- ann[!duplicated(ann$gene),]
rownames(ann) <- ann$gene
head(ann)

func_4 <- function(x){
    
    vec_x <- strsplit(x,split = ',')[[1]]
    gene <- ann[vec_x,'gene'] %>% na.omit()
    id <- ann[vec_x,'genename'] %>% na.omit()
    paste(gene,id,sep = '|') %>% paste(collapse = ' , ')

}

df_allscore$genes <- df_allscore %>%
                    select('genes') %>%
                    apply(MARGIN = 1,FUN = func_4)
df_allscore

df_allscore[order(df_allscore$rank_all,decreasing = F),]

write.table(df_allscore[order(df_allscore$rank_all,decreasing = F),],'3_stage_DNB_moudles_top5_ann.csv',sep = ',',quote = T)

b <- readRDS('3_stage_DNB_b_8_24.rds')

DNBplot(b, resource = 'mf_55', show = TRUE, save_pdf = FALSE, meta_levels = c('ef','mf','lf'))

df_pt <- ScoreExtract(b, resource = 'mf_55')
df_pt

 pp <- function(df, y, meta_levels = NULL,red=2) {
        Names <- "Names"
        df_score <- df[, c(Names, y)]
        rownames(df_score) <- df_score$Names
        df_score <- df_score[meta_levels,]
        colnames(df_score)[2] <- 'y'
        df_score$Names <- factor(df_score$Names, levels = meta_levels)
     
        ggplot(df_score, aes(x = Names, y = y, group = 1)) +
        geom_point() +
        geom_point(aes(x=Names[red],y=y[red], size = 1),colour="red")+
        geom_line() +
        theme_classic(base_size = 15) +
        ggtitle(y) +
        theme(axis.text.x  = element_text(face ="bold", size = 12, color = "black"),
              axis.text.y  = element_text(face ="bold", size = 12, color = "black"),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.position = "none")
    }

#pp(df_pt,'SCORE',meta_levels = c('ef','mf','lf'))

p1 <- pp(df_pt, 'SCORE',meta_levels = c('ef','mf','lf'))
p2 <- pp(df_pt, 'PCC_IN',meta_levels = c('ef','mf','lf'))
p3 <- pp(df_pt, 'PCC_OUT',meta_levels = c('ef','mf','lf'))
p4 <- pp(df_pt, 'SD',meta_levels = c('ef','mf','lf'))
p <- cowplot::plot_grid(p1, p2, p3, p4)

p

ggsave(p,filename = '3_stage_DNB_mf55.pdf',width = 10,height = 10)

fliter <- readRDS('all_3_stage_delbf_all_filtered_blood(2).rds')

rna_df <- as.matrix(fliter@assays$RNA@data)

tmp <- pair_orth[rownames(rna_df),'SYMBOL']

rna_df <- rna_df[!is.na(tmp),]

rownames(rna_df) <- tmp[!is.na(tmp)]

rna_df <- rna_df[!duplicated(rownames(rna_df)),]


fliter@meta.data$cells <- rownames(fliter@meta.data)

meta_data <- fliter@meta.data %>% dplyr::select(c('cells','seurat_clusters')) %>% as.matrix

write.table(rna_df, '/all/cellphonedb_count.txt', sep='\t', quote=F)
write.table(meta_data, '/all/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)



library(Seurat)
library(iTALK)
library(Matrix)
library(dplyr)
library(igraph)
library(ggplot2)
library(circlize)

fliter <- readRDS('all_3_stage_delbf_all_filtered_blood(2).rds')

rna_df <- as.matrix(fliter@assays$RNA@data)

tmp <- pair_orth[rownames(rna_df),'SYMBOL']

rna_df <- rna_df[!is.na(tmp),]

rownames(rna_df) <- tmp[!is.na(tmp)]

rna_df <- rna_df[!duplicated(rownames(rna_df)),]


fliter@meta.data$cells <- rownames(fliter@meta.data)

meta_data <- fliter@meta.data %>% dplyr::select(c('cells','cell_type')) %>% as.matrix

#write.table(rna_df, './cpdb_9_4/all/cellphonedb_count.txt', sep='\t', quote=F)
write.table(meta_data, './cpdb_9_4/all/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

meta_data <- fliter[,fliter@meta.data$group == 'mf']@meta.data %>% 
dplyr::select(c('cells','cell_type')) %>% 
as.matrix
write.table(meta_data, './cpdb_9_4/mf/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

meta_data <- fliter[,fliter@meta.data$group == 'ef']@meta.data %>% 
dplyr::select(c('cells','cell_type')) %>% 
as.matrix
write.table(meta_data, './cpdb_9_4/ef/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

meta_data <- fliter[,fliter@meta.data$group == 'lf']@meta.data %>% 
dplyr::select(c('cells','cell_type')) %>% 
as.matrix
write.table(meta_data, './cpdb_9_4/lf/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

deconvoluted <- read.table('./cpdb_9_4/all/out/deconvoluted.txt',sep='\t',header=T,check.names = F)
count_network <- read.table('./cpdb_9_4/all/out/count_network.txt',sep='\t',header=T,check.names = F)
interaction_count <- read.table('./cpdb_9_4/all/out/interaction_count.txt',sep='\t',header=T,check.names = T)
means <- read.table('./cpdb_9_4/all/out/means.txt',sep='\t',header=T,check.names = F)
pvalues <- read.table('./cpdb_9_4/all/out/pvalues.txt',sep='\t',header=T,check.names = F)
significant_means <- read.table('./cpdb_9_4/all/out/significant_means.txt',sep='\t',header=T,check.names = F)

deconvoluted[grep('complex',deconvoluted$complex_name),]

interaction_count$X <- factor(interaction_count$X,levels = sort(interaction_count$X))

count_network_mtx <- spread(data = count_network,key = 'TARGET',value = 'count')
rownames(count_network_mtx) <- count_network_mtx$SOURCE
count_network_mtx <- count_network_mtx[,2:ncol(count_network_mtx)] %>% as.matrix()

max(count_network_mtx)

library(ComplexHeatmap)

pdf('3_stage_all_cellphonedb_count_heatmap.pdf',width = 11,height = 10)

col_fun = circlize::colorRamp2(c(0, 30, 60), c("#6B8968", "#FBD786", "#91374B"))
ComplexHeatmap::Heatmap(count_network_mtx,
                        col = col_fun,
                        row_names_gp = gpar(fontsize = 12),
                        column_names_gp = gpar(fontsize = 12),
                        rect_gp = gpar(col = "white", lwd = 2),
                       column_names_rot = 45)

dev.off()

means[,5:12]

mean_mtx <- means[,c(2,8,grep('krt8',colnames(means)))]

pvalue_mtx <- pvalues[,c(2,8,grep('krt8',colnames(means)))]

write.table(mean_mtx,'3stage_all_mean_mtx.csv',sep = ',',quote = T,row.names = F)

mean_long  <- mean_mtx %>% tidyr::gather(cells,value,'Epi(gata1a+)|Epi(krt8+)':'smooth muscle cells|Epi(krt8+)')

plot_df <- mean_long %>%  select(-'receptor_a') %>% cbind(p_long$value)
colnames(plot_df) <- c('interacting_pair','x','value','pvalue')
plot_df$logpvalue <- -log10(plot_df$pvalue)
plot_df$logvalue <- log2(plot_df$value)
plot_df[plot_df$pvalue == 0,]$logpvalue <- 3 

plot_df <- plot_df %>% filter(logpvalue > 2)

plot_df$interacting_pair <- stringr::str_replace(string = plot_df$interacting_pair,pattern = '_',replacement = ' | ') 
plot_df$x <- stringr::str_replace(string = plot_df$x,pattern = '\\|',replacement = ' | ') 

cells <- plot_df$x  %>% unique()
order_cells <- c(grep('^Epi\\(krt',cells))
order_cells <-c(cells[order_cells],cells[-order_cells])
plot_df$x <- factor(plot_df$x,levels = order_cells)

change_to_italk_cell <- function(x){
    
    ls <- list()
    
    for (i in unique(x)){
        #message(i)
        cells <- stringr::str_split_fixed(i,pattern = '\\|',n = 2)
        cell_from <- cells[1]
        cell_to <- cells[2]

        ls[[i]] <- c(cell_from,cell_to)
        #message(ls[[i]])

    }
    return(ls)
}

chang_list_cell <- mean_long$cells %>% change_to_italk_cell() %>% as.data.frame(check.names = F) %>% t() %>% as.data.frame()
colnames(chang_list_cell) <- c('cell_from','cell_to')
chang_list_cell <- chang_list_cell %>% mutate(cells=rownames(chang_list_cell))

chang_list_cell

#修正了受配体标记问题
change_to_italk_pair <- function(x){
    
    unqi_x <- x[!duplicated(x$interacting_pair),]
    
    ls <- list()
    
    for (i in 1:nrow(unqi_x)){
        #message(i)
        str_n <- unqi_x$interacting_pair[i]
        str_n <- stringr::str_replace(str_n ,pattern = '_complex',replace = ' complex')
        cells <- stringr::str_split(str_n,pattern = '_',n = 2)[[1]]
        if (length(cells) >2 ){
            message(length(cells))
        }
        if (unqi_x$receptor_a[i] == 'True'){
            cell_from <- cells[2]
            cell_to <- cells[1]
        }else{
            cell_from <- cells[1]
            cell_to <- cells[2]
        }
        ls[[unqi_x$interacting_pair[i]]] <- c(cell_from,cell_to)
        #message(ls[[i]])

    }
    return(ls)
}

chang_list_pair <- mean_long %>% change_to_italk_pair() %>% as.data.frame(check.names = F) %>% t() %>% as.data.frame()
colnames(chang_list_pair) <- c('ligand','receptor')
chang_list_pair <- chang_list_pair %>% mutate(interacting_pair=rownames(chang_list_pair))

chang_list_pair

mean_long <- mean_long %>% 
    merge(chang_list_cell,by='cells') %>% 
    merge(chang_list_pair,by='interacting_pair')

mean_long <- mean_long %>%
    mutate(cell_from_mean_exprs = mean_long$value) %>%
    mutate(cell_to_mean_exprs = mean_long$value)

head(mean_long)

cp_05 <- mean_long %>% filter(value > 0.5)
cp_05

library(igraph)

NetView2<-function(data,color.use = NULL,title.name = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, remove.isolate = FALSE, top = 1,
                            weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex=1,vertex.label.color= "black",
                            edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                            edge.curved=0.2,shape='circle',layout=in_circle(), margin=0.2, vertex.size = NULL,
                            arrow.width=1,arrow.size = 0.2){
  net<-data %>% group_by(cell_from,cell_to) %>% dplyr::summarize(n=n())
  net<-as.data.frame(net,stringsAsFactors=FALSE)
  g<-graph.data.frame(net,directed=TRUE,)
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5

  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  igraph::E(g)$weight <- E(g)$n
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }

  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  igraph::E(g)$loop.angle <- rep(0, length(igraph::E(g)))

  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
  if (!is.null(title.name)) {
    text(0,1.5,title.name, cex = 1.1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  #gg <- recordPlot()
  #return(gg)
}

my10colors <- my36colors <-ggsci::pal_igv(alpha = 1)(17)
net_col= my10colors

scales::show_col(p_cols)

NetView2(cp_05,color.use=p_cols,vertex.label.cex=1,vertex.size.max=6,weight.scale=T)

pdf('3_stage_celltype_cellphonedb_netview_circle_all.pdf',width = 10,height = 10)
NetView2(cp_05,color.use=p_cols,vertex.label.cex=1,vertex.size.max=6,weight.scale=T)
dev.off()

library(circlize)
LRPlot<-function(data,datatype,gene_col=NULL,transparency=0.5,link.arr.lwd=1,link.arr.lty=NULL,link.arr.col=NULL,link.arr.width=NULL,
                 link.arr.type=NULL,facing='clockwise',cell_col=NULL,print.cell=TRUE,track.height_1=uh(2,'mm'),track.height_2=uh(12,'mm'),
                 annotation.height_1=0.01,annotation.height_2=0.01,text.vjust = '0.4cm',...){
  cell_group<-unique(c(data$cell_from,data$cell_to))
  genes<-c(structure(data$ligand,names=data$cell_from),structure(data$receptor,names=data$cell_to))
  genes<-genes[!duplicated(paste(names(genes),genes))]
  genes<-genes[order(names(genes))]
  if(is.null(link.arr.lty)){
    if(datatype=='mean count'){
      link.arr.lty='solid'
    }else if(datatype=='DEG'){
      link.arr.lty=structure(ifelse(data$cell_from_logFC==0.0001,'dashed','solid'),names=paste(data$cell_from,data$receptor))
    }else{
      print('invalid datatype')
    }
  }
  if(is.null(link.arr.col)){
    if(datatype=='mean count'){
      data<-data %>% mutate(link_col='black')
    }else if(datatype=='DEG'){
      data<-data %>% mutate(link_col=ifelse(cell_from_logFC==0.0001,ifelse(cell_to_logFC>0,'#d73027','#00ccff'),
                                            ifelse(cell_to_logFC==0.0001,ifelse(cell_from_logFC>0,'#d73027','#00ccff'),
                                                   ifelse(cell_from_logFC>0,ifelse(cell_to_logFC>0,'#d73027','#dfc27d'),
                                                          ifelse(cell_to_logFC>0,'#9933ff','#00ccff')))))
    }else{
      print('invalid datatype')
    }
  }else{
    data$link_col=link.arr.col
  }
  if(is.null(link.arr.type)){
    if(datatype=='mean count'){
      link.arr.type='triangle'
    }else if(datatype=='DEG'){
      link.arr.type=structure(ifelse(data$cell_to_logFC==0.0001,'ellipse','triangle'),names=paste(data$cell_from,data$receptor))
    }else{
      print('invalid datatype')
    }
  }
  if(is.null(gene_col)){
    #改了配色
    comm_col<-structure(c('#A15600','#99ccff','#ff9999','#ffcc99'),names=c('other','cytokine','checkpoint','growth factor'))
    gene_col<-structure(c(comm_col[data$comm_type],rep('#073c53',length(data$receptor))),names=c(data$ligand,data$receptor))
  }
  if(is.null(cell_col)){
    cell_col<-structure(randomColor(count=length(unique(names(genes))),luminosity='dark'),names=unique(names(genes)))
  }
  if(is.null(link.arr.lwd)){
    data<-data %>% mutate(arr_width=1)
  }else if(max(abs(link.arr.lwd))-min(abs(link.arr.lwd))==0 && all(link.arr.lwd!=0.0001)){
    data<-data %>% mutate(arr_width=ifelse(abs(link.arr.lwd<5),abs(link.arr.lwd),5))
  }else{
    #加了宽度
    data<-data %>% mutate(arr_width=ifelse(link.arr.lwd==0.0001,2,(1+5/(max(abs(link.arr.lwd))-min(abs(link.arr.lwd)))*(abs(link.arr.lwd)-min(abs(link.arr.lwd))))*2))
  }
  if(length(cell_group)!=1){
    gap.degree <- do.call("c", lapply(table(names(genes)), function(i) c(rep(1, i-1), 8)))
  }else{
    gap.degree <- do.call("c", lapply(table(names(genes)), function(i) c(rep(1, i))))
  }
  circos.par(gap.degree = gap.degree)
  if(length(gene_col)==1){
    grid.col=gene_col
  }else{
    grid.col=gene_col[genes]
    names(grid.col)<-paste(names(genes),genes)
  }
  if(is.null(link.arr.width)){
    data<-data %>% mutate(link.arr.width=data$arr_width/10)
  }else if(max(abs(link.arr.width))-min(abs(link.arr.width))==0 && all(link.arr.width!=0.0001)){
    data<-data %>% mutate(link.arr.width=ifelse(abs(link.arr.width)<0.5,abs(link.arr.width),0.5))
  }else{
    data<-data %>% mutate(link.arr.width=ifelse(link.arr.width==0.0001,0.2,(1+5/(max(abs(link.arr.width))-min(abs(link.arr.width)))*(abs(link.arr.width)-min(abs(link.arr.width))))/10))
  }
  chordDiagram(as.data.frame(cbind(paste(data$cell_from,data$ligand),paste(data$cell_to,data$receptor))), order=paste(names(genes),genes),
               grid.col=grid.col,transparency=transparency,directional=1,direction.type='arrows',link.arr.lwd=data$arr_width,link.arr.lty=link.arr.lty,
               link.arr.type=link.arr.type,link.arr.width=data$link.arr.width,link.arr.col=data$link_col,col='#00000000',annotationTrack=c('grid'),preAllocateTracks = list(
                 list(track.height = track.height_1),list(track.height = track.height_2)),annotationTrackHeight = c(annotation.height_1,annotation.height_2),...)

  circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = genes[get.cell.meta.data("sector.numeric.index")]
    circos.text(mean(xlim),mean(ylim),sector.index, col = "black", cex = 0.7, facing = facing, niceFacing = TRUE)
  }, bg.border = 0)

  if(print.cell){
    for(c in unique(names(genes))) {
      gene = as.character(genes[names(genes) == c])
      highlight.sector(sector.index = paste(c,gene), track.index = 1, col = ifelse(length(cell_col)==1,cell_col,cell_col[c]), text = c, text.vjust = text.vjust, niceFacing = TRUE,lwd=1)
    }
  }
  circos.clear()
}

library(networkD3)
library(dplyr)

library(circlize)
LRPlot_t<-function(data,datatype,gene_col=NULL,transparency=0.5,link.arr.lwd=1,link.arr.lty=NULL,link.arr.col=NULL,link.arr.width=NULL,
                 link.arr.type=NULL,facing='clockwise',cell_col=NULL,print.cell=TRUE,track.height_1=uh(2,'mm'),track.height_2=uh(12,'mm'),
                 annotation.height_1=0.01,annotation.height_2=0.01,text.vjust = '0.4cm',...){
  cell_group<-unique(c(data$cell_from,data$cell_to))
  genes<-c(structure(data$ligand,names=data$cell_from),structure(data$receptor,names=data$cell_to))
  genes<-genes[!duplicated(paste(names(genes),genes))]
  genes<-genes[order(names(genes))]
  if(is.null(link.arr.lty)){
    if(datatype=='mean count'){
      link.arr.lty='solid'
    }else if(datatype=='DEG'){
      link.arr.lty=structure(ifelse(data$cell_from_logFC==0.0001,'dashed','solid'),names=paste(data$cell_from,data$receptor))
    }else{
      print('invalid datatype')
    }
  }
  if(is.null(link.arr.col)){
    if(datatype=='mean count'){
      data<-data %>% mutate(link_col='black')
    }else if(datatype=='DEG'){
      data<-data %>% mutate(link_col=ifelse(cell_from_logFC==0.0001,ifelse(cell_to_logFC>0,'#d73027','#00ccff'),
                                            ifelse(cell_to_logFC==0.0001,ifelse(cell_from_logFC>0,'#d73027','#00ccff'),
                                                   ifelse(cell_from_logFC>0,ifelse(cell_to_logFC>0,'#d73027','#dfc27d'),
                                                          ifelse(cell_to_logFC>0,'#9933ff','#00ccff')))))
    }else{
      print('invalid datatype')
    }
  }else{
    data$link_col=link.arr.col
  }
  if(is.null(link.arr.type)){
    if(datatype=='mean count'){
      link.arr.type='triangle'
    }else if(datatype=='DEG'){
      link.arr.type=structure(ifelse(data$cell_to_logFC==0.0001,'ellipse','triangle'),names=paste(data$cell_from,data$receptor))
    }else{
      print('invalid datatype')
    }
  }
  if(is.null(gene_col)){
    #改了配色
    comm_col<-structure(c('#A15600','#99ccff','#ff9999','#ffcc99'),names=c('other','cytokine','checkpoint','growth factor'))
    gene_col<-structure(c(comm_col[data$comm_type],rep('#073c53',length(data$receptor))),names=c(data$ligand,data$receptor))
  }
  if(is.null(cell_col)){
    cell_col<-structure(randomColor(count=length(unique(names(genes))),luminosity='dark'),names=unique(names(genes)))
  }
  if(is.null(link.arr.lwd)){
    data<-data %>% mutate(arr_width=1)
  }else if(max(abs(link.arr.lwd))-min(abs(link.arr.lwd))==0 && all(link.arr.lwd!=0.0001)){
    data<-data %>% mutate(arr_width=ifelse(abs(link.arr.lwd<5),abs(link.arr.lwd),5))
  }else{
    #加了宽度
    data<-data %>% mutate(arr_width=ifelse(link.arr.lwd==0.0001,2,(1+5/(max(abs(link.arr.lwd))-min(abs(link.arr.lwd)))*(abs(link.arr.lwd)-min(abs(link.arr.lwd))))*2))
  }
  if(length(cell_group)!=1){
    gap.degree <- do.call("c", lapply(table(names(genes)), function(i) c(rep(1, i-1), 8)))
  }else{
    gap.degree <- do.call("c", lapply(table(names(genes)), function(i) c(rep(1, i))))
  }
  circos.par(gap.degree = gap.degree)
  if(length(gene_col)==1){
    grid.col=gene_col
  }else{
    grid.col=gene_col[genes]
    names(grid.col)<-paste(names(genes),genes)
  }
  if(is.null(link.arr.width)){
    data<-data %>% mutate(link.arr.width=data$arr_width/10)
  }else if(max(abs(link.arr.width))-min(abs(link.arr.width))==0 && all(link.arr.width!=0.0001)){
    data<-data %>% mutate(link.arr.width=ifelse(abs(link.arr.width)<0.5,abs(link.arr.width),0.5))
  }else{
    data<-data %>% mutate(link.arr.width=ifelse(link.arr.width==0.0001,0.2,(1+5/(max(abs(link.arr.width))-min(abs(link.arr.width)))*(abs(link.arr.width)-min(abs(link.arr.width))))/10))
  }
    return(as.data.frame(cbind(paste(data$cell_from,data$ligand),paste(data$cell_to,data$receptor))))
  chordDiagram(as.data.frame(cbind(paste(data$cell_from,data$ligand),paste(data$cell_to,data$receptor))), order=paste(names(genes),genes),
               grid.col=grid.col,transparency=transparency,directional=1,direction.type='arrows',link.arr.lwd=data$arr_width,link.arr.lty=link.arr.lty,
               link.arr.type=link.arr.type,link.arr.width=data$link.arr.width,link.arr.col=data$link_col,col='#00000000',annotationTrack=c('grid'),preAllocateTracks = list(
                 list(track.height = track.height_1),list(track.height = track.height_2)),annotationTrackHeight = c(annotation.height_1,annotation.height_2),...)

  circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = genes[get.cell.meta.data("sector.numeric.index")]
    circos.text(mean(xlim),mean(ylim),sector.index, col = "black", cex = 0.7, facing = facing, niceFacing = TRUE)
  }, bg.border = 0)

  if(print.cell){
    for(c in unique(names(genes))) {
      gene = as.character(genes[names(genes) == c])
      highlight.sector(sector.index = paste(c,gene), track.index = 1, col = ifelse(length(cell_col)==1,cell_col,cell_col[c]), text = c, text.vjust = text.vjust, niceFacing = TRUE,lwd=1)
    }
  }
  circos.clear()
}

#colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);")
#R语言颜色相关函数https://www.jianshu.com/p/bc94d1d6a7a0
#R语言的一些配色的R包https://zhuanlan.zhihu.com/p/358189206
my_pal <- c(ggsci::pal_d3(palette = "category20",alpha = 1)(20),ggsci::pal_d3(palette = "category20c",alpha = 1)(7))
scales::show_col(my_pal)

pal_chr <- my_pal %>% paste('"',. ,'"') %>% paste(.,collapse = ',') 
node_chr <- nodet$name %>% paste('"',. ,'"') %>% paste(.,collapse = ',')

my_color <- paste0('d3.scaleOrdinal(d3.schemeCategory20) .domain(["cell_1", "cell_2",',node_chr,']) .range(["#FF7F0ECD", "#2CA02CCD",',pal_chr,']);')
my_color

adj <- linkt %>% select('source','target','value') %>% spread(key = target,value = value) 
rownames(adj) <- adj$source
adj <- adj[,2:ncol(adj)] %>% as.matrix()

adj[is.na(adj)] <- 0

circos.clear()
circlize::chordDiagram(adj, transparency = 0.5)
