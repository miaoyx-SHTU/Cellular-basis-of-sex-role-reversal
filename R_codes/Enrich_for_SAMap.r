library("biomaRt")
library(data.table)
library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)
#library(org.Mm.eg.db)
library(msigdbr)
library(ggplot2)

hsGenes <- fread('hc_S_TGC(cxcl14)_human.csv',col.names = c('index','gene'),header=T)$gene

musGenes <- fread('hc_S_TGC(cxcl14)_mouse.csv',col.names = c('index','gene'),header=T)$gene

#https://www.soinside.com/question/hYHMikTFfKdgzc5GVyfeyP
#biomaRt has encountered an unexpected server error.
#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
human =useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "www")
#mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mouse =useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "www")

hs_musGenes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                 values = musGenes, 
                 mart = mouse, 
                 attributesL = c("hgnc_symbol"), 
                 martL = human, uniqueRows=T)

hs_musGenes = toupper(musGenes)

gene <- genes %>% bitr( fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")

go <- enrichGO(gene     = gene$ENTREZID,
              OrgDb        = org.Hs.eg.db,
              keyType       = 'ENTREZID',
              ont          = "BP",
              minGSSize    = 0,
              #maxGSSize    = 500,
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.1,
              #eps = eps,
              #verbose      = FALSE
              )

kk <- enrichKEGG(gene         = gene$ENTREZID,
                 #keyType      = 'uniprot',
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

saveRDS(go,'hc_S_TGC(cxcl14)_GO.rds')
saveRDS(kk,'hc_S_TGC(cxcl14)_KEGG.rds')

result <- go@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p

ggplot2::ggsave(p,filename = 'hc_S_TGC(cxcl14)_GO_top20.pdf',width = 10,height = 10)

result <- kk@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p

ggplot2::ggsave(p,filename = 'hc_S_TGC(cxcl14)_KEGG_top20.pdf',width = 10,height = 10)

hsGenes <- fread('hc_S_TGC(muc5)_human.csv',col.names = c('index','gene'),header=T)$gene
musGenes <- fread('hc_S_TGC(muc5)_mouse.csv',col.names = c('index','gene'),header=T)$gene
hs_musGenes = toupper(musGenes)

genes <- union(hsGenes,hs_musGenes)
length(genes)
gene <- genes %>% bitr( fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")

go <- enrichGO(gene     = gene$ENTREZID,
              OrgDb        = org.Hs.eg.db,
              keyType       = 'ENTREZID',
              ont          = "BP",
              minGSSize    = 0,
              #maxGSSize    = 500,
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.1,
              #eps = eps,
              #verbose      = FALSE
              )

kk <- enrichKEGG(gene         = gene$ENTREZID,
                 #keyType      = 'uniprot',
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

saveRDS(go,'hc_S_TGC(muc5)_GO.rds')
saveRDS(kk,'hc_S_TGC(muc5)_KEGG.rds')

result <- go@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p1 <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p1

result <- kk@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p2 <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p2

ggplot2::ggsave(p1,filename = 'hc_S_TGC(muc5)_GO_top20.pdf',width = 10,height = 10)
ggplot2::ggsave(p2,filename = 'hc_S_TGC(muc5)_KEGG_top20.pdf',width = 10,height = 10)



hsGenes <- fread('hc_S_TGC(nucb2)_human.csv',col.names = c('index','gene'),header=T)$gene
#musGenes <- fread('hc_S_TGC(nucb2)_mouse.csv',col.names = c('index','gene'),header=T)$gene
#hs_musGenes = toupper(musGenes)

#genes <- union(hsGenes,hs_musGenes)
length(hsGenes)
gene <- hsGenes %>% bitr( fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")

go <- enrichGO(gene     = gene$ENTREZID,
              OrgDb        = org.Hs.eg.db,
              keyType       = 'ENTREZID',
              ont          = "BP",
              minGSSize    = 0,
              #maxGSSize    = 500,
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.1,
              #eps = eps,
              #verbose      = FALSE
              )

kk <- enrichKEGG(gene         = gene$ENTREZID,
                 #keyType      = 'uniprot',
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

saveRDS(go,'hc_S_TGC(nucb2)_GO.rds')
saveRDS(kk,'hc_S_TGC(nucb2)_KEGG.rds')

result <- go@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p1 <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p1

result <- kk@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p2 <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p2

ggplot2::ggsave(p1,filename = 'hc_S_TGC(nucb2)_GO_top20.pdf',width = 10,height = 10)
ggplot2::ggsave(p2,filename = 'hc_S_TGC(nucb2)_KEGG_top20.pdf',width = 10,height = 10)



hsGenes <- fread('hc_S_decidual_human.csv',col.names = c('index','gene'),header=T)$gene
musGenes <- fread('hc_S_decidual_mouse.csv',col.names = c('index','gene'),header=T)$gene
hs_musGenes = toupper(musGenes)

genes <- union(hsGenes,hs_musGenes)
length(genes)
gene <- genes %>% bitr( fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")

go <- enrichGO(gene     = gene$ENTREZID,
              OrgDb        = org.Hs.eg.db,
              keyType       = 'ENTREZID',
              ont          = "BP",
              minGSSize    = 0,
              #maxGSSize    = 500,
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.1,
              #eps = eps,
              #verbose      = FALSE
              )

kk <- enrichKEGG(gene         = gene$ENTREZID,
                 #keyType      = 'uniprot',
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

saveRDS(go,'hc_S_decidual_GO.rds')
saveRDS(kk,'hc_S_decidual_KEGG.rds')

result <- go@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p1 <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p1

result <- kk@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p2 <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p2

ggplot2::ggsave(p1,filename = 'hc_S_decidual_GO_top20.pdf',width = 10,height = 10)
ggplot2::ggsave(p2,filename = 'hc_S_decidual_KEGG_top20.pdf',width = 10,height = 10)



hsGenes <- fread('hc_S_immune_human.csv',col.names = c('index','gene'),header=T)$gene
musGenes <- fread('hc_S_immune_mouse.csv',col.names = c('index','gene'),header=T)$gene
hs_musGenes = toupper(musGenes)

genes <- union(hsGenes,hs_musGenes)
length(genes)
gene <- genes %>% bitr( fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")

go <- enrichGO(gene     = gene$ENTREZID,
              OrgDb        = org.Hs.eg.db,
              keyType       = 'ENTREZID',
              ont          = "BP",
              minGSSize    = 0,
              #maxGSSize    = 500,
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.1,
              #eps = eps,
              #verbose      = FALSE
              )

kk <- enrichKEGG(gene         = gene$ENTREZID,
                 #keyType      = 'uniprot',
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

saveRDS(go,'hc_S_immune_GO.rds')
saveRDS(kk,'hc_S_immune_KEGG.rds')

result <- go@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p1 <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p1

result <- kk@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p2 <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p2

ggplot2::ggsave(p1,filename = 'hc_S_immune_GO_top20.pdf',width = 10,height = 10)
ggplot2::ggsave(p2,filename = 'hc_S_immune_KEGG_top20.pdf',width = 10,height = 10)



hsGenes <- fread('hc_S_FB_human.csv',col.names = c('index','gene'),header=T)$gene
musGenes <- fread('hc_S_FB_mouse.csv',col.names = c('index','gene'),header=T)$gene
hs_musGenes = toupper(musGenes)

genes <- union(hsGenes,hs_musGenes)
length(genes)
gene <- genes %>% bitr( fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")

go <- enrichGO(gene     = gene$ENTREZID,
              OrgDb        = org.Hs.eg.db,
              keyType       = 'ENTREZID',
              ont          = "BP",
              minGSSize    = 0,
              #maxGSSize    = 500,
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.1,
              #eps = eps,
              #verbose      = FALSE
              )

kk <- enrichKEGG(gene         = gene$ENTREZID,
                 #keyType      = 'uniprot',
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

saveRDS(go,'hc_S_FB_GO.rds')
saveRDS(kk,'hc_S_FB_KEGG.rds')

result <- go@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p1 <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p1

result <- kk@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p2 <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p2

ggplot2::ggsave(p1,filename = 'hc_S_FB_GO_top20.pdf',width = 10,height = 10)
ggplot2::ggsave(p2,filename = 'hc_S_FB_KEGG_top20.pdf',width = 10,height = 10)



hsGenes <- fread('hc_S_endo_human.csv',col.names = c('index','gene'),header=T)$gene
musGenes <- fread('hc_S_endo_mouse.csv',col.names = c('index','gene'),header=T)$gene
hs_musGenes = toupper(musGenes)

genes <- union(hsGenes,hs_musGenes)
length(genes)
gene <- genes %>% bitr( fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")

go <- enrichGO(gene     = gene$ENTREZID,
              OrgDb        = org.Hs.eg.db,
              keyType       = 'ENTREZID',
              ont          = "BP",
              minGSSize    = 0,
              #maxGSSize    = 500,
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.1,
              #eps = eps,
              #verbose      = FALSE
              )

kk <- enrichKEGG(gene         = gene$ENTREZID,
                 #keyType      = 'uniprot',
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

saveRDS(go,'hc_S_endo_GO.rds')
saveRDS(kk,'hc_S_endo_KEGG.rds')

result <- go@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p1 <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p1

result <- kk@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p2 <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p2

ggplot2::ggsave(p1,filename = 'hc_S_endo_GO_top20.pdf',width = 10,height = 10)
ggplot2::ggsave(p2,filename = 'hc_S_endo_KEGG_top20.pdf',width = 10,height = 10)



hsGenes <- fread('hc_S_basal&progenitor_human.csv',col.names = c('index','gene'),header=T)$gene
musGenes <- fread('hc_S_basal&progenitor_mouse.csv',col.names = c('index','gene'),header=T)$gene
hs_musGenes = toupper(musGenes)

genes <- union(hsGenes,hs_musGenes)
length(genes)
gene <- genes %>% bitr( fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")

go <- enrichGO(gene     = gene$ENTREZID,
              OrgDb        = org.Hs.eg.db,
              keyType       = 'ENTREZID',
              ont          = "BP",
              minGSSize    = 0,
              #maxGSSize    = 500,
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.1,
              #eps = eps,
              #verbose      = FALSE
              )

kk <- enrichKEGG(gene         = gene$ENTREZID,
                 #keyType      = 'uniprot',
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

saveRDS(go,'hc_S_basal&progenitor_GO.rds')
saveRDS(kk,'hc_S_basal&progenitor_KEGG.rds')

result <- go@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p1 <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p1

result <- kk@result[1:20,]
result$logp <- -log10(result$p.adjust)
result$Description <- factor(result$Description,levels = result$Description[order(result$logp,decreasing = F)])

p2 <- ggplot(data = result ,aes(x=Description, y=logp)) +
    geom_bar(aes(fill=Count), stat="identity", alpha=1, width=.8) +
    geom_text(aes(label = Description, hjust = 0,y=0),colour = "black",size=5)+
    scale_fill_gradient(low='#259CA2BB',high = '#A1C7BB')+
    #geom_label(aes(label=Description,y=3))+
    coord_flip() +
    xlab("Biological processes") +
    ylab('-log10 (adjusted p-value)')+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=24),
          axis.line.y = element_blank(),
          axis.ticks.y= element_blank())
p2

ggplot2::ggsave(p1,filename = 'hc_S_basal&progenitor_GO_top20.pdf',width = 10,height = 10)
ggplot2::ggsave(p2,filename = 'hc_S_basal&progenitor_KEGG_top20.pdf',width = 10,height = 10)


