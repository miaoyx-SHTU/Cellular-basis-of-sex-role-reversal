options(stringsAsFactors = F)
library("BSgenome")
library("Seahorse.bs.genome")
library("ArchR")
library("stringr")

#########################################################################
################################ clustering and annotation
options(ArchR.genome = "Seahorse.bs.genome")
addArchRChrPrefix(chrPrefix = FALSE)
addArchRThreads(threads = 1) 
fofn <- read.table("fofn")[,1]
chromSize <- read.table("Seahorse.genome.size",sep="\t",header=F)
chromSize <- chromSize[which(chromSize[,1] %in% fofn),]

gff <- read.table("seahorse.gff.txt",sep='\t',heade=F)
tx  <- gff[which(gff$V3=="mRNA"),]
tx[,9] <- gsub(";","",str_split(tx[,9],"=",simplify=T)[,3])
filter.list <-apply(tx,1,function(x)
{
if(x[7]=="+")
{
s = as.numeric(x[4])
} else {
s = as.numeric(x[5])
}
size = as.numeric(chromSize[which(chromSize[,1]==x[1]),2])
if(s<=2000 | size-s <= 2000)
{
return(x[9])
} else{
return("0")
}
})
filter.list <- filter.list[which(filter.list!="0")]

index <- c()
for(i in filter.list)
{
print(i)
index <- c(index,grep(i,gff[,9]))

}

gff <- gff[-unique(index),]
chromSize <- chromSize[which(chromSize[,1] %in% gff[,1]),]
chromSize.gr <- GRanges(chromSize[,1],IRanges(start=1,end=chromSize[,2]))
genomeAnnotation <- createGenomeAnnotation(genome = Seahorse.bs.genome,
                                           chromSizes = chromSize.gr
)

gene <- gff[which(gff$V3=="gene"),]
gene.gr <- GRanges(seqnames=gene[,1],
                   ranges=IRanges(gene[,4],gene[,5]),
                   strand = gene[,7],
                   gene_id=gsub("ID=","",gsub(";","",gene[,9])),
                   symbol=gsub("ID=","",gsub(";","",gene[,9]))
                   )
tx  <- gff[which(gff$V3=="mRNA"),]
tx_F <- tx[which(tx$V7=="+"),]
tx_R <- tx[which(tx$V7=="-"),]
temp <- str_split(c(tx_F[,9],tx_R[,9]),";",simplify=T)
tx.gr <- GRanges(seqnames=c(tx_F[,1],tx_R[,1]),
                 ranges=IRanges(c(tx_F[,4],tx_R[,5])),
                 strand = c(tx_F[,7],tx_R[,7]),
                 tx_id=1:length(c(tx_F[,7],tx_R[,7])),
                 symbol=gsub("Parent=","",temp[,2])
)



exon <- gff[which(gff$V3=="exon"),]
temp <- str_split(exon[,9],"\\.",simplify=T)
exon.gr <- GRanges(seqnames=exon[,1],
                   ranges=IRanges(exon[,4],exon[,5]),
                   strand = exon[,7],
                   gene_id=gsub("Parent=","",temp[,1]),
                   symbol=gsub("Parent=","",temp[,1]) )

geneAnnotation <- createGeneAnnotation(
  genome = "Seahorse.bs.genome",
  genes = gene.gr,
  TSS = tx.gr, 
  exons = exon.gr
)
inputFiles <- c(        "./1_cellranger_count/BP1/BP1/outs/fragments.tsv.gz",
				        "./1_cellranger_count/BP2/BP2/outs/fragments.tsv.gz",
						"./1_cellranger_count/HY1/HY1/outs/fragments.tsv.gz",
						"./1_cellranger_count/HY2/HY2/outs/fragments.tsv.gz"
						)
				
				
names(inputFiles) <- c("BP1","BP2","HY1","HY2")


ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 5, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  TileMatParams = list(tileSize = 5000)
  
)

#####################################
###########doubScores
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)



#######################
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Seahorse",
  copyArrows = TRUE, #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)


projHeme2 <- filterDoublets(projHeme1,filterRatio = 3)

ZQBP <- readRDS("./Seahorse_with_peak_motif/Save-ArchR-Project.rds")
projHeme2_filter <- projHeme2[which(projHeme2@cellColData$Sample %in% c("HY1","HY2") | projHeme2@cellColData@rownames  %in% ZQBP@cellColData@rownames)]



projHeme3 <- addIterativeLSI(
  ArchRProj = projHeme2_filter,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 4, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 50000, 
  dimsToUse = 1:30,
  force = TRUE
)

projHeme3 <- addHarmony(
  ArchRProj = projHeme3,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)


projHeme3 <- addClusters(
  input = projHeme3,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force = TRUE
)

projHeme3 <- addUMAP(
  ArchRProj = projHeme3, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)


saveArchRProject(ArchRProj = projHeme3, outputDirectory = "Seahorse", load = FALSE)

