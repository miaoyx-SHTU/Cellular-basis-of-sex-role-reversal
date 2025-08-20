options(stringsAsFactors = F)
library("BSgenome")
library("Seahorse.bs.genome")
library("ArchR")
library("stringr")

#########################################################################
#################### clustering and annotation


options(ArchR.genome = "Seahorse.bs.genome")
addArchRChrPrefix(chrPrefix = FALSE)
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
inputFiles <- c("ZQ1_fragments.tsv.gz",
				        "ZQ2_fragments.tsv.gz",
				        "BP1_fragments.tsv.gz",
				        "BP2_fragments.tsv.gz")
				
				
names(inputFiles) <- c("ZQ1","ZQ2","BP1","BP2")


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

saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Seahorse", load = FALSE)


projHeme2 <- filterDoublets(projHeme1,filterRatio = 3)

projHeme2 <- addIterativeLSI(
  ArchRProj = projHeme2,
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

projHeme2 <- addHarmony(
  ArchRProj = projHeme2,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)


projHeme2 <- addClusters(
  input = projHeme2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force = TRUE
)

projHeme2 <- addUMAP(
  ArchRProj = projHeme2, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")

plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "DoubletScore", embedding = "UMAP")
plotEmbedding(ArchRProj = projHeme3, colorBy = "cellColData", name = "predictedScore_Un", embedding = "UMAP")

markersGS <- getMarkerFeatures(
  ArchRProj = projHeme2,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

ATAC.marker <- c()
for(i in names(markerList))
{
  print(i)
  ATAC.marker <- c(ATAC.marker,markerList[[i]]$name[1:150])
}

GeneScoreMatrix <- getMatrixFromProject(
  ArchRProj = projHeme2,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

gene.mat <- GeneScoreMatrix@assays@data@listData$GeneScoreMatrix
rownames(gene.mat) <- GeneScoreMatrix@elementMetadata@listData$name
colnames(gene.mat) <- GeneScoreMatrix@colData@rownames
var <- FindVariableFeatures(gene.mat,nfeatures =2000)
ATAC.top.2000 <- rownames(var)[order(var$vst.variance.standardized,decreasing = T)[1:2000]]


library(Seurat)

RNA <- readRDS("all_3_stage_delbf_all_filtered_blood.rds")
RNA2 <- RenameIdents(object = RNA,
                    "0" = "Epi(krt8+)",
                    "1" = "Epi(tfa+)",
                    "2" = "fibroblast",
                    "3" = "Epi(sox4+)",
                    "4" = "endothelial",
                    "5" = "ionocyteI",
                    "6" = "ionocyteII",
                    "7" = "immune cell",
                    "8" = "Epi(gata1a+)",
                    "9" = "basal cell",
                    "10" = "smooth muscle cells",
                    "11" = "mucus cell",
                    "12" = "Epi(krt8+)",
                    "13" = "ionocyte III",
                    "14" = "Epi(pastn3+)",
                    "15" = "Epi(notch3+)",
                    "16" = "Epi(notch3+)",
                    "17" = "Epi(sox4+)")

markers <- FindAllMarkers(RNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
marker.use <- c()

for(i in 0:17)
{
  cluster.marker <- rownames(markers)[which(markers$cluster==i)]
  marker.use <- c(marker.use,cluster.marker[1:100])
}
marker.use <- unique(marker.use)

index <- marker.use
index <- index[which(index %in% RNA@assays$RNA@counts@Dimnames[[1]] )]

marker.use <- unique(c(marker.use,ATAC.marker))
marker.use2 <- unique(c(VariableFeatures(RNA),ATAC.marker))



projHeme3 <- addGeneIntegrationMatrix(
  ArchRProj = projHeme2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = RNA,
  addToArrow = TRUE,
  groupRNA = "seurat_clusters",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
  #genesUse = VariableFeatures(RNA) ,               
  genesUse = marker.use2,
  force = TRUE
  )

saveArchRProject(ArchRProj = projHeme3, outputDirectory = "Seahorse.projHeme3", load = FALSE)

p1 <- plotEmbedding(ArchRProj = projHeme3, colorBy = "cellColData", name = "DoubletScore", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme3, colorBy = "cellColData", name = "predictedScore_Un", embedding = "UMAP")
ggAlignPlots(p1, p2 ,  type = "h")




cM <- as.matrix(confusionMatrix(projHeme3$Clusters, projHeme3$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
Assignments <- cbind(preClust, rownames(cM)) #Assignments

projHeme3_filter <- projHeme3[-which(projHeme3@cellColData$Clusters %in% c("C9","C22","C2","C13"))] 

plotEmbedding(ArchRProj = projHeme3_filter, colorBy = "cellColData", name = "DoubletEnrichment", embedding = "UMAP")
plotEmbedding(ArchRProj = projHeme3_filter, colorBy = "cellColData", name = "predictedScore_Un", embedding = "UMAP")

p1 <- plotEmbedding(ArchRProj = projHeme3_filter, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme3_filter, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  projHeme3_filter, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un", 
  size=-0.05,plotAs="points"
)

ggAlignPlots(p1, p2 , p3, type = "h")


projHeme3_filter@cellColData$Annotation_manual <- "unKnown"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C3","C4","C5","C6","C7"))] <- "Epi(krt8+)"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C1"))] <- "Epi(notch3+)"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C8"))] <- "endothelial"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C10","C11"))] <- "Epi(pastn3+)"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C12","C14"))] <- "Epi(tfa+)"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C15"))] <- "ionocyte II"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C16"))] <- "ionocyte I"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C17"))] <- "muscle cell"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C18","C19"))] <- "fibroblast"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C20") & projHeme3_filter@cellColData$predictedGroup_Un == 10)] <- "smooth muscle"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C20") & projHeme3_filter@cellColData$predictedGroup_Un == 13)] <- "ionocyte III"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C23"))] <- "Epi(gata1a+)"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C24","C25"))] <- "immune cell"
projHeme3_filter <- projHeme3_filter[which(projHeme3_filter@cellColData$Annotation_manual != "unKnown")]

plotEmbedding(
  projHeme3_filter, 
  colorBy = "cellColData", 
  name = "Annotation_manual", 
  size=-0.05,plotAs="points"
)


pal <- paletteDiscrete(values = projHeme3_filter@cellColData$Annotation_manual)    
pal <- c("#f69775","#F092E5","#fab26b","#6CA3FF","#73cfe0",
         "#927bb8","#32c2f2","#FFA9AA","#e2c936","#00C282", 
         "#dce67f","#59bd87","#556C9A")
names(pal) <- unique(projHeme3_filter@cellColData$Annotation_manual)

p1 <- plotEmbedding(
  projHeme3_filter, 
  colorBy = "cellColData", 
  name = "Annotation_manual", 
  plotAs="points",
  pal = pal
)

print(p1)

#####################################################################
################### peak calling and motif enrichment

options(stringsAsFactors = F)
library("BSgenome")
library("Seahorse.bs.genome")
library("ArchR")
library("stringr")

addArchRThreads(threads = 1) 
#############################
#####peak calling 
projHeme3 <- readRDS("./Seahorse.projHeme3/Save-ArchR-Project.rds")
projHeme3@projectMetadata@listData$outputDirectory <- "./Seahorse.projHeme3/"
projHeme3@sampleColData@listData$ArrowFiles <- c("./Seahorse.projHeme3/ArrowFiles/BP1.arrow",
												 "./Seahorse.projHeme3/ArrowFiles/BP2.arrow",
												 "./Seahorse.projHeme3/ArrowFiles/ZQ1.arrow",
												 "./Seahorse.projHeme3/ArrowFiles/ZQ2.arrow")
projHeme3_filter <- projHeme3[-which(projHeme3@cellColData$Clusters %in% c("C9","C22","C2","C13"))] # remove low-quality annotated clusters
projHeme3_filter@cellColData$Annotation_manual <- "unKnown"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C3","C4","C5","C6","C7"))] <- "Epi(krt8+)"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C1"))] <- "Epi(notch3+)"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C8"))] <- "endothelial"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C10","C11"))] <- "Epi(pastn3+)"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C12","C14"))] <- "Epi(tfa+)"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C15"))] <- "ionocyte II"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C16"))] <- "ionocyte I"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C17"))] <- "muscle cell"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C18","C19"))] <- "fibroblast"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C20") & projHeme3_filter@cellColData$predictedGroup_Un == 10)] <- "smooth muscle"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C20") & projHeme3_filter@cellColData$predictedGroup_Un == 13)] <- "ionocyte III"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C23"))] <- "Epi(gata1a+)"
projHeme3_filter@cellColData$Annotation_manual[which(projHeme3_filter@cellColData$Clusters %in% c("C24","C25"))] <- "immune cell"
projHeme3_filter <- projHeme3_filter[which(projHeme3_filter@cellColData$Annotation_manual != "unKnown")]



projHeme4 <- addGroupCoverages(
  ArchRProj = projHeme3_filter, 
  groupBy = "Annotation_manual"								
)


projHeme4 <- addReproduciblePeakSet(
    ArchRProj = projHeme4, 
    groupBy = "Annotation_manual",
	pathToMacs2	= "/opt/service/miniconda3/envs/rex_env/bin/macs2",
	genomeAnnotation = getGenomeAnnotation(projHeme4),
	geneAnnotation = getGeneAnnotation(projHeme4),
	genomeSize = 419030441
)

# save.image(file="temp1.Rdata")

peak <- getPeakSet(projHeme4)
index <- names(peak)

names(peak) <- NULL
peak <- as.data.frame(peak)
peak$cluster_id <- index
peak <- peak[order(peak$cluster),]
write.csv(peak,file="all.cluster.peaks.csv",row.names=F,quote=F)

projHeme5 <- addPeakMatrix(projHeme4)
getAvailableMatrices(projHeme5)

#############################################
### motif enrichment

markersPeaks <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix", 
  groupBy = "Annotation_manual",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
pdf("marker.peaks.pdf")
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


# markerTest <- getMarkerFeatures(
  # ArchRProj = projHeme5, 
  # useMatrix = "PeakMatrix",
  # groupBy = "Annotation_manual",
  # testMethod = "wilcoxon",
  # bias = c("TSSEnrichment", "log10(nFrags)"),
  # useGroups = "6",
  # bgdGroups = as.character(c(0:5,7:13))
# )


projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "homer", name = "Motif",species="mus musculus", force = TRUE )


motifsUp <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projHeme5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = projHeme5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = TRUE)
pdf( "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

dev.off()
# save.image(file="temp2.Rdata")



projHeme5 <- addBgdPeaks(projHeme5)

projHeme5 <- addDeviationsMatrix(
  ArchRProj = projHeme5, 
  peakAnnotation = "Motif",
  force = TRUE
)

motifs <- c("Grhl1","Gata5") #"ARE.NR_6", AR.halfsite.NR_7
markerMotifs <- getFeatures(projHeme5, select = c("ARE.NR_6"), useMatrix = "MotifMatrix")
markerMotifs
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% "z:Terf2_648"]
markerMotifs


pdf("ARE.pdf")
p <- plotEmbedding(
    ArchRProj = projHeme5, 
    colorBy = "MotifMatrix", 
    name = markerMotifs, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme5)
)
print(p)
dev.off()


p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})

do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
dev.off()

saveArchRProject(ArchRProj = projHeme5, outputDirectory = "Seahorse_with_peak_motif", load = FALSE)
##########################################################################################################
############################# pseudotime analysis

options(stringsAsFactors = F)
library("BSgenome")
library("Seahorse.bs.genome")
library("ArchR")
library("stringr")
library("Seurat")

addArchRThreads(threads = 1) 
#############################

projHeme5 <- readRDS("./Seahorse_with_peak_motif/Save-ArchR-Project.rds")

projHeme5@cellColData$Annotation_manual[which(projHeme5@cellColData$Clusters %in% c("C3"))] <- "Epi(krt8+)3"
projHeme5@cellColData$Annotation_manual[which(projHeme5@cellColData$Clusters %in% c("C4"))] <- "Epi(krt8+)4"
projHeme5@cellColData$Annotation_manual[which(projHeme5@cellColData$Clusters %in% c("C5"))] <- "Epi(krt8+)5"
projHeme5@cellColData$Annotation_manual[which(projHeme5@cellColData$Clusters %in% c("C6"))] <- "Epi(krt8+)6"
projHeme5@cellColData$Annotation_manual[which(projHeme5@cellColData$Clusters %in% c("C7"))] <- "Epi(krt8+)7"

trajectory <- c("C7", "C1", "C10", "C11", "C12", "C14")

projHeme5 <- addTrajectory(
    ArchRProj = projHeme5, 
    name = "ZQBP", 
    groupBy = "Clusters",
    trajectory = trajectory, 
    embedding = "UMAP", 
    force = TRUE
)


p <- plotTrajectory(projHeme5, trajectory = "ZQBP", colorBy = "cellColData", name = "ZQBP",addArrow=F)
pdf("ZQ_BP.pdf")
print(p)
dev.off()



trajGSM <- getTrajectory(ArchRProj = projHeme5, name = "ZQBP", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
trajGIM <- getTrajectory(ArchRProj = projHeme5, name = "ZQBP", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)

corGSM_MM <- correlateTrajectories(trajGSM, trajGIM)
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajGIM2 <- trajGIM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined,withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajGIM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))

genes <- read.table("pse.txt",sep="\t",header=T)

for(i in 1:nrow(genes))
{
	print(i)
	rownames(trajGSM2)[grep(genes[i,1],rownames(trajGSM2))] <- genes[i,2]
	rownames(trajGIM2)[grep(genes[i,1],rownames(trajGIM2))] <- genes[i,2]
}

# col1="coolwarm"
# col2="zissou"

# col1 <- c("#98a9ba","#fff6fe","#ad6f73","#c6c0a3","#fff6fe","#615373")
# col2 <- c("#98a9ba","#a29496","#ad6f73","#c6c0a3","#a29496","#615373")

# col1 <- c("#98a9ba","#fff6fe","#ad6f73")
# col2 <- c("#c6c0a3","#fff6fe","#615373")

col1 <- c("#98a9ba","#a29496","#ad6f73")
col2 <- c("#c6c0a3","#a29496","#615373")

ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = col1,  varCutOff = 0, rowOrder = rowOrder, labelMarkers = genes[,2],labelTop=0)     #zissou sambaNight fireworks
ht2 <- plotTrajectoryHeatmap(trajGIM2,  pal = col2, varCutOff = 0, rowOrder = rowOrder, labelMarkers = genes[,2],labelTop=0)       #greenBlue whiteBlue beach

pdf(paste0("GeneScoreMatrix.GeneIntegrationMatrix.",col1,"_",col2,".pdf"),height=10,width=15)
print(ht2 + ht1)
dev.off()

library(stringr)
geneList <- rownames(combinedMat)
gene_id <- str_split(geneList,":",simplify=T)[,2]
write.table(gene_id,file="pesudotime.gene.list",quote=F,row.names=F,col.names=F)



for( i in seq(100,500,by=100))
{
print(i)

col1="coolwarm"
col2="zissou"
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = col1),  varCutOff = 0, rowOrder = rowOrder[seq(285,350,by=1)])     #zissou sambaNight fireworks
ht2 <- plotTrajectoryHeatmap(trajGIM2, pal = paletteContinuous(set = col2), varCutOff = 0, rowOrder = rowOrder[seq(285,350,by=1)])       #greenBlue whiteBlue beach

pdf(paste0("GeneScoreMatrix.GeneIntegrationMatrix.",i,".pdf"),height=10,width=15)
print(ht2 + ht1)
dev.off()
}


geneList2 <- rownames(combinedMat)[285:320]
gene_id2 <- str_split(geneList2,":",simplify=T)[,2]
write.table(gene_id2,file="pesudotime.gene2.list",quote=F,row.names=F,col.names=F)

#########################################################################################
################################ ARE motif analysis

options(stringsAsFactors = F)
library("BSgenome")
library("Seahorse.bs.genome")
library("ArchR")
library("stringr")
library("Seurat")


options(ArchR.genome = "Seahorse.bs.genome")
addArchRThreads(threads = 1) 
addArchRChrPrefix(chrPrefix = FALSE)
#############################

load("temp3.Rdata")


peak.gr <- GRanges(peak[,1],IRanges(as.numeric(peak[,2]),as.numeric(peak[,3])))
motifPositions <- getPositions(projHeme5)
ARE_location <- motifPositions[["ARE.NR_6"]]

hits <- as.data.frame(findOverlaps(peak.gr,ARE_location))
peak$ARE <- 0
temp <- table(hits[,1])
peak$ARE[as.numeric(names(temp))] <- temp

#peak_candidnate <- peak[which(peak$peakType=="Promoter" & peak$ARE >0 ),]
peak_candidnate2 <- peak[which( peak$ARE >0 ),]

#write.csv(peak_candidnate,file="AR.target.gene.csv",quote=F,row.names=F)
write.csv(peak_candidnate2,file="AR.target.gene.v2.csv",quote=F,row.names=F)


GeneScoreMatrix <- getMatrixFromProject(
  ArchRProj = projHeme5,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

Motif_matrix <- getMatrixFromProject(
  ArchRProj = projHeme5,
  useMatrix = "MotifMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

GeneExpression <- getMatrixFromProject(
  ArchRProj = projHeme5,
  useMatrix = "GeneIntegrationMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)


Motif_matrix_dataframe <- as.data.frame(Motif_matrix@assays@data$deviations)
GeneScoreMatrix_mat <- as.data.frame(GeneScoreMatrix@assays@data$GeneScoreMatrix)
rownames(GeneScoreMatrix_mat) <- GeneScoreMatrix@elementMetadata@listData$name
GeneScoreMatrix_mat <- GeneScoreMatrix_mat[which(rownames(GeneScoreMatrix_mat) %in% peak_candidnate2$nearestTSS),]
expressionMatrix <- as.data.frame(GeneExpression@assays@data@listData$GeneIntegrationMatrix)
rownames(expressionMatrix) <- GeneExpression@elementMetadata@listData$name

motif_score <- t(Motif_matrix_dataframe["ARE.NR_6",])
motif_ATAC_cor <- cor(motif_score,t(GeneScoreMatrix_mat))
motif_geneExp_cor <- cor(motif_score,t(expressionMatrix))






