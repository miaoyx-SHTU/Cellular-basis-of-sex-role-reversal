library(Seurat)#main analyze package
library(harmony)# use to remove batch effect
library(dplyr)
library(data.table)

#load 10x
rawdatalf1 <- Read10X(data.dir = "F:/singlecell/lesson8/BP1/")
rawdatalf2 <- Read10X(data.dir = "F:/singlecell/lesson8/BP2/")
# Create Seurat Object from matrix 
lf1 <- CreateSeuratObject(rawdatalf1,
                          min.cells = 1,
                          min.features = 0,
                          project = '10x_lf1')

lf2 <- CreateSeuratObject(rawdatalf2,
                          min.cells = 1,
                          min.features = 0,
                          project = '10x_lf2')

# add metadata original identity
lf1@meta.data$orig.ident <- "lf1"
lf2@meta.data$orig.ident <- "lf2"

#merge Seurat Object in different batch
lf <- merge(lf1,
            y = lf2,
            add.cell.ids = c("lf1.","lf2."),
            project = "10x_lf")
saveRDS(bp,file = "./pipline_reanalyze/raw_bp_rna.rds")

#load 10x
rawdataef1 <- Read10X(data.dir = "F:/singlecell/lesson8/ZQ1/")

rawdataef2 <- Read10X(data.dir = "F:/singlecell/lesson8/ZQ2/")

# Create Seurat Object from matrix 
ef1 <- CreateSeuratObject(rawdataef1,
                          min.cells = 1,
                          min.features = 0,
                          project = '10x_ef1')

ef2 <- CreateSeuratObject(rawdataef2,
                          min.cells = 1,
                          min.features = 0,
                          project = '10x_ef2')

# add metadata original identity
ef1@meta.data$orig.ident <- "ef1"
ef2@meta.data$orig.ident <- "ef2"

#merge Seurat Object in different batch
ef <- merge(ef1,
            y = ef2,
            add.cell.ids = c("ef1.","ef2."),
            project = "10x_ef")
saveRDS(zq,file = "./pipline_reanalyze/raw_zq_rna.rds")

#load 10x
rawdataep1 <- Read10X(data.dir = "F:/singlecell/lesson8/HY1/")

rawdataep2 <- Read10X(data.dir = "F:/singlecell/lesson8/HY2/")

# Create Seurat Object from matrix 
ep1 <- CreateSeuratObject(rawdataep1,
                          min.cells = 1,
                          min.features = 0,
                          project = '10x_ep1')

ep2 <- CreateSeuratObject(rawdataep2,
                          min.cells = 1,
                          min.features = 0,
                          project = '10x_ep2')

# add metadata original identity
ep1@meta.data$orig.ident <- "ep1"
ep2@meta.data$orig.ident <- "ep2"

#merge Seurat Object in different batch
ep <- merge(ep1,
            y = ep2,
            add.cell.ids = c("ep1.","ep2."),
            project = "10x_ep")
saveRDS(hy,file = "./pipline_reanalyze/raw_hy_rna.rds")

#load 10x(bf)
rawdatabf <- Read10X(data.dir = "F:/singlecell/lesson8/new_seahorse/Male_filtered_feature_bc_matrix/")

# Create Seurat Object from matrix 
bf <- CreateSeuratObject(rawdatabf,
                         min.cells = 1,
                         min.features = 0,
                         project = '10x_bf')
# add metadata original identity
bf@meta.data$orig.ident <- "bf"

#load 10x(mf)
rawdatamf <- Read10X(data.dir = "F:/singlecell/lesson8/new_seahorse/3m_filtered_feature_bc_matrix/")

# Create Seurat Object from matrix 
mf <- CreateSeuratObject(rawdatamf,
                         min.cells = 1,
                         min.features = 0,
                         project = '10x_mf')
# add metadata original identity
mf@meta.data$orig.ident <- "mf"

#load 10x(mp)
rawdatamp <- Read10X(data.dir = "F:/singlecell/lesson8/new_seahorse/p2_filtered_feature_bc_matrix/")

# Create Seurat Object from matrix 
mp <- CreateSeuratObject(rawdatamp,
                         min.cells = 1,
                         min.features = 0,
                         project = '10x_mp')
# add metadata original identity
mp@meta.data$orig.ident <- "mp"

#load 10x(lp)
rawdatalp <- Read10X(data.dir = "F:/singlecell/lesson8/new_seahorse/p3_filtered_feature_bc_matrix/")

# Create Seurat Object from matrix 
lp <- CreateSeuratObject(rawdatalp,
                         min.cells = 1,
                         min.features = 0,
                         project = '10x_lp')
# add metadata original identity
lp@meta.data$orig.ident <- "lp"
# add original identity in to metadata
lf1 <- lf[ ,lf@meta.data$orig.ident == 'lf1']
lf2 <- lf[ ,lf@meta.data$orig.ident == 'lf2']

ep1 <- ep[ ,ep@meta.data$orig.ident == 'ep1']
ep2 <- ep[ ,ep@meta.data$orig.ident == 'ep2']

ef1 <- ef[ ,ef@meta.data$orig.ident == 'ef1']
ef2 <- ef[ ,ef@meta.data$orig.ident == 'ef2']

bf <- bf[ ,bf@meta.data$orig.ident == 'bf']
mf <- mf[ ,mf@meta.data$orig.ident == 'mf']
mp <- mp[ ,mp@meta.data$orig.ident == 'mp']
lp <- lp[ ,lp@meta.data$orig.ident == 'lp']

#add group information into metadata
lf1@meta.data$group <- 'lf'
lf2@meta.data$group <- 'lf'
ep1@meta.data$group <- 'ep'
ep2@meta.data$group <- 'ep'
ef1@meta.data$group <- 'ef'
ef2@meta.data$group <- 'ef'

bf@meta.data$group <- 'bf'
mf@meta.data$group <- 'mf'
mp@meta.data$group <- 'mp'
lp@meta.data$group <- 'lp'

# merge Seurat Object
allRNA  <- merge(lf1,
                 y = c(lf2, ep1, ep2, ef1, ef2, bf, mf, mp, lp),
                 #add.cell.ids = c("BP.","HY.",'ZQ.'),
                 project = "10x_RNA")
#add batch information into metadata
allRNA@meta.data$batch = allRNA@meta.data$orig.ident
#show statistics of metadata
table(allRNA$group)
table(allRNA$orig.ident)
table(allRNA$batch)

#save in to rds
allRNA
saveRDS(allRNA, ,file = "./pipline_reanalyze/raw_all_rna_harmony.rds")

# general steps of scRNA-seq analyze in Seurat
allRNA <- NormalizeData(allRNA, verbose = FALSE)

allRNA <- FindVariableFeatures(object = allRNA, mean.cutoff = c(0.125,8), dispersion.cutoff = c(1,Inf))

allRNA <- ScaleData(allRNA, verbose = FALSE)

allRNA <- RunPCA(allRNA, verbose = FALSE, features = VariableFeatures(object = allRNA))

ElbowPlot(allRNA)

allRNA# just show Seurat Object in jupyter notebook(will not work in python script)

# remove batch effect use Harmony
allRNA <- RunHarmony(allRNA,
                     "batch",
                     plot_convergence = TRUE,
                     #assay.use = "SCT", # use it when use SCTranform before
                     assay.use = "RNA"
)

#use harmony in below analyze
allRNA <- RunTSNE(allRNA, reduction = "harmony", verbose = FALSE, dims = 1:20)
allRNA <- RunUMAP(allRNA, reduction = "harmony", verbose = FALSE, dims = 1:20)

allRNA <- FindNeighbors(allRNA, reduction = "harmony", dims = 1:20)
allRNA <- FindClusters(object = allRNA, resolution = 0.125)

allRNA

#draw T-SNE plot
DimPlot(object =allRNA,
        reduction = "tsne",
        label = T,
)
DimPlot(object =allRNA,
        group.by = 'group',
        reduction = "tsne")
DimPlot(object =allRNA,
        group.by = 'batch',
        reduction = "tsne")

#draw UMAP plot
DimPlot(object =allRNA,
        reduction = "umap",
        label = T,
)
DimPlot(object =allRNA,
        group.by = 'group',
        reduction = "umap")
DimPlot(object =allRNA,
        group.by = 'batch',
        reduction = "umap")

saveRDS(allRNA ,file = "./pipline_reanalyze/all_rna_harmony.rds")

#find markers for each cluster
all_markers <- FindAllMarkers(allRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(all_markers ,file = "./pipline_reanalyze/all_rna_markers_harmony.rds") 


FeaturePlot(allRNA, features = 'EVM0012341', 
            reduction = "umap", 
            label=T)
#remove cls0(blood cells)
#remove blood cells
filtered_lf1 <- CreateSeuratObject(allRNA@assays$RNA[  ,allRNA@meta.data$seurat_clusters!=1 & allRNA@meta.data$orig.ident == 'lf1'],
                                   min.cells = 1,
                                   min.features = 0,
                                   project = '10x_filtered_lf1')
filtered_lf2 <- CreateSeuratObject(allRNA@assays$RNA[  ,allRNA@meta.data$seurat_clusters!=1 & allRNA@meta.data$orig.ident == 'lf2'],
                                   min.cells = 1,
                                   min.features = 0,
                                   project = '10x_filtered_lf2')

filtered_ep1 <- CreateSeuratObject(allRNA@assays$RNA[  ,allRNA@meta.data$seurat_clusters!=1 & allRNA@meta.data$orig.ident == 'ep1'],
                                   min.cells = 1,
                                   min.features = 0,
                                   project = '10x_filtered_1')
filtered_ep2 <- CreateSeuratObject(allRNA@assays$RNA[  ,allRNA@meta.data$seurat_clusters!=1 & allRNA@meta.data$orig.ident == 'ep2'],
                                   min.cells = 1,
                                   min.features = 0,
                                   project = '10x_filtered_hy2')

filtered_ef1 <- CreateSeuratObject(allRNA@assays$RNA[  ,allRNA@meta.data$seurat_clusters!=1 & allRNA@meta.data$orig.ident == 'ef1'],
                                   min.cells = 1,
                                   min.features = 0,
                                   project = '10x_filtered_zq1')
filtered_ef2 <- CreateSeuratObject(allRNA@assays$RNA[  ,allRNA@meta.data$seurat_clusters!=1 & allRNA@meta.data$orig.ident == 'ef2'],
                                   min.cells = 1,
                                   min.features = 0,
                                   project = '10x_filtered_zq2')
filtered_bf <- CreateSeuratObject(allRNA@assays$RNA[  ,allRNA@meta.data$seurat_clusters!=1 & allRNA@meta.data$orig.ident == 'bf'],
                                   min.cells = 1,
                                   min.features = 0,
                                   project = '10x_filtered_bf')
filtered_mf <- CreateSeuratObject(allRNA@assays$RNA[  ,allRNA@meta.data$seurat_clusters!=1 & allRNA@meta.data$orig.ident == 'mf'],
                                  min.cells = 1,
                                  min.features = 0,
                                  project = '10x_filtered_mf')
filtered_mp <- CreateSeuratObject(allRNA@assays$RNA[  ,allRNA@meta.data$seurat_clusters!=1 & allRNA@meta.data$orig.ident == 'mp'],
                                  min.cells = 1,
                                  min.features = 0,
                                  project = '10x_filtered_mp')
filtered_lp <- CreateSeuratObject(allRNA@assays$RNA[  ,allRNA@meta.data$seurat_clusters!=1 & allRNA@meta.data$orig.ident == 'lp'],
                                  min.cells = 1,
                                  min.features = 0,
                                  project = '10x_filtered_lp')

filtered_lf1@meta.data$orig.ident <- "lf1"
filtered_lf2@meta.data$orig.ident <- "lf2"
filtered_ep1@meta.data$orig.ident <- "ep1"
filtered_ep2@meta.data$orig.ident <- "ep2"
filtered_ef1@meta.data$orig.ident <- "ef1"
filtered_ef2@meta.data$orig.ident <- "ef2"
filtered_bf@meta.data$orig.ident <- "bf"
filtered_mf@meta.data$orig.ident <- "mf"
filtered_mp@meta.data$orig.ident <- "mp"
filtered_lp@meta.data$orig.ident <- "lp"

filtered_lf1@meta.data$group <- "lf"
filtered_lf2@meta.data$group <- "lf"
filtered_ep1@meta.data$group <- "ep"
filtered_ep2@meta.data$group <- "ep"
filtered_ef1@meta.data$group <- "ef"
filtered_ef2@meta.data$group <- "ef"
filtered_bf@meta.data$group <- "bf"
filtered_mf@meta.data$group <- "mf"
filtered_mp@meta.data$group <- "mp"
filtered_lp@meta.data$group <- "lp"

filtered_lf1@meta.data$batch =
filtered_hy2@meta.data$batch <- "B2"
filtered_bp1@meta.data$batch <- "B1"
filtered_bp2@meta.data$batch <- "B2"
filtered_zq1@meta.data$batch <- "B1"
filtered_zq2@meta.data$batch <- "B2"

filtered <- merge(filtered_lf1,
                  y = c(filtered_lf2,filtered_ep1,filtered_ep2,filtered_ef1,filtered_ef2,filtered_bf, filtered_mf,filtered_mp, filtered_lp),
                  project = "10x_rna_filtered")
filtered@meta.data$batch=filtered@meta.data$orig.ident
table(filtered$group)
table(filtered$orig.ident)
table(filtered$batch)

filtered <- NormalizeData(filtered, verbose = FALSE)

filtered <- FindVariableFeatures(object = filtered, mean.cutoff = c(0.125,8), dispersion.cutoff = c(1,Inf))

filtered <- ScaleData(filtered, verbose = FALSE)

filtered <- RunPCA(filtered, verbose = FALSE, features = VariableFeatures(object = filtered))

ElbowPlot(filtered)
filtered

filtered <- RunHarmony(filtered,
                       "batch",
                       plot_convergence = TRUE,
                       #assay.use = "SCT",
                       assay.use = "RNA"
)

filtered <- RunTSNE(filtered, reduction = "harmony", verbose = FALSE, dims = 1:18)
filtered <- RunUMAP(filtered, reduction = "harmony", verbose = FALSE, dims = 1:18)

filtered <- FindNeighbors(filtered, reduction = "harmony", dims = 1:18)
filtered <- FindClusters(object = filtered, resolution = 0.125)