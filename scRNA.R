tung_counts <- read.table("data/tung/molecules.txt", sep = "\t")
tung_annotation <- read.table("data/tung/annotation.txt", sep = "\t", header = TRUE)
# load the library
library(SingleCellExperiment)

# note that the data passed to the assay slot has to be a matrix!
tung <- SingleCellExperiment(
  assays = list(counts = as.matrix(tung_counts)),
  colData = tung_annotation
)

# remove the original tables as we don't need them anymore
rm(tung_counts, tung_annotation)


class(colData(tung))
class(assay(tung))
counts(tung)[1:2,1:3]
base::table(tung$batch)
#add new column in colData
colData(tung)$total_counts <- colSums(counts(tung))
#add new assay
assay(tung, "logcounts") <- log2(counts(tung) + 1)
assay(tung, "cpm") <- counts(tung)/tung$total_counts/1e6
# first 10 rows and 4 columns of the logcounts assay
logcounts(tung)[1:10, 1:4] # or: assay(tung, "logcounts")[1:10, 1:5]

colMeans(counts(tung))
colData(tung)$mean_counts <- colMeans(counts(tung))
# subset by numeric index

# calculate the mean counts per gene
gene_means <- rowMeans(counts(tung))
tung[gene_means > 0.01, ]
# total number of detected genes per cell
total_detected_per_cell <- colSums(counts(tung) > 0)
tung[, total_detected_per_cell > 5000]

# print the first 10 values
total_detected_per_cell[1:10]
# print the first 10 values
gene_means[1:10]
#Create a new object called tung_filtered which contains:
#  cells with at least 25000 total counts
#genes that have more than 5 counts in at least half of the cells
filterCell = colSums(counts(tung)) >=25000
k=dim(counts(tung))[2]/2
filterGene = rowSums(counts(tung)>5) > k 
tung_filtered = tung[filterGene, filterCell]
dim(counts(tung_filtered))
#------------------------------
cell_info <- as.data.frame(colData(tung))
head(cell_info)
# load the library
library(ggplot2)

ggplot(data = cell_info, aes(x = batch, y = total_counts)) +
  geom_violin(fill = 'brown') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#-----------------------------------
library(scater)
ggcells(tung, aes(x = batch, y = total_counts)) + 
  geom_violin(fill = 'orange') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

ggcells(tung, aes(x = batch, y = ENSG00000198938), exprs_values = "logcounts") + 
  geom_violin(fill = 'coral2') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

colData(tung)$var_counts = colVars(counts(tung))
ggcells(tung, aes(x = mean_counts, y = var_counts,color = batch)) + 
  geom_point(fill = 'coral2') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
#-----------------------------
library(DropletUtils)

# importing the raw count data
sce <- read10xCounts("data/pbmc_1k_raw")

# importing the pre-filtered count data
sce <- read10xCounts("data/pbmc_1k_filtered")
colMeans(counts(tung))
#------------------------------
#chp 6 QC
#----------------------------
library(scater)
library(SingleCellExperiment)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
#----------------------------
molecules <- read.delim("data/tung/molecules.txt",row.names=1)
annotation <- read.delim("data/tung/annotation.txt",stringsAsFactors = T)
head(molecules[,1:3])
head(annotation)
#Here we set altExp to contain ERCC, removing ERCC features from the main object:
umi <- SingleCellExperiment(assays = list(counts = as.matrix(molecules)), colData = annotation)
altExp(umi,"ERCC") <- umi[grep("^ERCC-",rownames(umi)), ]
umi <- umi[grep("^ERCC-",rownames(umi),invert = T), ]
ensdb_genes <- genes(EnsDb.Hsapiens.v86)
MT_names <- ensdb_genes[seqnames(ensdb_genes) == "MT"]$gene_id
is_mito <- rownames(umi) %in% MT_names
table(is_mito)
#---------------------
#QC function from scater
umi_cell <- perCellQCMetrics(umi,subsets=list(Mito=is_mito))
umi_feature <- perFeatureQCMetrics(umi)
head(umi_cell)
umi <- addPerCellQC(umi, subsets=list(Mito=is_mito))
umi <- addPerFeatureQC(umi)
hist(
  umi$total,
  breaks = 100
)
abline(v = 25000, col = "red")
dev.off()
#-------------------------
#Be careful to specify if the correct direction of the deviation: indeed, low number of detected genes,
#but high MT gene percentage, are hallmarks of a low quality cell:
qc.lib2 <- isOutlier(umi_cell$sum, log=TRUE, type="lower")
attr(qc.lib2, "thresholds")
qc.nexprs2 <- isOutlier(umi_cell$detected, log=TRUE, type="lower")
attr(qc.nexprs2, "thresholds")
qc.spike2 <- isOutlier(umi_cell$altexps_ERCC_percent, type="higher")
attr(qc.spike2, "thresholds")
qc.mito2 <- isOutlier(umi_cell$subsets_Mito_percent, type="higher")
attr(qc.mito2, "thresholds")
discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2
DataFrame(LibSize=sum(qc.lib2), NExprs=sum(qc.nexprs2), SpikeProp=sum(qc.spike2), MitoProp=sum(qc.mito2), Total=sum(discard2))
#All the actions performed above could be done in one scater command, 
reasons <- quickPerCellQC(umi_cell, sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent"))
colSums(as.matrix(reasons))
umi$discard <- reasons$discard
plotColData(umi, x="sum", y="subsets_Mito_percent", colour_by="discard")
plotColData(umi, x="sum", y="detected", colour_by="discard")
plotColData(umi, x="altexps_ERCC_percent", y="subsets_Mito_percent",colour_by="discard")
library(scales)
plotColData(umi, x="sum", y="detected", colour_by="discard", other_fields = "individual") + 
  facet_wrap(~individual) + scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3))
plotHighestExprs(umi, exprs_values = "counts", 
                 feature_names_to_plot = "SYMBOL", colour_cells_by="detected")

keep_feature <- nexprs(umi,byrow = TRUE,detection_limit = 1) >= 2
rowData(umi)$discard <- ! keep_feature
table(rowData(umi)$discard)
assay(umi, "logcounts_raw") <- log2(counts(umi) + 1)
saveRDS(umi, file = "data/tung/umi.rds")
#By default only the top 500 most variable genes are used by scater to calculate the PCA. 
#This can be adjusted by changing the ntop argument.
umi.qc <- umi[! rowData(umi)$discard,! colData(umi)$discard]
umi <- runPCA(umi, exprs_values = "counts")
dim(reducedDim(umi, "PCA"))
umi.qc <- runPCA(umi.qc, exprs_values = "logcounts_raw")
dim(reducedDim(umi.qc, "PCA"))
plotPCA(umi.qc, colour_by = "batch", size_by = "detected", shape_by = "individual")
dev.off()

umi.qc <- runPCA(umi.qc,ntop=14154, exprs_values = "logcounts_raw")
plotPCA(umi.qc, colour_by = "batch", size_by = "detected", shape_by = "individual")

umi.qc <- runPCA(umi.qc,ntop=50, exprs_values = "logcounts_raw")
plotPCA(umi.qc, colour_by = "batch", size_by = "detected", shape_by = "individual")
dev.off()

logcounts(umi.qc) <- assay(umi.qc, "logcounts_raw")
getExplanatoryPCs(umi.qc,variables = "sum")
plotExplanatoryPCs(umi.qc,variables = "sum") 
dev.off()
logcounts(umi.qc) <- NULL

plotExplanatoryVariables(umi.qc,exprs_values = "logcounts_raw",
                         variables = c("detected","sum","batch",
                                       "individual","altexps_ERCC_percent","subsets_Mito_percent"))
dev.off()
#------------------------------------------------------------
#counfounder
library(scRNA.seq.funcs)
library(scater)
library(scran)
library(sva)
library(batchelor)
library(kBET)
set.seed(1234567)

umi    <- readRDS("data/tung/umi.rds")
umi.qc <- umi[! rowData(umi)$discard, ! colData(umi)$discard]
qclust <- quickCluster(umi.qc, min.size = 30)
umi.qc <- computeSumFactors(umi.qc, clusters = qclust)
umi.qc <- logNormCounts(umi.qc)
assay(umi.qc, "combat") <- ComBat(logcounts(umi.qc),batch = umi.qc$replicate)
assay(umi.qc, "combat_tf") <- ComBat(logcounts(umi.qc),batch = umi.qc$detected)
mnn_out <- fastMNN(umi.qc,batch = umi.qc$replicate)
assay(umi.qc, "mnn") <- assay(mnn_out,'reconstructed')
#1.PCA
for(n in assayNames(umi.qc)) {
  tmp <- runPCA(umi.qc, exprs_values = n, ncomponents = 20)
  
  print(
    plotPCA(
      tmp,
      colour_by = "batch",
      size_by = "detected",
      shape_by = "individual"
    ) +
      ggtitle(n)
  )
}
dev.off()
#2.RLE
res <- list()
for(n in assayNames(umi.qc)) {
  res[[n]] <- suppressWarnings(calc_cell_RLE(assay(umi.qc, n)))
}
par(mar=c(6,4,1,1))
boxplot(res, las=2)
dev.off()
#3.kBET
compare_kBET_results <- function(sce){
  sce <- umi.qc
  indiv <- unique(as.character(sce$individual))
  norms <- assayNames(sce) # Get all normalizations
  results <- list()
  for (i in indiv){ 
    for (j in norms){
      tmp <- kBET(
        df = t(assay(sce[,sce$individual== i], j)), 
        batch = sce$batch[sce$individual==i], 
        heuristic = TRUE, 
        verbose = FALSE, 
        addTest = FALSE, 
        plot = FALSE)
      results[[i]][[j]] <- tmp$summary$kBET.observed[1]
    }
  }
  return(do.call(rbind.data.frame, results))
}

eff_debatching <- compare_kBET_results(umi.qc)
eff_debatching

library("reshape2")
library("RColorBrewer")
# Plot results
dod <- melt(as.matrix(eff_debatching),  value.name = "kBET")
colnames(dod)[1:2] <- c("Normalisation", "Individual")

colorset <- c('gray', brewer.pal(n = 9, "RdYlBu"))

ggplot(dod, aes(Normalisation, Individual, fill=kBET)) +  
  geom_tile() +
  scale_fill_gradient2(
    na.value = "gray",
    low = colorset[2],
    mid=colorset[6],
    high = colorset[10],
    midpoint = 0.5, limit = c(0,1)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(
    axis.text.x = element_text(
      angle = 45, 
      vjust = 1, 
      size = 12, 
      hjust = 1
    )
  ) + 
  ggtitle("Effect of batch regression methods per individual")
dev.off()
#-----------------------------------------------------
#clustering
library(pcaMethods)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)
set.seed(1234567)
#-------------------------------------------------
#DE analysis.
library(scRNA.seq.funcs)
library(edgeR)
library(monocle)
library(MAST)
library(ROCR)
set.seed(1)

DE <- read.table("data/tung/TPs.txt")
notDE <- read.table("data/tung/TNs.txt")
GroundTruth <- list(
  DE = as.character(unlist(DE)), 
  notDE = as.character(unlist(notDE))
)

molecules <- read.table("data/tung/molecules.txt", sep = "\t")
anno <- read.table("data/tung/annotation.txt", sep = "\t", header = TRUE)
keep <- anno[,1] == "NA19101" | anno[,1] == "NA19239"
data <- molecules[,keep]
group <- anno[keep,1]
batch <- anno[keep,4]
# remove genes that aren't expressed in at least 6 cells
gkeep <- rowSums(data > 0) > 5;
counts <- data[gkeep,]
# Library size normalization
lib_size = colSums(counts)
norm <- t(t(counts)/lib_size * median(lib_size)) 
# Variant of CPM for datasets with library sizes of fewer than 1 mil molecules
#-----------------------------------------------
#Seurat
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)

adj.matrix <- Read10X("data/update/soupX_pbmc10k_filt")
srat <- CreateSeuratObject(adj.matrix,project = "pbmc10k") 
srat
adj.matrix <- NULL
str(srat)

meta <- srat@meta.data
dim(meta)
summary(meta$nCount_RNA)
summary(meta$nFeature_RNA)
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")
doublets <- read.table("data/update/scrublet_calls.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat <- AddMetaData(srat,doublets)
head(srat[[]])
VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.rb")
FeatureScatter(srat, feature1 = "percent.rb", feature2 = "percent.mt")
FeatureScatter(srat, feature1 = "nFeature_RNA", feature2 = "Doublet_score")
srat[['QC']] <- ifelse(srat@meta.data$Is_doublet == 'True','Doublet','Pass')
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC == 'Pass','Low_nFeature',srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$percent.mt > 15 & srat@meta.data$QC == 'Pass','High_MT',srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'High_MT',paste('High_MT',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
table(srat[['QC']])
VlnPlot(subset(srat, subset = QC == 'Pass'), 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()

#normalization
#In order to do further analysis, we need to normalize the data to account for sequencing depth.
#Conventional way is to scale it to 10,000 (as if all cells have 10k UMIs overall),
#and log2-transform the obtained values. 
#Normalized data are stored in srat[['RNA']]@data of the 'RNA' assay.
#LogNormalize: Feature counts for each cell are divided by the total counts for that cell 
#and multiplied by the scale.factor=10000. 
#This is then natural-log transformed using log1p.
srat <- NormalizeData(srat)
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(srat), 10)
top10 
plot1 <- VariableFeaturePlot(srat)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
#ScaleData converts normalized gene expression to Z-score (values centered at 0 and with variance of 1). 
#It's stored in srat[['RNA']]@scale.data and used in following PCA. Default is to run scaling only on variable genes.
all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)
srat <- RunPCA(srat, features = VariableFeatures(object = srat))

VizDimLoadings(srat, dims = 1:9, reduction = "pca") & 
  theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))

DimHeatmap(srat, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
DimPlot(srat, reduction = "pca")
ElbowPlot(srat)
srat <- FindNeighbors(srat, dims = 1:10)
srat <- FindClusters(srat, resolution = 0.5)
srat <- RunUMAP(srat, dims = 1:10, verbose = F)
table(srat@meta.data$seurat_clusters)
DimPlot(srat,label.size = 4,repel = T,label = T)
FeaturePlot(srat, features = c("LILRA4", "TPM2", "PPBP", "GP1BB"))
FeaturePlot(srat, features = "Doublet_score") & theme(plot.title = element_text(size=10))
FeaturePlot(srat, features = "percent.mt") & theme(plot.title = element_text(size=10))
FeaturePlot(srat, features = "nFeature_RNA") & theme(plot.title = element_text(size=10))
DimPlot(srat,label.size = 4,repel = T,label = T)
dev.off()

srat <- subset(srat, subset = QC == 'Pass')
DimPlot(srat,label.size = 4,repel = T,label = T)
#cell cycle genes
cc.genes.updated.2019
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes)
table(srat[[]]$Phase)
FeaturePlot(srat,features = "percent.mt",label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
VlnPlot(srat,features = "percent.mt") & theme(plot.title = element_text(size=10))
FeaturePlot(srat,features = "percent.rb",label.size = 4,repel = T,label = T) & theme(plot.title = element_text(size=10))
VlnPlot(srat,features = "percent.rb") & theme(plot.title = element_text(size=10))
VlnPlot(srat,features = c("nCount_RNA","nFeature_RNA")) & 
  theme(plot.title = element_text(size=10))
FeaturePlot(srat,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
VlnPlot(srat,features = c("S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))
dev.off()
#Since we have performed extensive QC with doublet and empty cell removal, we can now apply SCTransform normalization, 
#that was shown to be beneficial for finding rare cell populations by improving signal/noise ratio.
#Single SCTransform command replaces NormalizeData, ScaleData, and FindVariableFeatures. 
#We will also correct for % MT genes and cell cycle scores using vars.to.regress variables; 
#our previous exploration has shown that neither cell cycle score nor MT percentage change very dramatically
#between clusters, so we will not remove biological signal, but only some unwanted variation.
srat <- SCTransform(srat, method = "glmGamPoi", ncells = 8824, 
                    vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F)
srat

FeaturePlot(srat,"PPBP") & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

FeaturePlot(srat,"LILRA4") & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()
srat <- FindNeighbors(srat, dims = 1:30, k.param = 15, verbose = F)
srat <- FindClusters(srat, verbose = F, algorithm = 4, resolution = 0.9)

FeaturePlot(srat,"MS4A1") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle("MS4A1: B cells")
FeaturePlot(srat,"LYZ") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle("LYZ: monocytes")
FeaturePlot(srat,"NKG7") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle("NKG7: natural killers")
FeaturePlot(srat,"CD8B") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle("CD8B: CD8 T cells")
FeaturePlot(srat,"IL7R") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle("IL7R: CD4 T cells")

DefaultAssay(srat) <- "RNA"
srat <- NormalizeData(srat)
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)
all.markers <- FindAllMarkers(srat, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
dim(all.markers)
table(all.markers$cluster)
top3_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
top3_markers
#cell type annotation
monaco.ref <- celldex::MonacoImmuneData()
# hpca.ref <- celldex::HumanPrimaryCellAtlasData()
# dice.ref <- celldex::DatabaseImmuneCellExpressionData()
sce <- as.SingleCellExperiment(DietSeurat(srat))
sce
monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
# hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
# hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
# dice.main <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
# dice.fine <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)
table(monaco.fine$pruned.labels)
srat@meta.data$monaco.main <- monaco.main$pruned.labels
srat@meta.data$monaco.fine <- monaco.fine$pruned.labels

# srat@meta.data$hpca.main   <- hpca.main$pruned.labels
# srat@meta.data$dice.main   <- dice.main$pruned.labels
# srat@meta.data$hpca.fine   <- hpca.fine$pruned.labels
# srat@meta.data$dice.fine   <- dice.fine$pruned.labels
srat <- SetIdent(srat, value = "monaco.fine")
DimPlot(srat, label = T , repel = T, label.size = 3) + NoLegend()
FeaturePlot(srat,"CD38") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
table(srat[[]]$seurat_clusters)

srat <- RunPCA(srat, verbose = F)
srat <- RunUMAP(srat, dims = 1:30, verbose = F)
srat <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat <- FindClusters(srat, verbose = F)
table(srat[[]]$seurat_clusters)
DimPlot(srat, label = T)