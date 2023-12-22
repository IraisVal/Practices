#Human DLPFC workflow
#This workflow analyzes one sample of human brain from the dorsolateral prefrontal cortex (DLPFC) region, 
#measured using the 10x Genomics Visium platform
#12 samples in total, from 3 individuals, with 2 pairs of spatially adjacent replicates (serial sections) per individual 
#(4 samples per individual).
#Workflow https://lmweber.org/BestPracticesST/chapters/workflow-Visium-humanDLPFC.html

library(SpatialExperiment)
library(STexampleData)

# load object
spe <- Visium_humanDLPFC()
spe

library(ggspavis)
# plot spatial coordinates (spots)
plotSpots(spe)

spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)

library(scater)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)

rowData(spe)$gene_name[is_mito]

# calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
head(colData(spe), 3)

# histograms of QC metrics
par(mfrow = c(1, 4))
hist(colData(spe)$sum, xlab = "sum", main = "UMIs per spot")
hist(colData(spe)$detected, xlab = "detected", main = "Genes per spot")
hist(colData(spe)$subsets_mito_percent, xlab = "percent mitochondrial", main = "Percent mito UMIs")
hist(colData(spe)$cell_count, xlab = "number of cells", main = "No. cells per spot")

par(mfrow = c(1, 1))

# select QC thresholds
qc_lib_size <- colData(spe)$sum < 600
qc_detected <- colData(spe)$detected < 400
qc_mito <- colData(spe)$subsets_mito_percent > 28
qc_cell_count <- colData(spe)$cell_count > 10

# number of discarded spots for each metric
apply(cbind(qc_lib_size, qc_detected, qc_mito, qc_cell_count), 2, sum)

# combined set of discarded spots
discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count
table(discard)

# store in object
colData(spe)$discard <- discard

# check spatial pattern of discarded spots
plotQC(spe, type = "spots", discard = "discard")

# filter low-quality spots
spe <- spe[, !colData(spe)$discard]
dim(spe)

library(scran)
# calculate library size factors
spe <- computeLibraryFactors(spe)

summary(sizeFactors(spe))

hist(sizeFactors(spe), breaks = 20)

hist(sizeFactors(spe), breaks = 20)

# calculate logcounts and store in object
spe <- logNormCounts(spe)

assayNames(spe)

# remove mitochondrial genes
spe <- spe[!is_mito, ]
dim(spe)

# fit mean-variance relationship
dec <- modelGeneVar(spe)

# visualize mean-variance relationship
fit <- metadata(dec)
plot(fit$mean, fit$var, 
     xlab = "mean of log-expression", ylab = "variance of log-expression")
curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

# select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)

library(nnSVG)

# subsample spots
n <- 100
set.seed(123)
ix <- sample(seq_len(n), n)

spe_nnSVG <- spe[, ix]

# filter low-expressed and mitochondrial genes
# using very stringent filtering parameters for faster runtime in this example
# note: for a full analysis, use alternative filtering parameters (e.g. defaults)
spe_nnSVG <- filter_genes(
  spe_nnSVG, filter_genes_ncounts = 10, filter_genes_pcspots = 3
)

# re-calculate logcounts after filtering
# using library size factors
spe_nnSVG <- logNormCounts(spe_nnSVG)

# run nnSVG
# using a single core for compatibility on build system
# note: for a full analysis, use multiple cores
set.seed(123)
spe_nnSVG <- nnSVG(spe_nnSVG, n_threads = 1)

# investigate results

# show results
#head(rowData(spe_nnSVG), 3)

# number of significant SVGs
table(rowData(spe_nnSVG)$padj <= 0.05)

# show results for top n SVGs
rowData(spe_nnSVG)[order(rowData(spe_nnSVG)$rank)[1:6], ]

# identify top-ranked SVG
rowData(spe_nnSVG)$gene_name[which(rowData(spe_nnSVG)$rank == 1)]

# compute PCA
set.seed(123)
spe <- runPCA(spe, subset_row = top_hvgs)

reducedDimNames(spe)

dim(reducedDim(spe, "PCA"))

# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")

reducedDimNames(spe)
#Next step Clustering
