# Load required libraries
library(NOISeq)   # Quality control and differential expression analysis for RNA-seq
library(plotly)   # Generating interactive plots, including Sankey diagrams
library(maSigPro) # Time-series differential expression analysis
library(dplyr)    # Data manipulation and filtering
library(writexl)  # Exporting data to Excel files
library(openxlsx) # Handling Excel workbooks
library(tidyr)    # Data transformation
library(BiocManager)  # Bioconductor package manager
library(AnnotationDbi)  # Functions for annotation databases
library(org.Sc.sgd.db)  # Yeast genome annotation database
library(clusterProfiler) # Gene Ontology (GO) enrichment analysis

# Load gene expression count data
COUNTS <- read.csv2("CountMatrix.csv", row.names = 1)  # Expression count matrix
gen_len <- read.table("sacCer3_genelenghts.txt")  # Gene length information
gen_biotypes <- read.table("sacCer3_biotypes.txt")  # Gene biotype annotations
coldata <- read.csv2("ColData.csv", row.names = 1)  # Metadata for experimental conditions

# Create NOISeq data object for quality control
mydata <- readData(data = COUNTS, factors = coldata, length = gen_len, biotype = gen_biotypes)

# Biotype detection analysis
mybiotype <- dat(mydata, k = 0, type = "biodetection", factor = NULL)
explo.plot(mybiotype)  # Plot biotype distribution

# RNA composition analysis without normalization
mycd <- dat(mydata, type = "cd", norm = FALSE, refColumn = 1)
explo.plot(mycd, samples = 1:12)  # Plot composition for the first 12 samples

# Saturation plot to evaluate sequencing depth
mysaturation <- dat(mydata, k = 0, ndepth = 7, type = "saturation")
explo.plot(mysaturation, toplot = 1, samples = 1:36)

# Length bias analysis to check for systematic bias
mylenbias <- dat(mydata, factor = "Rep", type = "lengthbias")
explo.plot(mylenbias, samples = NULL, toplot = "global")

# Principal Component Analysis (PCA) for batch effect exploration
myPCA <- dat(mydata, type = "PCA")
explo.plot(myPCA, factor = "Rep")

# Count per million (CPM) distribution
mycounts <- dat(mydata, factor = NULL, type = "countsbio")
explo.plot(mycounts, toplot = 1, samples = NULL, plottype = "barplot")

# TMM normalization and filtering of low-count genes
myTMM <- tmm(assayData(mydata)$exprs, long = gen_len)  # Perform TMM normalization
myfilt <- filtered.data(myTMM, norm = TRUE, method = 1, cpm = 1, factor = "Full")  # Filter low-count genes

# Apply batch correction using ARSyNseq
mydata2 <- readData(data = myfilt, factors = coldata)
mydata2 <- ARSyNseq(mydata2, factor = "Rep", batch = TRUE, norm = "n", logtransf = FALSE, Variability = 0.99)

# Export TMM normalized counts
write.table(assayData(mydata2)$exprs, row.names = TRUE, col.names = TRUE,
            sep = "\t", file = "Normalised_counts.txt")

# Function for GO enrichment analysis
get_GOs <- function(gene_list) {
  if (!requireNamespace("org.Sc.sgd.db", quietly = TRUE)) {
    stop("Package org.Sc.sgd.db is not installed. Please install it using BiocManager.")
  }
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package clusterProfiler is not installed. Please install it using BiocManager.")
  }
  
  library(clusterProfiler)
  library(org.Sc.sgd.db)
  
  valid_orfs <- gene_list[gene_list %in% keys(org.Sc.sgd.db, keytype = "ORF")]
  
  if (length(valid_orfs) == 0) {
    warning("No valid ORFs found for enrichment analysis.")
    return(NULL)
  }
  
  ego <- enrichGO(
    gene = valid_orfs,
    OrgDb = org.Sc.sgd.db,
    keyType = "ORF",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = FALSE
  )
  
  return(ego)
}

# Create a design matrix for differential expression analysis
design.mat <- model.matrix(~ (Time + I(Time^2)) * Memory * Strain, data = edesign)

design.df <- as.data.frame(design.mat)
names(design.df) <- make.names(names(design.df))

# Differential expression analysis using regression model
fit <- p.vector(
  data = normalised_counts,
  design = design.df,
  Q = 0.05,
  MT.adjust = "BH"
)

step.fit <- T.fit(
  fit,
  alfa = 0.05,
  step.method = "backward"
)

# Export GO enrichment results to an Excel file
wb <- createWorkbook()
categories <- c("induction_no change", "repression_no change", "repression_enhanced", "induction_enhanced",
                "repression_dampened", "induction_dampened", "induction_reverted", "repression_reverted")

for (category in categories) {
  if (!is.null(GO_results_filtered[[category]])) {
    df <- GO_results_filtered[[category]]@result
    addWorksheet(wb, category)
    writeData(wb, category, df)
  }
}

saveWorkbook(wb, "GO_results_filtered.xlsx", overwrite = TRUE)
