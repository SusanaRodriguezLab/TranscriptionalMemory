## LIBRARIES TO USE
library(NOISeq)
library(plotly)
library(maSigPro)

## Raw Counts
COUNTS <- read.csv2("CountMatrix.csv", row.names = 1)
gen_len <- read.table("sacCer3_genelenghts.txt")
gen_biotypes <- read.table("sacCer3_biotypes.txt")
coldata <- read.csv2("ColData.csv", row.names = 1)

## NOISeq QC ------------
mydata <- readData(data = COUNTS, factors = coldata, length = gen_len, biotype = gen_biotypes)

# Biotype detection
mybiotype <- dat(mydata, k=0, type = "biodetection", factor = NULL)
explo.plot(mybiotype)

# RNA composition
mycd<- dat(mydata, type = "cd", norm = F, refColumn = 1)
explo.plot(mycd, samples = 1:12)

# Saturation plot
mysaturation <- dat(mydata, k = 0, ndepth = 7, type = "saturation")
explo.plot(mysaturation, toplot = 1, samples = 1:36)

# Length bias
mylenbias <- dat(mydata, factor = "Rep", type = "lengthbias")
explo.plot(mylenbias, samples = NULL, toplot="global")

# Batch exploration
myPCA <- dat(mydata, type = "PCA")
explo.plot(myPCA, factor = "Rep")

# CPM 
mycounts <- dat(mydata, factor = NULL, type = "countsbio")
explo.plot(mycounts, toplot=1, samples = NULL, plottype = "barplot")

## TMM normalisation and filtering
myTMM <- tmm(assayData(mydata)$exprs, long = mylength_1)
myfilt <- filtered.data(myTMM, norm = T, method = 1, cpm = 1, factor = "Full")

# Batch correction
mydata2 <- readData(data = myfilt, factors = coldata)
mydata2 <- ARSyNseq(mydata2, factor = "Rep", batch = T, norm = "n", logtransf = F,  Variability = 0.99)

## Export TMM normalised counts
write.table(assayData(mydata2)$exprs, row.names = T, col.names = T,
            sep = "\t", file = "Normalised_counts.txt")

## Pairwise comparisons with NOISeq
res15_WT_no <- noiseqbio(mydata2, factor = "Full", conditions = c("WT_15_No", "WT_0_No"), norm = "n",random.seed = 123 )
degs15_WT_no <- degenes(res15_WT_no, q = 0.95, M= NULL)
res15_WT_no_df <- res15_WT_no@results[[1]]

res20_WT_no <- noiseqbio(mydata2, factor = "Full", conditions = c("WT_20_No", "WT_0_No"), norm = "n",random.seed = 123 )
degs20_WT_no <- degenes(res20_WT_no, q = 0.95, M= NULL)
res20_WT_no_df <- res20_WT_no@results[[1]]

res15_WT_yes <- noiseqbio(mydata2, factor = "Full", conditions = c("WT_15_Yes", "WT_0_Yes"), norm = "n",random.seed = 123 )
degs15_WT_yes <- degenes(res15_WT_yes, q = 0.95, M= NULL)
res15_WT_yes_df <- res15_WT_yes@results[[1]]

res20_WT_yes <- noiseqbio(mydata2, factor = "Full", conditions = c("WT_20_Yes", "WT_0_Yes"), norm = "n",random.seed = 123 )
degs20_WT_yes <- degenes(res20_WT_yes, q = 0.95, M= NULL)
res20_WT_yes_df <- res20_WT_yes@results[[1]]

res15_mip6_yes <- noiseqbio(mydata2, factor = "Full", conditions = c("mip6_15_Yes", "mip6_0_Yes"), norm = "n",random.seed = 123 )
degs15_mip6_yes <- degenes(res15_mip6_yes, q = 0.95, M= NULL)
res15_mip6_yes_df <- res15_mip6_yes@results[[1]]

res20_mip6_yes <- noiseqbio(mydata2, factor = "Full", conditions = c("mip6_20_Yes", "mip6_0_Yes"), norm = "n",random.seed = 123 )
degs20_mip6_yes <- degenes(res20_mip6_yes, q = 0.95, M= NULL)
res20_mip6_yes_df <- res20_mip6_yes@results[[1]]

res15_mip6_no <- noiseqbio(mydata2, factor = "Full", conditions = c("mip6_15_No", "mip6_0_No"), norm = "n",random.seed = 123 )
degs15_mip6_no <- degenes(res15_mip6_no, q = 0.95, M= NULL)
res15_mip6_no_df <- res15_mip6_no@results[[1]]

res20_mip6_no <- noiseqbio(mydata2, factor = "Full", conditions = c("mip6_20_No", "mip6_0_No"), norm = "n",random.seed = 123 )
degs20_mip6_no <- degenes(res20_mip6_no, q = 0.95, M= NULL)
res20_mip6_no_df <- res20_mip6_no@results[[1]]

## Function to generate Sankey plot for two given DESeq2 results objects
generate_sankey_plot <- function(res1, res2) {
  # Ensure gene identifiers are included
  res1$gene <- rownames(res1)
  res2$gene <- rownames(res2)
  
  # Categorize fold changes
  categorize_fc <- function(df) {
    df %>%
      mutate(category = case_when(
        1 - prob >= 0.05 ~ 'unchanged',
        1 - prob < 0.05 & log2FC > 0 & log2FC <= log(1.2)/log(2) ~ 'unchanged',
        1 - prob < 0.05 & log2FC < -0 & log2FC >= -log(1.2)/log(2) ~ 'unchanged',
        1 - prob < 0.05 & log2FC > log(1.2)/log(2) & log2FC <= 1 ~ 'up',
        1 - prob < 0.05 & log2FC < -log(1.2)/log(2) & log2FC >= -1 ~ 'down',
        1 - prob < 0.05 & log2FC > 1 ~ 'strong up',
        1 - prob < 0.05 & log2FC < -1 ~ 'strong down'
      ))
  }
  
  res1 <- categorize_fc(res1)
  res2 <- categorize_fc(res2)
  
  # Extract relevant columns
  res1 <- res1[, c("category", "gene")]
  colnames(res1) <- c("category_1", "gene")
  
  res2 <- res2[, c("category", "gene")]
  colnames(res2) <- c("category_2", "gene")
  
  # Merge results by gene
  merged_res <- merge(res1, res2, by = "gene")
  
  # Add a default category for the initial time point
  merged_res <- merged_res %>%
    mutate(category_0 = 'unchanged')
  
  # Create transitions for initial to time point 1 and time point 1 to time point 2
  transitions_0_1 <- merged_res %>%
    group_by(category_0, category_1) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    mutate(time = "0_1")
  
  transitions_1_2 <- merged_res %>%
    group_by(category_1, category_2) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    mutate(time = "1_2")
  
  # Combine the transitions into one dataframe
  transitions <- bind_rows(transitions_0_1, transitions_1_2)
  
  # Define unique nodes and their positions
  nodes <- c("unchanged_0", "unchanged_1", "up_1", "down_1", "strong up_1", "strong down_1", "unchanged_2", "up_2", "down_2", "strong up_2", "strong down_2")
  node_map <- setNames(seq_along(nodes), nodes)
  
  # Add the node names to transitions dataframe
  transitions <- transitions %>%
    mutate(source = case_when(
      time == "0_1" & category_0 == "unchanged" ~ node_map["unchanged_0"],
      time == "1_2" & category_1 == "unchanged" ~ node_map["unchanged_1"],
      time == "1_2" & category_1 == "up" ~ node_map["up_1"],
      time == "1_2" & category_1 == "down" ~ node_map["down_1"],
      time == "1_2" & category_1 == "strong up" ~ node_map["strong up_1"],
      time == "1_2" & category_1 == "strong down" ~ node_map["strong down_1"]
    ),
    target = case_when(
      time == "0_1" & category_1 == "unchanged" ~ node_map["unchanged_1"],
      time == "0_1" & category_1 == "up" ~ node_map["up_1"],
      time == "0_1" & category_1 == "down" ~ node_map["down_1"],
      time == "0_1" & category_1 == "strong up" ~ node_map["strong up_1"],
      time == "0_1" & category_1 == "strong down" ~ node_map["strong down_1"],
      time == "1_2" & category_2 == "unchanged" ~ node_map["unchanged_2"],
      time == "1_2" & category_2 == "up" ~ node_map["up_2"],
      time == "1_2" & category_2 == "down" ~ node_map["down_2"],
      time == "1_2" & category_2 == "strong up" ~ node_map["strong up_2"],
      time == "1_2" & category_2 == "strong down" ~ node_map["strong down_2"]
    ))
  
  # Calculate the number of genes in each category for each time point
  node_counts <- merged_res %>%
    select(starts_with("category")) %>%
    gather(key = "time", value = "category") %>%
    group_by(time, category) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    mutate(node_label = paste0(category, "_", sub("category_", "", time)))
  
  # Create a lookup table for the node labels
  node_labels <- sapply(nodes, function(node) {
    count <- node_counts %>% filter(node_label == node) %>% pull(n)
    if (length(count) == 0) {
      return(node)
    } else {
      # Use HTML tags to increase the size of the labels
      return(paste0("<b><span style='font-size:16px;'>", as.character(count), "</span></b>"))
    }
  })
  
  # Define colors for each node
  color_map <- c(
    'unchanged_0' = 'grey',
    'unchanged_1' = 'grey',
    'up_1' = 'lightgreen',
    'down_1' = 'orange',
    'strong up_1' = 'green',
    'strong down_1' = 'red',
    'unchanged_2' = 'grey',
    'up_2' = 'lightgreen',
    'down_2' = 'orange',
    'strong up_2' = 'green',
    'strong down_2' = 'red'
  )
  
  # Create Sankey plot
  fig <- plot_ly(
    type = "sankey",
    arrangement = "snap",
    node = list(
      label = node_labels,
      pad = 15,
      thickness = 20,
      color = unname(color_map[nodes]),
      # Add this line for enabling HTML in labels
      hoverlabel = list(font = list(size = 20)),
      line = list(color = "black", width = 0.5)
    ),
    link = list(
      source = transitions$source - 1,  # Plotly uses 0-based indexing
      target = transitions$target - 1,
      value = transitions$n,
      color = "#ededed"
    )
  )
  
  fig <- fig %>%
    layout(
      title = "Sankey Diagram of Gene Expression Changes",
      font = list(size = 12),
      xaxis = list(
        title = "",
        showticklabels = FALSE
      ),
      yaxis = list(
        title = "",
        showticklabels = FALSE
      )
    )
  
  # Display the plot
  fig
}

## Generate Sankey plots for pairwise comparisons
generate_sankey_plot(res15_WT_no_df, res20_WT_no_df)
generate_sankey_plot(res15_WT_yes_df, res20_WT_yes_df)
generate_sankey_plot(res15_mip6_no_df, res20_mip6_no_df)
generate_sankey_plot(res15_mip6_yes_df, res20_mip6_yes_df)

## MaSigPro analysis

# Memory vs No Memory in WT
wt.cols <- grep("wt", colnames(normalised_counts))
counts_soloWT <- normalised_counts[,wt.cols]
Time <- rep(c(0,15,20),6) ; Time
Replicates <- rep(c(1,2,3,4,5,6), 3) ; Replicates
me<- rep(c(rep(1,3), rep(0,3)),3) ; me
sin <- rep(c(rep(0,3), rep(1,3)),3) ; sin
design.wt <- data.frame(Time = Time, Replicates = Replicates, me = me, sin = sin)
rownames(design.wt) <- colnames(counts_soloWT)
design2 <- make.design.matrix(design.wt, degree = 2)
fit_soloWT <- p.vector(counts_soloWT, design = design2, Q = 0.05, MT.adjust = "BH", min.obs = 6)
tstep_soloWT <- T.fit(fit_soloWT, step.method = "backward", alfa = 0.05)
sigs.groups_soloWT <- get.siggenes(tstep_soloWT, rsq = 0.6, vars = "groups")
sigs.all_soloWT <- get.siggenes(tstep_soloWT, rsq = 0.6, vars = "all")

# Memory vs No Memory in mip6
counts_solomip6 <- normalised_counts[,-wt.cols]
Time <- rep(c(0,15,20),6) ; Time
Replicates <- rep(c(1,2,3,4,5,6), 3) ; Replicates
me<- rep(c(rep(1,3), rep(0,3)),3) ; me
sin <- rep(c(rep(0,3), rep(1,3)),3) ; sin
design.mip <- data.frame(Time = Time, Replicates = Replicates, me = me, sin = sin)
rownames(design.mip) <- colnames(counts_solomip6)
design3 <- make.design.matrix(design.mip, degree = 2)
fit_solomip6 <- p.vector(counts_solomip6, design = design3, Q = 0.05, MT.adjust = "BH", min.obs = 6)
tstep_solomip6 <- T.fit(fit_solomip6, step.method = "backward", alfa = 0.05)
sigs.groups_solomip6 <- get.siggenes(tstep_solomip6, rsq = 0.6, vars = "groups")
sigs.all_solomip6 <- get.siggenes(tstep_solomip6, rsq = 0.6, vars = "all")

# Model with all groups for cluster visualisation

design <- read_excel("edesign.xlsx")

edesign <- design %>%
  select(-...1)
rownames(edesign) <- design$...1

edesign <- as.matrix(edesign)
design <- make.design.matrix(edesign, degree = 2  )
design$groups.vector
fit <- p.vector(normalised_counts, design, Q = 0.05, MT.adjust = "BH", min.obs = 6)
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
sigs.groups <- get.siggenes(tstep, rsq = 0.6, vars = "groups")
sigs.all <- get.siggenes(tstep, rsq = 0.6, vars = "all")
memvsnomem <- see.genes(sigs.groups$sig.genes$WT_con_memoriavsWT_sin_memoria, show.fit = T, dis =design$dis,
                        cluster.method="hclust", legend = F,show.lines=TRUE, k = 4)
lista_clusters_memvsnomem <- data.frame(memvsnomem$cut)
lista_clusters_memvsnomem_1 <- lista_clusters_memvsnomem
lista_clusters_memvsnomem_1$names <- rownames(lista_clusters_memvsnomem_1)
for (i in 1:4){
  cluster <- rownames(subset(lista_clusters_memvsnomem, lista_clusters_memvsnomem$memvsnomem.cut == as.character(i)))
  write.table(data.frame(cluster), paste0("cluster",as.character(i),".txt"), row.names = F, col.names = F, quote = F)
}