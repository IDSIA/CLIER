if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require("recount3", quietly = TRUE)) {
  BiocManager::install("recount3")
  library(recount3)
}

if (!require("recount", quietly = TRUE)) {
  BiocManager::install("recount")
  library(recount)
}

if (!require("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
  library(biomaRt)
}

if(!require(devtools, quietly = TRUE)){
  install.packages("devtools")
  library(devtools)
}

if(!require(PLIER, quietly = TRUE)){
  install_github("wgmao/PLIER")
  library(PLIER)
}

if(!require(stringr, quietly = TRUE)){
  install.packages("stringr")
  library(stringr)
}

if(!require(readxl, quietly = TRUE)){
  install.packages("readxl")
  library(readxl)
}

if(!require(writexl, quietly = TRUE)){
  install.packages("writexl")
  library(writexl)
}

if(!require(data.table, quietly = TRUE)){
  install.packages("data.table")
  library(data.table)
}

if(!require(ggplot2, quietly = TRUE)){
  install.packages("ggplot2")
  library(ggplot2)
}

if(!require(networkD3, quietly = TRUE)){
  install.packages("networkD3")
  library(networkD3)
}

if(!require(htmlwidgets, quietly = TRUE)){
  install.packages("htmlwidgets")
  library(htmlwidgets)
}


###### DOWNLOADING TCGA

#' Download the TCGA dataset in TPM
#'
#' @param save The file path to save the dataset.
#' @return A data.table containing the TCGA dataset.
downloadTcga <- function(save = NULL) {
  
  # Select TCGA studies among all studies available in recount3
  human_projects <- available_projects()
  projects <- human_projects[human_projects$project_home %like% "tcga", ]$project
  #print(projects)
  
  # Iteratively build the full gene expression matrix
  full_mat <- NULL
  for (pr in projects) {
    # Project being retrieved
    proj_info <- subset(human_projects, project_type == "data_sources" & project == pr)
    print(paste("Retrieving", pr))
    
    # Create an RSE for each study and get counts and TPM
    rse_gene <- create_rse(proj_info)
    assays(rse_gene)$counts <- transform_counts(rse_gene)
    assays(rse_gene)$tpm <- recount::getTPM(rse_gene, length_var = "bp_length")
    full_mat <- cbind(full_mat, assays(rse_gene)$tpm)
  }
  
  # Create a dataset that contains all the TPM
  dt_tcga <- as.data.table(full_mat, keep.rownames = "ensembl_gene_id")[order(ensembl_gene_id)]
  
  # Save the dataset
  if (!is.null(save)) {
    saveRDS(dt_tcga, save)
  }
  
  return(dt_tcga)
}

###### PREPROCESSING

### Main Functions

#' Preprocess a TPM dataset for training KCLIER
#'
#' @param dt The TPM dataset.
#' @param atlas The atlas to be used.
#' @param convert Logical; whether to convert gene IDs.
#' @param normalize Logical; whether to normalize the data.
#' @return A processed dataset.
trainPreprocess <- function(dt, atlas, convert = TRUE, normalize = FALSE) {
  dt1 <- initialFilters(dt)
  
  dt2 <- if (convert) {
    hgncConversion(dt1)
  } else {
    dt1
  }
  
  dt3 <- if (normalize) {
    renorm(dt2)
  } else {
    dt2
  }
  
  mat4 <- trainFurtherFilter(dt2)
  mat5 <- trainSignaturesFilter(mat4, atlas)
  out6 <- trainLogStd(mat5)
  
  return(out6)
}

#' Preprocess a TPM dataset for applying KCLIER
#'
#' @param dt The TPM dataset.
#' @param m_train The mean values from training data.
#' @param sd_train The standard deviation values from training data.
#' @param convert Logical; whether to convert gene IDs.
#' @param normalize Logical; whether to normalize the data.
#' @return A processed dataset.
testPreprocess <- function(dt, dt_stats, convert = TRUE, normalize = FALSE) {
  dt1 <- initialFilters(dt)
  
  dt2 <- if (convert) {
    hgncConversion(dt1)
  } else {
    dt1
  }
  
  dt3 <- if (normalize) {
    renorm(dt2)
  } else {
    dt2
  }
  
  mat4 <- testFilter(dt2, dt_stats)
  mat5 <- testLogStd(mat4, dt_stats)
  
  return(mat5)
}

### Train and Test Common Functions

#' Apply simple filters on NAs and genes that are identically zero
#'
#' @param dt The dataset to be filtered.
#' @return The filtered dataset.
initialFilters <- function(dt) {
  print(paste(nrow(dt), "genes and", ncol(dt) - 1, "samples initially."))
  
  # Remove NAs
  dt <- dt[complete.cases(dt), ]
  print(paste(nrow(dt), "genes and", ncol(dt) - 1, "samples after NA removal."))
  
  # Remove genes which are all zero
  #filter_col <- rowSums(abs(as.matrix(dt[, -1]))) > 0
  #dt <- dt[filter_col, ]
  #print(paste(nrow(dt), "genes and", ncol(dt) - 1, "samples after all-zero genes removal."))
  
  return(dt)
}

#' Fetch a dataset for conversion from ensembl_gene_id to hgnc_symbol
#'
#' @return A data.table with ensembl_gene_id and hgnc_symbol.
getConvDataset <- function() {
  # Connect to the Ensembl BioMart database
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Define the attributes to retrieve
  attributes <- c("ensembl_gene_id", "hgnc_symbol")
  
  # Retrieve the data
  dt_conv <- data.table(getBM(attributes = attributes, mart = ensembl))
  
  # Basic cleaning and sorting
  dt_conv <- dt_conv[hgnc_symbol != ""]
  dt_conv <- dt_conv[order(ensembl_gene_id, hgnc_symbol)]
  
  # Check for non-univocity and make conversion one-to-one
  dt_conv <- dt_conv[, .SD[1], by = ensembl_gene_id]
  dt_conv <- dt_conv[, .SD[1], by = hgnc_symbol]
  
  return(dt_conv)
}

#' Perform conversion from ensembl id to hgnc symbol
#'
#' @param dt The dataset with ensembl_gene_id.
#' @return The dataset with hgnc_symbol.
hgncConversion <- function(dt) {
  # Remove .number from geneID (if present)
  dt$ensembl_gene_id <- gsub("\\..*", "", dt$ensembl_gene_id)
  
  # Keep only one gene among those identical up to the point
  dt <- dt[, .SD[1], by = ensembl_gene_id]
  
  dt_conv <- getConvDataset()
  dt_conv <- dt_conv[ensembl_gene_id %in% dt$ensembl_gene_id]
  
  # Merge
  dt_hgnc <- merge(dt_conv, dt, by = "ensembl_gene_id")[, !c("ensembl_gene_id"), with = FALSE]
  
  print(paste(nrow(dt_hgnc), "genes and", ncol(dt_hgnc) - 1, "samples after conversion from ensembl_gene_id to hgnc_symbol."))
  
  return(dt_hgnc)
}

# Normalizes dataset at samples level
#'
#' @param dt The dataset to normalize.
#' @return A normalized dataset.
renorm <- function(dt) {
  samples <- setdiff(names(dt), "hgnc_symbol")
  dt_norm <- copy(dt)
  dt_norm[, (samples) := lapply(.SD, function(x) x / sum(x) * 10^6), .SDcols = samples]
  
  return(dt_norm)
}

### Train Specific Functions

# Applies filters to training dataset, e.g., removing low variance genes
#'
#' @param dt The training dataset.
#' @param thr The threshold for filtering.
#' @return A filtered dataset.
trainFurtherFilter <- function(dt, thr = 0.95) {
  genes_charact <- intersect(names(dt), c("hgnc_symbol", "ensembl_gene_id"))
  mat <- as.matrix(dt, rownames = genes_charact)
  
  # Remove rows for which more than 95% of samples have zero counts
  n_zeros <- rowSums(mat == 0)
  mat <- mat[(n_zeros < ncol(mat) * thr), ]
  print(paste(nrow(mat), "genes and", ncol(mat), "samples after removal of genes with large number of zeros."))
  
  return(mat)
}

# Restricts to genes contained both in the training dataset and in the signatures matrix
#'
#' @param mat_genes The matrix of genes in the training dataset.
#' @param mat_atlas The matrix containing all signatures.
#' @return A filtered matrix of genes.
trainSignaturesFilter <- function(mat_genes, mat_atlas) {
  genes_paths <- rownames(mat_atlas)
  genes_data <- rownames(mat_genes)
  
  genes_common <- intersect(genes_paths, genes_data)
  
  mat_genes <- mat_genes[genes_common, ]
  
  print(paste(nrow(mat_genes), "genes and", ncol(mat_genes), "samples after filtering genes in no signature."))
  
  return(mat_genes)
}

# Performs standardization on training dataset
#'
#' @param mat The training dataset matrix.
#' @return A list containing the standardized matrix, means, and standard deviations.
trainLogStd <- function(mat) {
  mat_log <- log2(mat + 0.5)
  
  m <- apply(mat_log, 1, mean)
  sd <- apply(mat_log, 1, sd)
  mat_log_zscore <- (mat_log - m) / sd
  
  dt_stats <- data.table(Gene = rownames(mat_log), Mean = m, StdDev = sd)
  
  return(list(Y = mat_log_zscore, dt_stats = dt_stats))
}

### Test Specific Functions

# Restricts dataset to which to apply KCLIER to selection of genes
#'
#' @param dt The dataset to be filtered.
#' @param v A vector of gene names.
#' @return A filtered matrix.
testFilter <- function(dt, dt_stats) {
  
  nonzero <- rowSums(as.matrix(dt[,-1]))!=0
  dt <- dt[nonzero]
  dt_filter <- dt[hgnc_symbol %in% dt_stats$Gene]
  
  mat <- as.matrix(dt_filter, rownames = "hgnc_symbol")
  
  print(paste(nrow(mat), "genes and", ncol(mat), "samples after taking only genes common to train."))
  
  return(mat)
}

# Performs standardization on dataset to apply KCLIER
#'
#' @param mat The matrix to be standardized.
#' @param m The mean values from the training data.
#' @param sd The standard deviation values from the training data.
#' @return The standardized matrix.
testLogStd <- function(mat, dt_stats) {
  common <- sort(intersect(dt_stats$Gene, rownames(mat)))
  mat <- mat[common, ]
  dt_stats <- dt_stats[Gene %in% common][order(Gene)]
  
  mat_log <- log2(mat + 0.5)
  mat_log_zscore <- (mat_log - dt_stats$Mean) / dt_stats$StdDev
  
  return(mat_log_zscore)
}

##### POST-PROCESSING

# Renames LVs in objects in PLIER decomposition (and pathway to Signature in summary)
#'
#' @param plier The PLIER object.
#' @return The PLIER object with renamed LVs.
renameLvs <- function(plier) {
  n_lvs <- nrow(plier$B)
  n_digits <- floor(log10(n_lvs)) + 1
  lvs_names <- paste0("LV", str_pad(1:n_lvs, n_digits, pad = "0"))
  
  rownames(plier$B) <- lvs_names
  colnames(plier$Z) <- lvs_names
  colnames(plier$U) <- lvs_names
  colnames(plier$Uauc) <- lvs_names
  if (length(plier$summary$pathway) > 0){ 
    setnames(plier$summary, "pathway", "Signature")
  }
  plier$summary$`LV name` <- paste0("LV", str_pad(plier$summary$`LV index`, n_digits, pad = "0"))
  
  return(plier)
}

# Extracts U matrix elements and complements them with info such as AUC, FDR, p-value
#'
#' @param plier The PLIER object.
#' @param auc_thr The AUC threshold.
#' @param fdr_thr The FDR threshold.
#' @param save The file path to save the information.
#' @return A data.table with U matrix information.
reshapeUInfo <- function(plier, auc_thr = 0.7, fdr_thr = 0.05, save = NULL) {
  dt_plier <- data.table(plier$summary)
  dt_plier <- dt_plier[AUC > auc_thr & FDR < fdr_thr]
  
  sel_lvs <- unique(dt_plier$`LV name`)
  U <- plier$U[, sel_lvs]
  dt_U <- data.table(U, keep.rownames = "Signature")
  dt_melted <- melt(dt_U, measure.vars = sel_lvs, variable.name = "LV name", value.name = "U")
  dt_dict <- merge(dt_plier, dt_melted[U > 0], by = c("Signature", "LV name"))
  dt_dict <- dt_dict[, .(`LV name`, Signature, U, AUC, FDR, `p-value`)][order(`LV name`, -U)]
  
  if (!is.null(save)) {
    write_xlsx(dt_dict, save)
  }
  
  return(dt_dict)
}

##### APPLICATION AND CHECKS

# Applies PLIER to new data
#'
#' @param plier The PLIER object.
#' @param newdata The new dataset.
#' @param save The file path to save the results.
#' @return The B matrix for the new data.
getB <- function(plier, newdata, save = NULL) {
  common <- intersect(rownames(plier$Z), rownames(newdata))
  
  print(dim(plier$Z))  
  Z <- plier$Z[common, ]
  
  L2 <- plier$L2
  newB <- solve(t(Z) %*% Z + L2 * diag(ncol(Z))) %*% t(Z) %*% newdata[common, ]
  rownames(newB) <- rownames(plier$B)
  
  if (!is.null(save)) {
    saveRDS(newB, save)    
  }
  
  return(newB)
}

# Computes correlations of reconstructed gene expressions with original values
#'
#' @param plier The PLIER object.
#' @param Y The original gene expression matrix.
#' @param B The B matrix (optional). If not provided, the one carried by the PLIER object.
#' @return A data.table with gene-wise correlations.
extractCor <- function(plier, Y, B = plier$B) {
  common <- intersect(rownames(plier$Z), rownames(Y))
  Z <- plier$Z[common, ]
  Y <- Y[common, sort(colnames(Y))]
  B <- B[, sort(colnames(B))]
  
  Y_rec <- Z %*% B
  Y_rec <- Y_rec #[common, ]
  
  cor_val <- sapply(1:nrow(Y), function(i) cor(Y[i, ], Y_rec[i, ], method = "spearman"))
  dt_cor <- data.table(Gene = rownames(Y), Correlation = cor_val)
  
  return(dt_cor)
}

# Builds a reconstruction plot with correlations between original and reconstructed genes
#'
#' @param plier The PLIER object.
#' @param Y The original gene expression matrix.
#' @param B The B matrix (optional). If not provided, the one carried by the PLIER object.
#' @param save The file path to save the plot.
#' @return The reconstruction plot.
reconstrPlot <- function(plier, Y, B = plier$B, save = NULL) {
  dt_cor <- extractCor(plier, Y, B)
  
  plot <- ggplot(dt_cor, aes(x = Correlation)) + 
    geom_density(fill = "blue", alpha = 0.4) + 
    xlim(-1, 1) +
    theme_bw() +
    labs(x = "Spearman Correlation", title = "Reconstruction Quality") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if (!is.null(save)) {
    ggsave(save, plot = plot, width=8, height=4, dpi = 300)    
  }
  
  return(plot)
}

##### USE OF LV DATASET

# Performs Mann-Whitney U test on latent variables (LVs) and adjusts p-values
#'
#' @param dt_lv A data.table contining the LV description of each sample.
#' @param dt_clin A data.table containing clinical data with a column that specifies the group of interest.
#' @param on A character string specifying the name of the column in `dt_clin` that indicates the group of interest (e.g., "Yes" or "No").
#' @return A data.table with p-values for each LV, adjusted p-values with BH correction, and regulation status.
lvsMannWhitney <- function(dt_lv, dt_clin, on) {
  # Remove rows with NA in the specified clinical column
  dt_clin <- dt_clin[!is.na(get(on))]
  
  # Get the 'Sample' identifiers for "Yes" and "No" cases
  yes_cases <- dt_clin[get(on) == "Yes"]$Sample
  no_cases <- dt_clin[get(on) == "No"]$Sample
  
  # Get the names of the latent variables (LVs)
  lvs_names <- setdiff(names(dt_lv), "Sample")
  
  # Initialize vectors to store p-values and regulation types
  wilc <- c()
  reg <- c()
  
  # Perform the Mann-Whitney U test for each LV
  for (lv in lvs_names) {
    # Get values for "Yes" and "No" cases for the current LV
    yes_vals <- dt_lv[Sample %in% yes_cases][[lv]]
    no_vals <- dt_lv[Sample %in% no_cases][[lv]]
    
    # Perform the Mann-Whitney U test
    w_pval <- wilcox.test(yes_vals, no_vals, alternative = "two.sided")$p.value
    
    # Determine if the LV is upregulated or downregulated
    reg_type <- ifelse(median(yes_vals) > median(no_vals), "Upregulated", "Downregulated")
    
    # Store the p-value and regulation type
    wilc <- c(wilc, w_pval)
    reg <- c(reg, reg_type)
  }
  
  # Create a data table with LV names and p-values
  dt_wilc <- data.table(`LV name` = lvs_names, `p-value` = wilc)
  
  # Adjust p-values using the BH method
  dt_wilc$`Adj. p-value` <- p.adjust(wilc, method = "BH")
  
  # Add regulation type to the data table
  dt_wilc$Regulation <- reg
  
  # Set regulation type to NA for non-significant adjusted p-values
  dt_wilc$Regulation[dt_wilc$`Adj. p-value` > 0.05] <- NA
  
  # Order the data table by p-value
  dt_wilc <- dt_wilc[order(`p-value`)]
  
  # Return the final data table
  return(dt_wilc)
}

#' Plot Significant Latent Variables (LVs)
#'
#' This function generates a heatmap showing the regulation (up or down) and significance (adjusted p-value) of latent variables (LVs).
#'
#' @param dt_signif A data.table containing the significant LVs, with columns `LV name`, `Regulation`, `Adj. p-value`, and `p-value`.
#' @param save Character. The file path to save the resulting plot, or NULL to not save the plot.
#'
#' @return A ggplot object showing a heatmap of LVs with their regulation (up or down) and significance.
#'
#' @details
#' This function processes the input data to categorize LVs by their regulation and adjusted p-value, then plots them using `ggplot2`.
#' The p-values are divided into three categories: p < 0.001, 0.001 <= p < 0.01, and 0.01 <= p < 0.05. The resulting plot can be optionally saved to a file.
plotSignificantLVs <- function(dt_signif, save = NULL) {
  
  # Create a simplified column for regulation (Up/Down)
  dt_signif[Regulation == "Upregulated", Reg := "Up"]
  dt_signif[Regulation == "Downregulated", Reg := "Down"]
  
  # Format p-values for display and fill color categories
  dt_signif$text <- ifelse(dt_signif$`Adj. p-value` < 0.001, "<0.001", round(dt_signif$`Adj. p-value`, 3))
  dt_signif[`Adj. p-value` < 0.001, fill := "p<0.001"]
  dt_signif[(`Adj. p-value` >= 0.001 & `Adj. p-value` < 0.01), fill := "0.001<=p<0.01"]
  dt_signif[(`Adj. p-value` >= 0.01 & `Adj. p-value` < 0.05), fill := "0.01<=p<0.05"]

  # Convert fill column to a factor with ordered levels
  dt_signif$fill <- as.factor(dt_signif$fill)
  dt_signif$fill <- factor(dt_signif$fill, levels = c("p<0.001", "0.001<=p<0.01", "0.01<=p<0.05"))

  # Order the LVs by significance and set the levels for plotting
  dt_signif <- dt_signif[order(-`p-value`)]
  dt_signif$`LV name` <- as.factor(dt_signif$`LV name`)
  dt_signif$`LV name` <- factor(dt_signif$`LV name`, levels = rev(dt_signif$`LV name`))

  # Additional columns for plotting aesthetics
  dt_signif$Dataset <- ""  # Empty dataset column for visual purposes
  dt_signif$classes <- paste(dt_signif$Reg, dt_signif$fill)

  # Define levels for the classes (Up/Down with p-value ranges)
  vals <- c("p<0.001", "0.001<=p<0.01", "0.01<=p<0.05")
  levels <- rev(c(paste("Down", vals), "Not significant", paste("Up", rev(vals))))
  dt_signif$classes <- factor(dt_signif$classes, levels = levels)

  # Generate the heatmap plot
  plot <- ggplot(dt_signif, aes(`LV name`, Dataset, fill = classes)) + 
    geom_tile(color = "black", alpha = 0.7) + 
    geom_text(aes(x = `LV name`, y = Dataset, label = text), color = "black", size = 4) + 
    theme_bw() + theme(panel.grid.major = element_blank()) + 
    theme(legend.position = 'left') + 
    scale_fill_brewer(type = "seq", palette = "PiYG", direction = -1, drop = FALSE) + 
    coord_fixed() + ylab("") + 
    ggtitle("Regulation and significance (adj. p-value)") + 
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    theme(legend.title = element_blank()) + 
    theme(axis.text = element_text(size = 12))

  # Save the plot if a file path is provided
  if (!is.null(save)) {
    ggsave(save, plot = plot, width = 10, height = 4, dpi = 300)
  }

  return(plot)
}


##### CONVERSION FROM RAW COUNTS TO TPM

#' Convert raw counts to TPM
#'
#' This function converts raw gene counts from STAR output (.tab files) to TPM (Transcripts Per Million).
#'
#' @param folder Character. The folder containing the .tab files with raw counts.
#' @param lengths Character. The path to the file containing gene lengths in base pairs.
#' @param save Character. The file path to save the resulting TPM data.table, or NULL to not save.
#'
#' @return A data.table with ensembl_gene_id and TPM values for each sample.
convertCountsToTpm <- function(folder, lengths = "genelength.txt", save = NULL) {
  
  # Prepare empty matrix to input raw counts from STAR output (.tab files)
  tab_files <- list.files(folder, pattern = ".tab", full.names = TRUE)
  n_genes <- nrow(read.table(tab_files[1], header = FALSE, row.names = 1, stringsAsFactors = FALSE, skip = 4))
  n_samples <- length(tab_files)
  ensembl_ids_raw <- rownames(read.table(tab_files[1], header = FALSE, row.names = 1, stringsAsFactors = FALSE, skip = 4))
  ensembl_ids <- grep("ENSG", ensembl_ids_raw, value = TRUE)
  matrix_raw_counts <- matrix(data = NA, ncol = n_samples, nrow = n_genes)
  
  # Load raw counts into the matrix
  pb <- txtProgressBar(min = 0, max = length(tab_files), initial = 0, char = "=", width = NA, style = 3)
  for (i in 1:length(tab_files)) { 
    setTxtProgressBar(pb, i)
    matrix_raw_counts[, i] <- read.table(tab_files[i], header = FALSE, row.names = 1, stringsAsFactors = FALSE, skip = 4)[, 1]
  }
  close(pb)
  colnames(matrix_raw_counts) <- sub(".*(SRR[0-9]+)\\.Reads.*", "\\1", tab_files)
  rownames(matrix_raw_counts) <- ensembl_ids
  
  # Calculate RPK (reads per kilobase) and TPM (transcripts per million)
  
  # Load gene length data
  gene_lengths <- read.table(file = lengths, sep = "\t", header = TRUE, row.names = 1)
  gene_lengths$length_kbp <- gene_lengths$Length / 1000
  matrix_rpk <- data.matrix(matrix_raw_counts / gene_lengths$length_kbp)
  matrix_tpm <- t(t(matrix_rpk) * 1e6 / colSums(matrix_rpk))
  
  # Verify and save TPM data
  col_sums <- colSums(matrix_tpm)  # All columns should sum to 1e06
  
  dt_tpm <- as.data.table(matrix_tpm, keep.rownames = "ensembl_gene_id")
  
  # Save the dataset if save path is provided
  if (!is.null(save)) {
    saveRDS(dt_tpm, save)
  }
  
  return(dt_tpm)
}

##### SANKEY PLOT FOR LVS

# These functions script create a Sankey diagram describing LV002 and LV270 in terms of association to signatures.
# It can be customized to plot the LVs of interest.
# You can plot one or more LVs.

#' Create Sankey input data
#'
#' This function creates input data for generating a Sankey diagram using a data.table with values related to different latent variables (LVs) and cell types.
#'
#' @param dt A data.table containing latent variables (LVs), cell types, and other relevant columns.
#' @param LVs_to_plot A vector of latent variables to include in the Sankey diagram.
#'
#' @return A list with two data.tables: 'nodes' (characterizing the nodes) and 'links' (characterizing the links between nodes).
createSankeyInput <- function(dt, LVs_to_plot) {
  
  # Preprocessing and LV selection
  dt[, U := U / sum(U), .(`LV name`)]
  dt[, LVsum := sum(U), .(`LV name`)]
  dt[, gsum := sum(U), .(`LV name`, `Cell Type`)]
  
  dt_nodes_full <- data.table()
  dt_links_full <- data.table()
  
  c <- 0
  for (LV_sel in LVs_to_plot) {
    
    # Formatting
    st <- paste(rep(" ", c), collapse = "")
    
    # Selection
    dt_red <- dt[`LV name` == LV_sel]
    
    # Creating dataset characterizing nodes
    nodes2 <- unique(dt_red[order(-gsum)]$`Cell Type`)
    nodes3 <- unique(dt_red[order(-gsum, -U)]$`Signature Name`)
    nodes <- c(LV_sel, nodes2, nodes3)
    dt_nodes <- data.table(name = nodes)
    dt_nodes$group <- "nodes"
    dt_nodes$name <- paste0(dt_nodes$name, st)
    
    # Creating dataset characterizing links
    dt_links <- rbind(
      dt_red[, .(value = sum(U)), .(source = paste0(`LV name`, st), target = paste0(`Cell Type`, st))],
      dt_red[, .(value = sum(U)), .(source = paste0(`Cell Type`, st), target = paste0(`Signature Name`, st))]
    )
    nodes <- dt_nodes$name
    dt_links$source2 <- match(dt_links$source, nodes) - 1 + nrow(dt_nodes_full)
    dt_links$target2 <- match(dt_links$target, nodes) - 1 + nrow(dt_nodes_full)
    dt_links$group <- "links"
    dt_links <- dt_links[, .(source = source2, target = target2, value, group)]
    
    # Binding to global datasets
    dt_nodes_full <- rbind(dt_nodes_full, dt_nodes)
    dt_links_full <- rbind(dt_links_full, dt_links)
    
    # Updating count
    c <- c + 1
  }
  
  # Creating object containing both nodes and links
  sankey_list <- list(nodes = dt_nodes_full, links = dt_links_full)
  
  return(sankey_list)
  
}


#' Plot Sankey diagram
#'
#' This function generates a Sankey diagram from latent variable (LV) data and atlas details.
#'
#' @param clier A PLEIR object.
#' @param atlas_details A data.table containing information about the cell types and signatures.
#' @param LVs_to_plot A vector of latent variables to include in the Sankey diagram.
#' @param fontSize Numeric. Font size for the node labels. Default is 36.
#' @param nodeWidth Numeric. Width of the nodes. Default is 40.
#' @param nodePadding Numeric. Padding between the nodes. Default is 35.
#' @param width Numeric. Width of the Sankey diagram. Default is 1750.
#' @param height Numeric. Height of the Sankey diagram. Default is 1200.
#' @param margin.left Numeric. Left margin for the Sankey diagram. Default is 1200.
#' @param iterations Numeric. Number of iterations for optimizing node positions. Default is 3.
#'
#' @return A Sankey network plot.
plotSankey <- function(clier, atlas_details, LVs_to_plot, fontSize = 36, nodeWidth = 40, nodePadding = 35, 
                       width = 1750, height = 1200, margin.left = 1200, iterations = 3) {
  
  dt_dict <- reshapeUInfo(clier)
  dt <- merge(dt_dict, atlas_details[, .(Signature = `Signature Name`, `Cell Type`)])
  dt <- dt[, .(`LV name`, `Signature Name` = Signature, U, AUC, FDR, `p-value`, `Cell Type`)][order(`LV name`, -U)]
  dt[`Cell Type` == "Proximal tubule", `Cell Type` := "PT"] # Shortening cell type names
  dt[`Cell Type` == "Not associated", `Cell Type` := "Not ass."]
  dt[`Cell Type` == "Loop of Henle", `Cell Type` := "LOH"]
  
  # Generate Sankey input data
  sankeyInput <- createSankeyInput(dt, LVs_to_plot)
  
  # Define color scheme for nodes and links
  my_color <- 'd3.scaleOrdinal().domain(["nodes", "links"]).range(["black", "#D3D3D3"])'
  
  # Plot the Sankey diagram
  sn <- sankeyNetwork(Links = sankeyInput$links, Nodes = sankeyInput$nodes,
                      Source = "source", Target = "target", Value = "value", NodeID = "name",
                      fontSize = fontSize, nodeWidth = nodeWidth, colourScale = my_color,
                      NodeGroup = "group", LinkGroup = "group",
                      nodePadding = nodePadding, iterations = iterations,
                      width = width, height = height,
                      # Provide space for newly aligned labels
                      margin = list("left" = margin.left)
  )
  
  # This step is necessary to adjust node text alignment
  onRender(
    sn,
    '
  function(el, x) {
    // Select all node text
    var node_text = d3.select(el)
      .selectAll(".node text")
      // Adjust to match new alignment
      .attr("x", 6 + x.options.nodeWidth)
      .attr("text-anchor", "start");
  }
  '
  )
}


#### AGREEMENT BETWEEN TWO PLIERS

#' Get correlation vector
#'
#' This function calculates a vector of correlation values following the methodology by Taroni et al.
#'
#' @param Z.ref A matrix representing the reference latent variable (LV) loadings (Z matrix) from the reference PLIER model.
#' @param Z.test A matrix representing the test latent variable (LV) loadings (Z matrix) from the test PLIER model.
#' @param B.ref A matrix representing the reference LV expression matrix (B matrix) from the reference PLIER model.
#' @param B.test A matrix representing the test LV expression matrix (B matrix) from the test PLIER model.
#'
#' @return A numeric vector containing the correlation values between the reference and test models.
#'
#' @details
#' The function computes the Spearman correlation between the latent variable loadings of the reference and test models, selecting the most correlated pairs between columns of the Z matrices. It then computes the correlations between the B matrix values corresponding to the selected columns.
getCorrVector <- function(Z.ref, Z.test, B.ref, B.test) {

  # Intersect the genes between reference and test models
  genes <- intersect(rownames(Z.ref), rownames(Z.test))
  Z.ref <- Z.ref[genes, ]
  Z.test <- Z.test[genes, ]
  print(length(genes))

  # Get the number of columns in Z.ref and Z.test
  n_a <- ncol(Z.ref)
  n_b <- ncol(Z.test)

  # Initialize a vector to store the index of the most correlated column of Z.test for each column of Z.ref
  best_correlated_columns <- numeric(n_a)

  # Loop through each column of Z.ref
  for (i in 1:n_a) {
    #print(i)

    # Get the i-th column of Z.ref
    col_a <- Z.ref[, i]

    # Initialize a vector to store correlation values for each column of Z.test
    correlations <- numeric(n_b)

    # Loop through each column of Z.test
    for (j in 1:n_b) {

      # Get the j-th column of Z.test
      col_b <- Z.test[, j]

      # Compute the Spearman correlation between col_a and col_b
      correlations[j] <- cor(col_a, col_b, method = "spearman")
    }

    # Find the index of the maximum correlation value
    best_correlated_columns[i] <- which.max(correlations)
  }

  # Initialize a vector to store the correlation between the reference and test B matrices
  corr_vector <- numeric(n_a)

  # Loop through each column of Z.ref
  for (i in 1:n_a) {
    # Compute the correlation between the i-th row of B.ref and the corresponding row of B.test
    corr_vector[i] <- cor(B.ref[i, ], B.test[best_correlated_columns[i], ], method = "spearman")
  }

  return(corr_vector)
}

#' Check agreement between reference and test PLIER models
#'
#' This function checks the agreement between the reference PLIER model and a test PLIER model by calculating the correlation between latent variable (LV) and gene loadings.
#'
#' @param plier.ref A reference PLIER object.
#' @param plier.test A test PLIER object.
#'
#' @return A ggplot object showing the density plot of correlation values between the reference and test models.
#'
#' @details
#' The function extracts the relevant latent variables (LVs) from both the reference and test models, computes the correlation between their loadings, and generates a density plot of the correlation values.
checkAgreement <- function(plier.ref, plier.test) {

  # Reference PLIER
  Z.ref <- plier.ref$Z
  B.ref <- plier.ref$B
  interp <- unique(reshapeUInfo(plier.ref)$`LV name`)
  rm(plier.ref)

  Z.ref <- Z.ref[, interp]
  B.ref <- B.ref[interp, ]

  # Test PLIER
  Z.test <- plier.test$Z
  B.test <- plier.test$B
  interp <- unique(reshapeUInfo(plier.test)$`LV name`)
  rm(plier.test)

  Z.test <- Z.test[, interp]
  B.test <- B.test[interp, ]

  # Check agreement by calculating correlation vector
  vec <- getCorrVector(Z.ref, Z.test, B.ref, B.test)
  dt.plot <- data.table(Model = "Test model", Values = vec)

  # Generate density plot
  p <- ggplot(dt.plot, aes(x = Values, fill = Model)) +
    geom_density() +
    theme_bw() +
    xlim(-1, 1)

  return(p)
}

