
# NOTE
# this script contains the code presented in the manuscript
# it can be entirely executed with the three following inputs (that are available on this Git):
# the file "kidney_atlas.rds"
# the file "DKD_tpm.rds"
# the file "DKD_clin.rds"

# Please note that, if you are interested in executing only specific parts, you will need intermediate files.
# You can find them here: https://drive.switch.ch/index.php/s/OpvMh1vGRgRmKKf


### TRAINING CLIER

# Load auxiliary functions 
source("aux_functions.R")

# Download TCGA dataset
dt_tcga <- downloadTcga(save = "tcga.rds")

# Read the signatures atlas
atlas <- readRDS("kidney_atlas_matrix.rds")  # customize 

# Preprocess data
res_prep <- trainPreprocess(dt_tcga, atlas)

# Train PLIER
clier <- PLIER(res_prep$Y, atlas, k = 700, allGenes = TRUE, trace = FALSE, scale = FALSE, frac = 0.5)

# Add preprocessing info to PLIER object
clier[["trainStats"]] <- res_prep$dt_stats

# Rename LVs
clier <- renameLvs(clier)

# Reshape information
dt_info <- reshapeUInfo(clier, save = "associations.xlsx")
head(dt_info)

# Reconstruction plot
reconstrPlot(clier, Y = res_prep$Y, B = clier$B, save = "train_rec.png")

# Save PLIER object
saveRDS(clier, "kclier.rds")  # customize


### APPLYING CLIER TO NEW DATASET

# Load auxiliary functions 
#source("aux_functions.R"). # Uncomment if you only execute the code from here on

# Load previously trained PLIER
clier <- readRDS("kclier.rds")  # customize

# Read new dataset
dt_new <- readRDS("./DKD_tpm.rds")  # customize

# Preprocess the new dataset
Y_new <- testPreprocess(dt_new, clier$trainStats)

# Apply PLIER to new data
B_new <- getB(clier, Y_new, save = "DKD_B.rds")  # customize

# Reconstruction plot
reconstrPlot(clier, Y = Y_new, B = B_new, save = "new_dataset_rec.png")


### USING THE NEW LV DATASET

# Load auxiliary functions 
#source("aux_functions.R"). # Uncomment if you only execute the code from here on

# Loading clinical dataset
dt_clin <- readRDS("DKD_clin.rds")  # customize

# Loading LV matrix and convert it to a data.table
B <- readRDS("DKD_B.rds")  # customize
dt_lv <- as.data.table(t(B), keep.rownames = "Sample")
  
# restrict to interpretable LVs
dt_info <- read_excel("associations.xlsx")
interp_lvs <- unique(dt_info$`LV name`)
dt_lv <- dt_lv[, c("Sample", interp_lvs), with = FALSE]

# Conduct Mann-Whitney tests
dt_mw <- lvsMannWhitney(dt_lv, dt_clin, on = "Fibrosis")  # customize
head(dt_mw, 10)

# Plot the results of Mann-Whitney tests
plotSignificantLVs(head(dt_mw, 10), save="significant_LVs.png")

# Inspect the structure of the most significant LVs
clier <- readRDS("kclier.rds") # customize
atlas_details <- read_excel("kidney_atlas_info.xlsx") %>% as.data.table() # customize
plotSankey(clier, atlas_details, LVs_to_plot=c("LV002","LV270")) # customize