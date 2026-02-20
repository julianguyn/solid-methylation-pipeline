# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(minfi)
    library(ComplexHeatmap)
    library(circlize)
    library(umap)
    library(ggplot2)
})

source("workflow/scripts/utils/palettes.R")

set.seed(101)

###########################################################
# Load in data
###########################################################

message("Loading in data")

# load in RGSets and targets
load("data/procdata/RGSets_targets.RData")
targets$match <- gsub(".*/", "", targets$Basename)

# get beta values
bVals <- fread("data/results/data/bVals.tsv", data.table = FALSE)
rownames(bVals) <- bVals$V1
bVals$V1 <- NULL

# get m values
mVals <- fread("data/results/data/mVals.tsv", data.table = FALSE)
rownames(mVals) <- mVals$V1
mVals$V1 <- NULL

###########################################################
# Format dataframes for plotting
###########################################################

message("Identifying top 5000 CpGs")

# keep only top 5000 cpgs
cpg_var <- apply(mVals, 1, var)
top_cpgs <- names(sort(cpg_var, decreasing = TRUE)[1:5000])
bVals <- bVals[top_cpgs, ]
mVals <- mVals[top_cpgs, ]

# transpose matrices
bVals <- t(bVals) |> as.data.frame()
mVals <- t(mVals) |> as.data.frame()

###########################################################
# Heatmap clustering
###########################################################

message("Performing heatmap clustering")

col_fun <- colorRamp2(c(0, 0.5, 1), c("#457B9D", "white", "#9D3535"))

# colour palettes
pfs_pal <- colorRamp2(c(0, max(targets$PFS_Days, na.rm = TRUE)), c("white", "#E76F51"))
os_pal <- colorRamp2(c(0, max(targets$OS_days, na.rm = TRUE)), c("white", "#6A0572"))

ra <- rowAnnotation(
    RECIST = targets$Response_RANO[match(rownames(bVals), targets$Subject)],
    PFS_Days = targets$PFS_Days[match(rownames(bVals), targets$Subject)],
    OS_Days = targets$OS_days[match(rownames(bVals), targets$Subject)],
    col = list(
        RECIST = mrecist_pal,
        PFS_Days = pfs_pal,
        OS_Days = os_pal
    )
)

filename <- "data/results/figures/bvals_heatmap.pdf"
pdf(filename, width = 20, height = 9)
Heatmap(
    bVals,
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_dend = FALSE,
    show_column_names = FALSE,
    right_annotation = ra,
    heatmap_legend_param = list(title = "Beta\nValue")
)
dev.off()

###########################################################
# UMAP
###########################################################

message("Performing UMAP")

# perform umap
umap_df <- umap(mVals)$layout |> as.data.frame()
colnames(umap_df) <- paste0("UMAP", 1:2)
targets <- targets[match(rownames(umap_df), targets$Subject),]
umap_df$RECIST <- targets$Response_RANO
umap_df$PFS_Days <- log2(targets$PFS_Days)
umap_df$OS_Days <- log2(targets$OS_days)

# plot umaps
filename <- "data/results/figures/umap_RECIST.png"
png(filename, width = 5, height = 4, res = 600, units = "in")
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, fill = RECIST)) +
    geom_point(shape = 21, size = 3) +
    scale_fill_manual(values = mrecist_pal) +
    theme_minimal() +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.key.size = unit(0.5, 'cm')
    )
dev.off()

filename <- "data/results/figures/umap_PFS.png"
png(filename, width = 5, height = 4, res = 600, units = "in")
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, fill = PFS_Days)) +
    geom_point(shape = 21, size = 3) +
    scale_fill_gradient(high = "#1D3557", low = "#A8DADC") +
    theme_minimal() +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.key.size = unit(0.5, 'cm')
    )
dev.off()

filename <- "data/results/figures/umap_OS.png"
png(filename, width = 5, height = 4, res = 600, units = "in")
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, fill = OS_Days)) +
    geom_point(shape = 21, size = 3) +
    scale_fill_gradient(high = "#1D3557", low = "#A8DADC") +
    theme_minimal() +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.key.size = unit(0.5, 'cm')
    )
dev.off()

print("done")