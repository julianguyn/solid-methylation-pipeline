# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(limma)
    library(ggplot2)
    library(ggpubr)
    library(RColorBrewer)
    library(viridis)
    library(GenomicRanges)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(rGREAT)
    library(clusterProfiler)
    library(org.Hs.eg.db)
})

source("workflow/scripts/utils/palettes.R")

set.seed(101)

###########################################################
# Load in data
###########################################################

mat <- readRDS("data/rawdata/cholangio/MeDIP_Solid_cholangiocarcinoma_log2CPM.rds")

pheno1 <- fread("data/rawdata/cholangio/Pheno_cholangio_all_IDHwt_Vs_IDHmut.txt") #64
pheno2 <- fread("data/rawdata/cholangio/Pheno_cholangio_baseline_IDHwt_Vs_IDHmut.txt") #16

# keep only PD & SD
pheno <- pheno1[pheno1$response != "baseline",] # 38

###########################################################
# Quick PCA
###########################################################

# rank date for each sample
pheno$date <- as.Date(pheno$blood_collection_date, format = "%Y %b %d")
pheno <- pheno[order(pheno$solid_patient, pheno$date), ]
pheno$timepoint <- unlist(lapply(as.numeric(table(pheno$solid_patient)), seq_len))
pheno$timepoint <- factor(paste0("T", pheno$timepoint), levels = paste0("T", 1:12))

# perform PCA
mat <- mat[,match(pheno$MeDIP_ID, colnames(mat))] |> as.matrix()
pca_res <- prcomp(t(mat))$x[,1:2] |> as.data.frame()

p1 <- ggplot(pca_res, aes(x = PC1, y = PC2, color = pheno$solid_patient, shape = pheno$response)) +
    geom_point(size = 3) +
    scale_color_brewer(palette = "Set3") +
    theme_bw()

p2 <- ggplot(pca_res, aes(x = PC1, y = PC2, color = pheno$timepoint, shape = pheno$response)) +
    geom_point(size = 3) +
    scale_color_viridis_d(option = "rocket", begin = 0.2) +
    theme_bw()

p <- ggarrange(p1, p2, ncol=2)
ggsave("data/results/figures/cholangio/pca_dmrs_medip.png", p, width = 11.5, height = 4.5)

###########################################################
# Run linear model (SD-PD)
###########################################################

# set variables
labels <- factor(pheno$response)

# model
design <- model.matrix(~ labels) # 0 DMRs with FDR<0.05
fit <- lmFit(mat, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "labelsSD", number = Inf)
table(res$adj.P.Val < 0.1)

###########################################################
# Run linear model (patient + SD-PD)
###########################################################

# set variables
labels <- factor(pheno$response)
patientid <- factor(pheno$solid_patient)

# model
design <- model.matrix(~ patientid + labels) # 1 DMR with FDR<0.05
fit <- lmFit(mat, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "labelsSD", number = Inf)
table(res$adj.P.Val < 0.1)

###########################################################
# Run linear model (corfit)
###########################################################

# design matrix
design <- model.matrix(~ labels)

# estimate correlation within patients
corfit <- duplicateCorrelation(mat, design, block = patientid)

# model (0 DMRs)
fit <- lmFit(mat, design, block = patientid, correlation = corfit$consensus)
fit <- eBayes(fit)
res <- topTable(fit, coef = "labelsSD", number = Inf)
table(res$adj.P.Val < 0.1)

###########################################################
# Keep only first SD and first PD
###########################################################

table(pheno$response) # PD: 21, SD: 17

# keep first SDs
sd_pheno <- pheno[pheno$response == "SD",]
idx <- cumsum(as.numeric(table(sd_pheno$solid_patient)))
idx <- c(1, idx+1)[-c(length(idx)+1)]
sd_pheno <- sd_pheno[idx,]

# keep first PDs
pd_pheno <- pheno[pheno$response == "PD",]
idx <- cumsum(as.numeric(table(pd_pheno$solid_patient)))
idx <- c(1, idx+1)[-c(length(idx)+1)]
pd_pheno <- pd_pheno[idx,]

###########################################################
# Check PCA again
###########################################################

# subset matrix
subset_pheno <- rbind(sd_pheno, pd_pheno)
subset_mat <- mat[,match(subset_pheno$MeDIP_ID, colnames(mat))]

# add labels to PCA
pca_res$patientid <- pheno$solid_patient
pca_res$MeDIP_ID <- pheno$MeDIP_ID
pca_res$response <- pheno$response
pca_res$label <- ifelse(
    pca_res$MeDIP_ID %in% subset_pheno$MeDIP_ID,
    ifelse(pca_res$response == "SD", "SD", "PD"),
    "Not\nIncluded"
)

# plot pca
p3 <- ggplot() +
    geom_point(data = pca_res[pca_res$label == "Not\nIncluded",], aes(x = PC1, y = PC2), color = "gray", size = 3) +
    geom_point(data = pca_res[pca_res$label != "Not\nIncluded",], aes(x = PC1, y = PC2, color = patientid, shape = label), size = 3) +
    scale_color_brewer(palette = "Set3") +
    theme_bw()

p <- ggarrange(p1, p3, ncol=2)
ggsave("data/results/figures/cholangio/pca_dmrs_medip_subset.png", p, width = 11.5, height = 4.5)


###########################################################
# Run linear model (SD-PD)
###########################################################

# set variables
labels <- factor(subset_pheno$response)

# model
design <- model.matrix(~ labels) # 0 DMRs with FDR<0.05
fit <- lmFit(subset_mat, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "labelsSD", number = Inf)
table(res$adj.P.Val < 0.1) # 17 DMRs
table(res$adj.P.Val < 0.05) #17

###########################################################
# Run linear model (patient + SD-PD)
###########################################################

# set variables
labels <- factor(subset_pheno$response)
patientid <- factor(subset_pheno$solid_patient)

# model
design <- model.matrix(~ patientid + labels) # 1 DMR with FDR<0.05
fit <- lmFit(subset_mat, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "labelsSD", number = Inf)
table(res$adj.P.Val < 0.1) #479 DMRs
table(res$adj.P.Val < 0.05) #478

###########################################################
# Volcano plot
###########################################################

ggplot() +
    geom_point(data = res, aes(x = logFC, y = -log10(adj.P.Val)), color = "gray", size = 2) +
    geom_point(data = res[res$logFC > 2 & res$adj.P.Val < 0.05,], aes(x = logFC, y = -log10(adj.P.Val)), color = "#B55965", size = 2) +
    geom_point(data = res[res$logFC < -2 & res$adj.P.Val < 0.05,], aes(x = logFC, y = -log10(adj.P.Val)), color = "#777DA7", size = 2) +
    theme_bw()

###########################################################
# DMR analysis - region annotation
###########################################################

sig_res <- res[abs(res$logFC) > 2 & res$adj.P.Val < 0.05,]
hyper <- sig_res[sig_res$logFC > 2,]
hypo <- sig_res[sig_res$logFC < -2,]

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

###########################################################
# Create genomic ranges
###########################################################

# get coordinates
hyper_coords <- do.call(rbind, strsplit(rownames(hyper), "\\.")) |> as.data.frame()
hypo_coords <- do.call(rbind, strsplit(rownames(hypo), "\\.")) |> as.data.frame()

# make genomic ranges object
hyper_gr <- GRanges(seqnames = hyper_coords$V1, ranges = IRanges(as.numeric(hyper_coords$V2), as.numeric(hyper_coords$V3)))
hypo_gr <- GRanges(seqnames = hypo_coords$V1, ranges = IRanges(as.numeric(hypo_coords$V2), as.numeric(hypo_coords$V3)))

###########################################################
# Annotate genomic regions
###########################################################

# annotate regions
hyper_anno <- annotatePeak(hyper_gr, tssRegion=c(-250, 250), TxDb=txdb, annoDb="org.Hs.eg.db")
hypo_anno <- annotatePeak(hypo_gr, tssRegion=c(-250, 250), TxDb=txdb, annoDb="org.Hs.eg.db")

# create toPlot
hyper <- hyper_anno@annoStat
hypo <- hypo_anno@annoStat
hyper$Label <- "Hyper"
hypo$Label <- "Hypo"
toPlot <- rbind(hyper, hypo)

# plot gene features
filename <- paste0("data/results/figures/cholangio/sdpd_peakanno.png")
png(filename, width = 5, height = 5, res = 600, units = "in")
ggplot(toPlot, aes(fill=Feature, y=Frequency, x=Label)) + 
    geom_bar(position="fill", stat="identity", color = "black") +
    scale_fill_manual(values = genefeat_pal) +
    theme_minimal() + 
    labs(y = "Percentage (%)")
dev.off()

###########################################################
# Keep promoter regions
###########################################################

# 11 hypermethylated promoters
hyper_df <- as.data.frame(hyper_anno)
hyper_promoter <- hyper_df[grep("Promoter", hyper_df$annotation), ]

# 38 hypomethylated promoters
hypo_df <- as.data.frame(hypo_anno)
hypo_promoter <- hypo_df[grep("Promoter", hypo_df$annotation), ]

# make genomic ranges object
hyper_gr_promoter <- GRanges(seqnames = hyper_promoter$seqnames, ranges = IRanges(hyper_promoter$start, hyper_promoter$end))
hypo_gr_promoter <- GRanges(seqnames = hypo_promoter$seqnames, ranges = IRanges(hypo_promoter$start, hypo_promoter$end))

###########################################################
# Run GO analysis
###########################################################

# run GO
run_go <- function(df, label) {

    genes <- unique(df$geneId)
    ontologies <- c("BP", "MF", "CC")

    # run enrichGO
    go_list <- lapply(ontologies, function(o) {
        enrichGO(
            gene = genes,
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = o,
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.05,
            readable = TRUE
        ) |> as.data.frame()
    })
    names(go_list) <- ontologies

    # bind results
    bp <- go_list$BP
    mf <- go_list$MF
    cc <- go_list$CC

    if (nrow(bp) > 0) bp$Label <- "BP"
    if (nrow(mf) > 0) mf$Label <- "MF"
    if (nrow(cc) > 0) cc$Label <- "CC"

    res <- rbind(bp, mf, cc)
    filename <- paste0("data/results/data/cholangio/sdpd_", label, "_GO.tsv")
    write.table(res, file = filename, quote = F, sep = "\t", col.names = T, row.names = F)
    print(nrow(res))
    return(res)
}

# genes from all regions
hyper_go_all <- run_go(hyper_df, "hyper_all") # 2 results
hypo_go_all <- run_go(hypo_df, "hypo_all") # 0 results

# genes from promoter regions only
hyper_go_promoter <- run_go(hyper_promoter, "hyper_promoter") # no results
hypo_go_promoter <- run_go(hypo_promoter, "hypo_promoter") # no results


###########################################################
# GREAT analysis
###########################################################

# helper function to run GREAT anlaysis
run_great <- function(gr, label) {
    job <- submitGreatJob(gr, species = "hg19", genome = "hg19", help = FALSE)
    tbl <- getEnrichmentTables(job)

    # save results
    mf <- as.data.frame(tbl[1])
    bp <- as.data.frame(tbl[2])
    cc <- as.data.frame(tbl[3])

    # save column names
    cols <- gsub("GO.Molecular.Function.", "", colnames(mf))
    colnames(mf) <- colnames(bp) <- colnames(cc) <- cols

    mf$Label <- "MF"
    bp$Label <- "BP"
    cc$Label <- "CC"

    # combine results and filter
    res <- rbind(bp, mf, cc)
    res <- res[res$Hyper_Adjp_BH < 0.05,]
    res <- res[order(res$Hyper_Adjp_BH),]

    filename <- paste0("data/results/data/cholangio/sdpd_", label, ".tsv")
    write.table(res, file = filename, quote = F, sep = "\t", col.names = T, row.names = F)
    print(nrow(res))
    return(res)
}

hyper_res <- run_great(hyper_gr, "hyper_GREAT") # 0
#hyper_res_promoter <- run_great(hyper_gr_promoter, "hyper_promoter_GREAT")
hypo_res <- run_great(hypo_gr, "hypo_GREAT") # 51 all
hypo_res_promoter <- run_great(hypo_gr_promoter, "hypo_promoter_GREAT") #0

###########################################################
# Plot GREAT analysis
###########################################################

# helper function to plot GREAT analysis results
plot_GREAT <- function(great, label) {

    great <- great[great$Hyper_Total_Genes > 10,]
    great <- great[order(great$Hyper_Fold_Enrichment, decreasing = TRUE),]
    toPlot <- great[1:30,]
    toPlot$name <- factor(toPlot$name, levels=rev(toPlot$name))

    p <- ggplot(toPlot, aes(x = name, y = Hyper_Fold_Enrichment, color = Hyper_Adjp_BH)) +
        geom_point(aes(size = Hyper_Total_Genes, shape = Label)) +
        geom_segment(aes(x = name, xend = name, y = 0, yend = Hyper_Fold_Enrichment, color = Hyper_Adjp_BH), size = 1) +
        coord_cartesian(clip = "off") + 
        coord_flip() +
        scale_alpha(range = c(1, 0.2)) +
        guides(
            shape = guide_legend(title = "Ontology", override.aes = list(size = 4)),
            color = guide_legend(override.aes = list(shape = 19, size = 4))
        ) +
        scale_color_viridis_c(option = "rocket", direction = -1, begin = 0.15, end = 0.85) +
        theme_classic() + 
        theme(
            legend.key.size = unit(0.5, 'cm'), text = element_text(size = 10),
            legend.position = c(0.95, 0.05),   # bottom-right
            legend.justification = c("right", "bottom"),
            legend.background = element_rect(fill = alpha('white', 0.7), color = NA),
            plot.title = element_text(hjust = 0.5, size = 12), 
            panel.border = element_rect(color = "black", fill = NA, size = 0.5)
        ) + 
        labs(x = "GO Term", y = "Fold Enrichment", size = "Genes", color = "Adjusted\nP-Value", title = label)

    filename <- paste0("data/results/figures/cholangio/sdpd_", label, "_GREAT.png")
    png(filename, width = 8, height = 7, res = 600, units = "in")
    print(p)
    dev.off()
}

plot_GREAT(hypo_res, "hypomethylation")
