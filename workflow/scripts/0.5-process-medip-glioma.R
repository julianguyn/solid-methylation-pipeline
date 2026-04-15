# subset medip for glioma
# messy code!

suppressPackageStartupMessages({
    library(data.table)
    library(readxl)
})

###########################################################
# Load in data
###########################################################

# load in medip meta
medip_meta <- read_excel("metadata/solid_cfDNA_leftover_with_medip.xlsx", sheet = 1)

# glioma solid metadata
glioma_meta <- read.csv("metadata/clin_updated.csv")

# load in all medip
medip <- readRDS("data/rawdata/MeDIP_Solid_all_log2CPM.rds")

###########################################################
# Format metadata
###########################################################

# match solid IDs in medip metadata
for (i in 1:nrow(medip_meta)) {
    medip_meta$solid_name[i] <- ifelse(
        nchar(medip_meta$solid_name[i]) == 1,
        paste0("SOLID-00", medip_meta$solid_name[i]),
        paste0("SOLID-0", medip_meta$solid_name[i])
    )
}

# add medip availability to clin_updated
glioma_meta$MeDIP <- "No"
for (i in 1:nrow(glioma_meta)) {
    if (glioma_meta$Subject[i] %in% medip_meta$solid_name) glioma_meta$MeDIP[i] <- "Yes"
}
# send to Farnoosh
write.csv(glioma_meta, file = "metadata/clin_updated_medip.csv", quote = FALSE, row.names = FALSE)

# glioma patients with medip
have_medip <- glioma_meta$Subject[glioma_meta$MeDIP == "Yes"]

###########################################################
# Subset medip_meta and medip for glioma samples
###########################################################

# keep glioma samples
glioma_medip_meta <- medip_meta[medip_meta$solid_name %in% have_medip,] |> as.data.frame()

# screen samples
glioma_screen <- glioma_medip_meta[glioma_medip_meta$time == "screen", c(1:2, 5)] |> unique()

# subset medip
medip <- medip[,colnames(medip) %in% glioma_screen$MeDIP_ID]
colnames(medip) <- glioma_screen$solid_name[match(colnames(medip), glioma_screen$MeDIP_ID)]

save(medip, glioma_screen, file = "data/procdata/glioma_medip.RData")

# DROPPED

###########################################################
# Old code (not used, keep only screens)
###########################################################

# get cycle information
glioma_medip_meta$time[glioma_medip_meta$time == "screen"] <- "1d1"
glioma_medip_meta <- glioma_medip_meta[-which(glioma_medip_meta$time %in% c("eot", "fup")),]
glioma_medip_meta$cycle <- sub("c", "", sub("d.*", "", glioma_medip_meta$time))
glioma_medip_meta <- glioma_medip_meta[order(glioma_medip_meta$solid_name, as.numeric(glioma_medip_meta$cycle)),]

# take screen or first cycle
res <- data.frame(matrix(nrow=0, ncol=2))
for (patient in unique(glioma_medip_meta$solid_name)) {
    subset <- glioma_medip_meta[glioma_medip_meta$solid_name == patient,]
}

medip <- medip[,colnames(medip) %in% glioma_medip_meta$MeDIP_ID]



##### try mapping all medip

to_map <- colnames(medip) #174


# cholangio 
cho <- fread("data/rawdata/cholangio/Pheno_cholangio_all_IDHwt_Vs_IDHmut.txt")
cho_in_map <- cho$MeDIP_ID[which(cho$MeDIP_ID %in% to_map)] #56
num_patients <- length(unique(cho$solid_patient[cho$MeDIP_ID %in% to_map])) # 10

# cholangio not in medip
baseline <- cho$MeDIP_ID[-which(cho$MeDIP_ID %in% to_map)] # 8 baseline samples


remaining_to_map <- to_map[-which(to_map %in% cho_in_map)] #118

# glioma
glio <- unique(glioma_medip_meta$MeDIP_ID) # 68 samples across 16 patients



# missing 
missing_map <- remaining_to_map[-which(remaining_to_map %in% glio)]
