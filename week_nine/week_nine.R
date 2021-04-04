##### Week 9 #####

# Preliminaries
library(tidyverse)
library(DESeq2)
library(EnsDb.Hsapiens.v86)

## Where to get the data when recount2 isn't enough? ##

# 1. GREIN: http://www.ilincs.org/apps/grein/
# # a. Find a study. Download the count table.
# # b. Read into R and prepare for DESeq2
counts <- read_csv('GSE101770_GeneLevel_Raw_data.csv')
# Get the gene IDs as the rownames
tidy_counts <- counts %>%
  select(-gene_symbol) %>%
  column_to_rownames(var = "X1")
# Generate metadata colData object
colData <- data.frame(
  "sampleID" = colnames(tidy_counts),
  "condition" = c(rep("organoid", 2), 'cells',
                  rep('organoid', 2), 'cells',
                  rep('tissue', 2), 'cells'),
  row.names = colnames(tidy_counts)
)
# Generate DDS -- what went wrong??
dds <- DESeqDataSetFromMatrix(countData = tidy_counts,
                              colData = colData,
                              design = ~condition)
# Convert from decimals to integers.
# DO NOT DO THIS UNLESS YOU KNOW THE DATA SOURCE.
# YOU SHOULD NEVER DO THIS WITH NORMALIZED COUNTS!
rounded_counts <- apply(tidy_counts, MARGIN = 1:2, FUN = round)
dds <- DESeqDataSetFromMatrix(countData = rounded_counts,
                              colData = colData,
                              design = ~condition)
rld <- rlog(dds)
plotPCA(rld)

# 2. GDC: https://portal.gdc.cancer.gov/
# a. Download the ACC data from the GDC

untar("gdc_download_20201027_201844.644029.tar.gz", exdir = "gdc_acc_cancer")
untar("clinical.cart.2020-10-27 (1).tar.gz")
sample_sheet <- read_tsv("gdc_sample_sheet.2020-10-27 (1).tsv")
clinical <- read_tsv('clinical.tsv')
sampleTablejoin <- inner_join(sample_sheet, clinical,
                              by = c("Case ID" = "case_submitter_id"))
sampleTable <- sampleTablejoin %>%
  select(`Sample ID`, `File Name`, ajcc_pathologic_t) %>%
  distinct(`Sample ID`, .keep_all = TRUE) %>%
  dplyr::rename(condition = ajcc_pathologic_t)

DESeqDataSetFromHTSeqCount(sampleTable = as.data.frame(sampleTable),
                           directory = "gdc_acc_cancer",
                           design = ~condition)

sampleTable <- sampleTablejoin %>%
  mutate(`File Name` = paste0(`File ID`, "/", `File Name`)) %>%
  select(`Sample ID`, `File Name`, ajcc_pathologic_t) %>%
  distinct(`Sample ID`, .keep_all = TRUE) %>%
  dplyr::rename(condition = ajcc_pathologic_t)

dds <- DESeqDataSetFromHTSeqCount(sampleTable = as.data.frame(sampleTable),
                                  directory = "gdc_acc_cancer",
                                  design = ~condition)

vsd <- vst(dds)
plotPCA(vsd)

dds <- DESeq(dds)
resdf <- as.data.frame(results(dds, contrast = c("condition", "T4", "T1")))
# -- convert ENSG to gene symbol
ens2sym <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = keys(EnsDb.Hsapiens.v86),
                                 columns = c("SYMBOL"))
resdf <- resdf %>%
  rownames_to_column() %>%
  mutate(GENEID = gsub(rowname, pattern = "\\..+", replacement = "")) %>%
  dplyr::select(-rowname) %>%
  inner_join(y = ens2sym, by = "GENEID")
View(resdf)
# Write over-express to CSV
resdf %>%
  dplyr::filter(padj < .05 & log2FoldChange > 2) %>%
  write_csv(file = "T4_vs_T1.overexp.csv")
# Get the ranking score
resdf %>%
  dplyr::mutate(GSEA = -log10(padj) * sign(log2FoldChange)) %>%
  write_csv(file = "T4_vs_T1.GSEA_ranks.csv")



# 3. Synapse: https://www.synapse.org/
# 4. Encode: https://www.encodeproject.org/
# 5. ARCHS4: https://amp.pharm.mssm.edu/archs4/
# 6. GEO/SRA: https://www.ncbi.nlm.nih.gov/geo/





