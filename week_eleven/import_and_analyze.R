#### Week 11 RScript ####

### Import the data from salmon quant to DESeq2 using tximport ###

# Get tximport if you don't already have it installed

BiocManager::install("tximport")

# Load libraries
library(tximport)
library(DESeq2)
library(EnsDb.Hsapiens.v86)
library(tidyverse)

## Get the mapping from transcript IDs to gene symbols ##

# What are the columns in the database?
columns(EnsDb.Hsapiens.v86)

# Get the TXID and SYMBOL columns for all entries in database
tx2gene <- AnnotationDbi::select(EnsDb.Hsapiens.v86, 
                                   keys = keys(EnsDb.Hsapiens.v86),
                                   columns = c('TXID', 'SYMBOL'))

# Remove the gene ID column
tx2gene <- dplyr::select(tx2gene, -GENEID)

## Get the quant files and metadata ##
# Collect the sample quant files
samples <- list.dirs('Salmon.out/', recursive = FALSE, full.names = FALSE)
quant_files <- file.path('Salmon.out', samples, 'quant.sf')
names(quant_files) <- samples
print(quant_files)

# Ensure each file actually exists
file.exists(quant_files)  # all should be TRUE

# Set up metadata frame
colData <- data.frame(
  row.names = samples,
  condition = rep(c('untreated', 'dex'), 4)
)

## Compile the tximport counts object and make DESeq dataset ##

# Get tximport counts object
txi <- tximport(files = quant_files, 
                type = 'salmon',
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE)

# Make DESeq dataset
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = colData,
                                design = ~condition)

### Do DESeq analysis ! ###

# PCA
vsd <- vst(dds)
plotPCA(vsd)

# DEG analysis
dds <- DESeq(dds)

# Get the results
resdf <- results(dds)



############ Hands on Activity (time permititng) ############ 

# Get the raw data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103843

# Specifically, get SRR6035978 SRR6035979 SRR6035982 SRR6035983

# Then import to R using tximport 

# Run typical DEG analysis comparing siFLI1 to siNEG




