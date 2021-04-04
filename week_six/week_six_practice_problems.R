# 1. Find a data set on this page: http://www.bioconductor.org/packages/release/data/experiment/

# 2. Install and load the data set

# 3. Convert the data to a DESeq dataset object with an appropriate design formula

# 4. Run DESeq

# 5. Get the results 

# 6. If not already done, convert the gene symbols to gene IDs

# 7. Use EnhancedVolcano to plot the results

# 8. Save the results to a .csv file 

# When finished, send your R script to Henry Miller (millerh1@uthscsa.edu)
# Make sure to explain what comparison you made with DESeq and why.

# For example -- to test the effect of dex on gene expression in the "airway" dataset, I would start with...
library(airway)
dds <- DESeqDataSet(airway, design = ~dex)
# And then continue the analysis from there...




