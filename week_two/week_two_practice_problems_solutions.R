### Problem One ###
#Create the vector fruit

fruit <- c('Apple', 'Orange', 'Passion fruit', 'Banana')
# prt 1
# Create the for statement
for ( i in fruit) {
  print(i)
}

# prt 2
for(var in 1:length(fruit)) {
  print(fruit[var])
}

# prt 3
# skip over the value "Orange"
for(var in 1:length(fruit)) {
  if (fruit[var] ==  'Orange'){
    next
  }
  print(fruit[var])
}

# prt 4
#Complete the following function that take in input X and outputs the cube of X which is equal to Y
Cube_It <- function(x){
  y <- x^3
  print(paste0("The Cube of ", x, " is ", y ))
}

# now use function on the nubmer 42
Cube_It(42)


### Problem Two ###

# You have a class with several students -- They just took an exam!
# Here are the results (run the code):

students <- c("Johnny", "Frank", "Timmy", "Alice", "Susan")
num_grades <- c(88, 90, 76, 98, 85)
names(num_grades) <- students 
num_grades

# Problem: Convert the numeric grades (e.g., 85) to letter grades (e.g., "B")
# Assume standard conversions... 90+ is "A", 80-89 is "B", etc
# Loop through the vector of numeric grades. Convert to letter grades.
# Save the results in a new vector, "letter_grades"
# Finally, make a dataframe, "student_grades", with columns made of "letter_grades", "num_grades", and "students"
# Save "student_grades" as a csv file.


# Solution #
letter_grades <- c()
for (i in 1:length(num_grades)) {
  grade <- num_grades[i]
  
  if (grade >= 90) {
    letter_grades[i] <- "A"
  } else if (grade >= 80) {
    letter_grades[i] <- "B"
  } else if (grade >= 70) {
    letter_grades[i] <- "C"
  } else if (grade >= 60) {
    letter_grades[i] <- "D"
  } else {
    letter_grades[i] <- "F"
  }
}
student_grades <- data.frame(letter_grades, num_grades, students)
write.csv(student_grades, file = "student_grades.csv")
############



### Problem Two ###



.....



### Challenge Problem ###

# You just discovered some genes which are associated with cancer progression. Cool!
# However... you only know the Ensembl gene IDs (e.g., ENSG00000141510) of these genes...
# These IDs are hard to read, so you need to convert them to gene symbols (e.g., TP53)

cancer_genes <- c("ENSG00000139618", "ENSG00000106462", "ENSG00000116288")
cancer_genes

# Problem: Make a function to convert Ensembl gene IDs to gene symbols.     
# You have also been provided with a data.frame has the mapping of ids to symbols

id2symbol <- data.frame(
  "Ensembl" = c("ENSG00000141510", "ENSG00000139618", "ENSG00000106462", "ENSG00000116288"),
  "gene_symbol" = c("TP53", "BRCA2", "EZH2", "PARK7")
)
id2symbol

# Requirements:                                                                  
#   1. The function, "gene_id_converter", takes one argument, "gene_id"                     
#   2. The function must return the corresponding gene symbol      

# Solution #
gene_id_converter <- function(gene_id) {
  
  # Get mapping between IDs and symbols
  id2symbol <- data.frame(
    "Ensembl" = c("ENSG00000141510", "ENSG00000139618", "ENSG00000106462", "ENSG00000116288"),
    "gene_symbol" = c("TP53", "BRCA2", "EZH2", "PARK7")
  )
  
  # Conversion code
  gene_symbol <- id2symbol$gene_symbol[id2symbol$Ensembl == gene_id]
  
  # Return symbol
  return(gene_symbol)
}
############

# Test -- this should output "TP53"
gene_id_converter(gene_id="ENSG00000141510")


### SUPER Challenge Problem ###

# Your gene_id_converter is working splendidly! But now you've got a new problem:
# Your colleague just sent you a ton of new genes to convert! None of these new
# genes are in your id2symbol data frame... 

# Problem: Extend you gene_id_converter function to convert ANY ensembl gene id to a gene symbol
# To assist you, your advisor sent you a csv file that contains the mapping of all ids to symbols. 

# This the mapping file and it's on Box
"gene_id_to_symbol.csv" 

# Solution #
gene_id_converter <- function(gene_id) {
  
  # Get mapping between IDs and symbols
  id2symbol <- read.csv(file = "gene_id_to_symbol.csv")
  
  # Conversion code
  result <- id2symbol$gene_symbol[which(id2symbol$Ensembl == gene_id)]
  
  # Return result
  return(result)
}
############

# Test -- this should output "BRCA1"
gene_id_converter(gene_id="ENSG00000012048")


### ULTRA Challenge Problem ###

# Good work! Your converter is working perfectly! The only problem:
# You just got a new list of genes IDs from a collaborator.. and they're in Entrez format, not Ensembl!
# Fortunately, the csv file your advisor gave you also contains Entrez IDs!

entrez_genes_from_collaborator <- c(8243, 23657, 472)
entrez_genes_from_collaborator

# Problem: Extend you gene_id_converter function to convert Entrez gene ids to gene symbols.
# Requirement: Add a new argument, "type", which is used to specify whether the gene ID is "Entrez" or "Ensembl" 

# Solution #
gene_id_converter <- function(gene_id, type) {
  
  # Get mapping between IDs and symbols
  id2symbol <- read.csv(file = "gene_id_to_symbol.csv")
  
  # Conversion code 
  if (type == "Entrez") {
    result <- id2symbol$gene_symbol[which(id2symbol$Entrez == gene_id)]
  } else if (type == "Ensembl") {
    result <- id2symbol$gene_symbol[which(id2symbol$Ensembl == gene_id)]
  }
  
  # Return result
  return(result)
}
############

# Test 1 -- this should output "BRCA1"
gene_id_converter(gene_id="ENSG00000012048", type = "Ensembl")

# Test 2 -- this should output "ATM"
gene_id_converter(gene_id=472, type = "Entrez")


### EXTREME Challenge Problem ###

# Perfecto! Only one problem: now we need to convert symbols to ids!

# Problem: extend the gene_id_converter function so it can convert symbols to gene IDs!
# Requirement: instead of "gene_id", it should now accept the argument "gene" as this 
# argument can be either a gene id or a symbol. 
# The "type" argument should now accept the option "gene_symbol"

# Solution #
gene_id_converter <- function(gene, type) {
  
  # Get mapping between IDs and symbols
  id2symbol <- read.csv(file = "gene_id_to_symbol.csv")
  
  # Conversion code 
  if (type == "Entrez") {
    result <- id2symbol$gene_symbol[which(id2symbol$Entrez == gene)]
  } else if (type == "Ensembl") {
    result <- id2symbol$gene_symbol[which(id2symbol$Ensembl == gene)]
  } else if (type == "gene_symbol") {
    result <- id2symbol[which(id2symbol$gene_symbol == gene), c(1, 2)]
  }
  
  # Return result
  return(result)
}
############

# Test 1 -- this should output "ENSG00000012048" and "672"
gene_id_converter(gene = "BRCA1", type = "gene_symbol")

# Test 2 -- this should output "BRCA1"
gene_id_converter(gene="ENSG00000012048", type = "Ensembl")

# Test 3 -- this should output "ATM"
gene_id_converter(gene=472, type = "Entrez")


### IMPOSSIBLE Challenge Problem ###

# Problem: extend the gene_id_converter function that it only needs one argument: "gene"

# Solution #
gene_id_converter <- function(gene) {
  
  # Get mapping between IDs and symbols
  id2symbol <- read.csv(file = "gene_id_to_symbol.csv")
  
  # Determine the type
  if (gene %in% id2symbol$Entrez) {
    type <- "Entrez"
  } else if (gene %in% id2symbol$Ensembl) {
    type <- "Ensembl"
  } else if (gene %in% id2symbol$gene_symbol) {
    type <- "gene_symbol"
  }
  
  # Conversion code 
  if (type == "Entrez") {
    result <- id2symbol$gene_symbol[which(id2symbol$Entrez == gene)]
  } else if (type == "Ensembl") {
    result <- id2symbol$gene_symbol[which(id2symbol$Ensembl == gene)]
  } else if (type == "gene_symbol") {
    result <- id2symbol[which(id2symbol$gene_symbol == gene), c(1, 2)]
  }
  
  # Return result
  return(result)
}
############


# Test 1 -- this should output "CDKN2A"
gene_id_converter(gene = "ENSG00000147889")

# Test 2 -- this should output "CDKN2A"
gene_id_converter(gene = 1029)

# Test 3 -- this should output ENSG00000147889 and 1029
gene_id_converter(gene = "CDKN2A")



### GIVE UP ###

# Problem: Extend the function to convert individual genes and a vector of genes!
# Requirement: to reflect this change, the function argument should now be "genes" instead of "gene"
# Requirement: the function must return results as a list with the names as the original "genes"

# Solution #
gene_id_converter <- function(genes) {
  
  # Get mapping between IDs and symbols
  id2symbol <- read.csv(file = "gene_id_to_symbol.csv")
  
  # Initialize the empty list
  result_list <- list()
  
  # Use a for-loop to convert the genes
  for (i in 1:length(genes)) {
    
    # Get the current gene
    gene <- genes[i]
    
    # Determine the type
    if (gene %in% id2symbol$Entrez) {
      type <- "Entrez"
    } else if (gene %in% id2symbol$Ensembl) {
      type <- "Ensembl"
    } else if (gene %in% id2symbol$gene_symbol) {
      type <- "gene_symbol"
    }
    
    # Conversion code 
    if (type == "Entrez") {
      result <- id2symbol$gene_symbol[which(id2symbol$Entrez == gene)]
    } else if (type == "Ensembl") {
      result <- id2symbol$gene_symbol[which(id2symbol$Ensembl == gene)]
    } else if (type == "gene_symbol") {
      result <- id2symbol[which(id2symbol$gene_symbol == gene), c(1, 2)]
    }
    
    # Add result to result_list
    result_list[[i]] <- result
    
  }
  
  # Add back in the original genes as the names
  names(result_list) <- genes
  
  # Return result
  return(result_list)
}
############

# Test -- this should output a list:
#
# $ENSG00000147889
# [1] "CDKN2A"
# 
# $`8243`
# [1] "SMC1A"
#
# $TP53
# Entrez         Ensembl
# 5973   7157 ENSG00000141510
#
gene_id_converter(genes = c("ENSG00000147889", 8243, "TP53"))


### GIVE UP challenge problem ###

# Problem: Same thing. But use lapply.

# Solution #
gene_id_converter <- function(genes) {
  
  # Get mapping between IDs and symbols
  id2symbol <- read.csv(file = "gene_id_to_symbol.csv")
  
  # Use lapply instead of a for-loop
  result_list <- lapply(genes, function(gene) {
    
    # Determine the type
    if (gene %in% id2symbol$Entrez) {
      type <- "Entrez"
    } else if (gene %in% id2symbol$Ensembl) {
      type <- "Ensembl"
    } else if (gene %in% id2symbol$gene_symbol) {
      type <- "gene_symbol"
    }
    
    # Conversion code 
    if (type == "Entrez") {
      result <- id2symbol$gene_symbol[which(id2symbol$Entrez == gene)]
    } else if (type == "Ensembl") {
      result <- id2symbol$gene_symbol[which(id2symbol$Ensembl == gene)]
    } else if (type == "gene_symbol") {
      result <- id2symbol[which(id2symbol$gene_symbol == gene), c(1, 2)]
    }
    
  })
  
  # Add back in the original genes as the names
  names(result_list) <- genes
  
  # Return result
  return(result_list)
}
############

# Test -- this should output a list:
#
# $ENSG00000147889
# [1] "CDKN2A"
# 
# $`8243`
# [1] "SMC1A"
#
# $TP53
# Entrez         Ensembl
# 5973   7157 ENSG00000141510
#
gene_id_converter(genes = c("ENSG00000147889", 8243, "TP53"))

