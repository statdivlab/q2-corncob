#!/usr/bin/env Rscript

cat(R.version$version.string, "\n")
args <- commandArgs(TRUE)

otu.file <- args[[1]] # FeatureTable[Frequency]
meta.file <- args[[2]] # Metadata file 
tax.file <- args[[3]] # Taxonomy file
variable <- args[[4]] #some kind of string
out.file <- args[[5]] # TSV file output? 


###############################
# PACKAGE INSTALLATION CHECKS #
###############################

# Phyloseq install check
if ("phyloseq" %in% installed.packages()[,"Package"]) {
  cat("Great! phyloseq is already available\n\n")
} else {
  cat("Error: Phyloseq needs to be installed!\n\n")
}
library(phyloseq)

# Devtools install check
if ("devtools" %in% installed.packages()[,"Package"]) {
  cat("Thank goodness! devtools is already available\n\n")
} else {
  cat("Error: Devtools needs to be installed!\n\n")
}

# Magrittr install check 
if ("magrittr" %in% installed.packages()[,"Package"]) {
  cat("You have magrittr installed!\n\n")
} else {
  cat("Error: Magrittr needs to be installed!\n\n")
}
library(magrittr)

# tidyr install check 
if ("tidyr" %in% installed.packages()[,"Package"]) {
  cat("You have tidyr installed!\n\n")
} else {
  cat("Error: tidyr needs to be installed!\n\n")
  
}
devtools::install_github("hadley/tidyr")
library(tidyr)


errQuit <- function(mesg, status=1) {
  message("Error: ", mesg)
  q(status=status)
}

################################
# VALIDATE & READ IN OTU TABLE #
################################

if(!file.exists(otu.file)) {
  errQuit("Input file does not exist.")
} else {
  the_otu_table <- read.table(file = otu.file,
                              skip = 0, 
                              header = F, 
                              row.names = 1
  )
}
colnames(the_otu_table) <- colnames(read.csv(otu.file, nrows=1, skip=1, sep = "\t"))[-1]



# 776 taxa and 34 samples  
#dim(the_otu_table)
################################
# VALIDATE & READ IN METADATA  #
################################

the_metadata <- read.csv(file = meta.file,
                         header = T, sep = "\t", row.names = 1)

#dim(the_metadata)
################################
# VALIDATE & READ IN TAXONOMY  #
################################
tax <- read.csv(file = tax.file, 
                         header = T, sep = "\t", row.names = 1)


library(tidyr)
the_taxonomy <- separate(data = tax, col = "Taxon", 
         into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species" ), 
         sep = ";")
the_taxonomy <- as.matrix(the_taxonomy)

#dim(the_taxonomy)

# 776 by 8 
##########################
# CREATE PHYLOSEQ OBJECT #
##########################

ps <- phyloseq(otu_table(the_otu_table, taxa_are_rows = TRUE), 
               sample_data(the_metadata), 
               tax_table(the_taxonomy)
)

print(ps)
####NOTE: REMOVE THIS OUT OF R-SCRIPT AND INTO README.md 

# Install corncob using:
devtools::install_github("bryandmartin/CORNCOB")
# To begin, we load our example data set as a phyloseq object.
library(corncob)
sample_data(the_metadata)

##########################
# TIME TO USE CORNCOB!!! #
##########################
print(variable)
set.seed(1)
fullAnalysis <- differentialTest(formula = ~ ReportedAntibioticUsage,
                                 phi.formula = ~ ReportedAntibioticUsage,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,
                                 data = ps,
                                 fdr_cutoff = 0.05,
                                 inits = rbind(rep(.01, 4)))

# write out table with fdr p-values 
fullAnalysisTable <- fullAnalysis$p_fdr
fulltable <- as.data.frame(fullAnalysisTable)
fulltable_tax <- merge(fulltable,the_taxonomy, by = 0)

colnames(fulltable_tax)[1] <- "Taxon" 

fulltable_tax <- fulltable_tax[c("Taxon","DA", "DV", "Kingdom", "Phylum","Class","Order", "Family","Genus","Species","Confidence")]

write.table(fulltable_tax, out.file , sep = "\t", 
            row.names = F, 
            quote = F)

q(status=0)
