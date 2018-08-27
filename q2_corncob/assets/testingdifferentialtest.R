the_otu_table <- read.table(file = "/Users/paulinetrinh/Downloads/corncobdata/6ee53529-3073-4bd6-8f0a-0b89025de08a/data/otu_table.txt",
                            skip = 0, 
                            header = F, 
                            row.names = 1)
dim(the_otu_table)
colnames(the_otu_table) <- colnames(read.csv(file = "/Users/paulinetrinh/Downloads/corncobdata/6ee53529-3073-4bd6-8f0a-0b89025de08a/data/otu_table.txt", nrows=1, skip=1, sep = "\t"))[-1]


the_metadata <- read.csv(file = "/Users/paulinetrinh/Downloads/corncobdata/metadata.tsv",
                         header = T, sep = "\t", row.names = 1)


tax <- read.csv(file = "/Users/paulinetrinh/Downloads/corncobdata/taxonomy.tsv", 
                header = T, sep = "\t", row.names = 1)

the_taxonomy <- separate(data = tax, col = "Taxon", 
                         into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species" ), 
                         sep = ";")
the_taxonomy <- as.matrix(the_taxonomy)

dim(the_taxonomy)
ps <- phyloseq(otu_table(the_otu_table, taxa_are_rows = TRUE), 
               sample_data(the_metadata), 
               tax_table(the_taxonomy)
)

ps2 <- ps %>% tax_glom("Family")

differential <- function(variable) {
          differentialTest(formula = ~ variable,
                                  phi.formula = ~ variable,
                                  formula_null = ~ 1,
                                  phi.formula_null = ~ 1,
                                  data = ps,
                                  fdr_cutoff = 0.05,
                                  inits = rbind(rep(.01, 4)))
}

fullAnalysis <- differentialTest(formula = ~ ReportedAntibioticUsage,
                                 phi.formula = ~ ReportedAntibioticUsage,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,
                                 data = ps,
                                 fdr_cutoff = 0.05,
                                 inits = rbind(rep(.01, 4)))
tax_table(ps)

# write out table with fdr p-values 
fullAnalysisTable <- fullAnalysis[[2]]
fulltable <- as.data.frame(fullAnalysisTable)
fulltable_tax <- merge(fulltable,the_taxonomy, by = 0)
head(fulltable_tax)
dim(fulltable_tax)

colnames(fulltable_tax)[1] <- "Taxon" 

write.table(fulltable_tax, "teststuff.tsv" , sep = "\t", 
            row.names = F,
            quote = F)
