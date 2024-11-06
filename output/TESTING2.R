# Attempt 2, using SNP database instead

pacman::p_load("tidyverse",
               "rentrez",
               "dplyr",
               "purrr",
               "progress",
               "msa",
               "DECIPHER",
               "Biostrings",
               "xml2", # for working with xml files from entrez
               "ggplot2", update = FALSE)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Get data ----
search_query <- '"Homo sapiens"[Organism] AND MT[Chromosome]'

# Query SNP database
SNP_search <- entrez_search(db = "snp",
                            term = search_query,
                            retmode = "xml",
                            retmax = 50,
                            use_history = TRUE)

# See links to clinvar, then links to pubmed
SNP_to_clinvar <- entrez_link(dbfrom = "snp", 
                              web_history = SNP_search$web_history,
                              db = "clinvar")

# SNP_to_pubmed
SNP_to_pubmed <- entrez_link(dbfrom = "snp",
                             web_history = SNP_search$web_history,
                             db = "pubmed")

# Search pubmed
test <- list()
test <- entrez_fetch(db = "pubmed",
                     id = SNP_to_pubmed$links,
                     web_history = TRUE,
                     rettype = NULL)
