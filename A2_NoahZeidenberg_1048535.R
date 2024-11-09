## ----error = FALSE, message=FALSE, warning=FALSE----------------------------
# some global settings
knitr::opts_chunk$set(echo = TRUE) # always include code in output doc
options(install.packages.ask = FALSE) # no need to ask to install packages/dependancies
options(crayon.enabled = FALSE) # could not knit with invalid characters from caret

library("pacman") # pacman automatically checks if a package is installed and (if not, installs it then) loads it into the current environment
pacman::p_load("dplyr", # for performing major dataframe selections/transformations/etc.
               "tidyverse", # has numerous helpful functions I may use
               "readr", # ""
               "rentrez", # for using EUtils, connecting to NCBIs dbs
               "DECIPHER", # for performing and  multiple sequence alignments
               "taxize", # works well with the bold package -> e.g. "downstream" fn reduces load on BOLD query
               "xml2", # for parsing certain entrez responses
               "randomForest", # using an RF model for comparison
               "caret", # for confusion matrices 
               "keras", "tensorflow", # for training the DNN
               "purrr", # I prefer map functions to the base::apply family
               "ggplot2", "patchwork", # better plotting visualization and control than base plots
               update = FALSE,
               dependencies = TRUE) 



## ---------------------------------------------------------------------------
# set rentrez API key from environment variable
tryCatch(rentrez::set_entrez_key(Sys.getenv("ENTREZ_KEY")))

# source functions
source("./functions.R")

# set the seed
set.seed(36) # for reproducability. 36 is 18 (lucky number) * 2 (for assignment #2 :)


## ---------------------------------------------------------------------------
# load from text file rather than re-run rentrez functions
if (!exists("gene_data") & file.exists("./doc/gene_data.txt")) {
  gene_data <- read_file(file = "./doc/gene_data.txt")
}

# again, only run if necessary
if (!exists("gene_data")) { 
  
  # search from NCBI's gene database for either gene in Aves returns ~ 300 entries
  genes_query <- rentrez::entrez_search(db = "gene", 
                                        term = '(BMP4 OR AHSG OR "Alpha 2-HS Glycoprotein" OR "bone morphogenetic protein 4") AND "Aves"[porgn] AND alive[prop]',
                                        use_history = TRUE) # using web_history object to mitigate size-issues in the subsequent API request
  
  # pull info for each entry using web_history object
  gene_data <- rentrez::entrez_fetch(db = "gene", web_history = genes_query$web_history, rettype = NULL, retmode = "default")
}

# save as text file rather than re-run rentrez functions
if (!file.exists("./doc/gene_data.txt")) {
  write(gene_data, file = "./doc/gene_data.txt")
}


## ----echo = FALSE, error = FALSE, eval = FALSE, include = FALSE-------------
## # not needed in output ^^. just to show that yes, I did look at my data before working with it
## gene_data %>% cat # entries are formatted consistently, can apply a vectorized function...


## ---------------------------------------------------------------------------
# split up each entry for parsing
entries <- str_split(gene_data, "\\n\\d+\\.\\s+")[[1]] %>% # split at each integer followed by a '.' and at least one whitespace character
  trimws # trim white space around each line

entries %>%
    strsplit(split = "\n") %>% # split into substrings by newline
    .[[20]] %>% # take the 20th entry, for example
    cat(sep = "\n") # concatenate the entry by newline

# Apply vectorized function to extract the ID, gene name, species name, accession ID and the start/stop bp numbers identifying where the gene resides within the fasta file for said accession ID.
df_spec <- purrr::map_dfr(entries, extract_info) # equivalent (though faster) to base::apply but output is a dataframe

# keep only the rows where the gene name is BMP4 or AHSG (some entries are from rentrez returning similar genes such as BMP7 or PAX6, which are also involved in bone development)
df_spec <- df_spec %>%
                filter(gene %in% c("BMP4", "AHSG")) %>% # still ~ 270 observations
                drop_na # remove rows with NAs using tidyr's drop_na fn

# Visualize
df_spec %>% head(3)


## ----warning = FALSE--------------------------------------------------------
# load from file if already run previously
if (file.exists("./doc/df_seqspec.tsv")) {
  df_seqspec <- read_tsv("./doc/df_seqspec.tsv", show_col_types = FALSE)
}

# "" otherwise:
if (!exists("df_seqspec")) {
  # apply vectorized function that leverages purrr's pmap ("parallel" map) to iterate over each row in df_spec and perform the entrez_fetch() using multiple column values in df_spec.
  ls_fastas <- pmap(df_spec[,c("accession_id", "bp_start", "bp_stop")],
                    function(accession_id, bp_start, bp_stop) {
                      return(entrez_fetch(db = "nuccore", rettype = "fasta", 
                                   id = accession_id, seq_start = bp_start, seq_stop = bp_stop))
                      Sys.sleep(0.1) # to avoid rate-limiting
  })
  
  # export fastas to file, for use with Jellyfish later 
  ls_fastas %>%
    unlist %>% # flatten to 
    writeLines(., con = "./doc/BMP4_AHSG_Aves.fasta")
  
  # clean fasta files (remove header)
  clean_fastas <- map(ls_fastas, clean_fasta)

  # convert to dataframe and bind to df_spec
  df_seqspec <- do.call(rbind, clean_fastas) %>% # bind rows
                as.data.frame %>% # convert to df
                setNames("sequence") %>% # set the column name to "sequence"
                cbind(df_spec, .) # bind to df_spec
  
  # save as file
  write_tsv(df_seqspec, file = "./doc/df_seqspec.tsv")
}

# write fasta files for each gene to two separate folders, to be used in counting k-mers and training classifiers
if (!file.exists("./data/AHSG/NC_021681.1.fasta")) { # test case, it will exist only if this code is run
  write_fasta_files(df_seqspec[,c("gene", "accession_id", "sequence")])
}
# drop unnecessary rows
df_seqspec <- df_seqspec[, c("gene", "species_name", "accession_id", "sequence"), drop = TRUE]

# show quartiles for sequence length in AHSG and BMP4, respectively
cat("Summary of quartiles in sequence lengths for AHSG \n")
df_seqspec %>%
  filter(gene == "AHSG") %>%
  pull(sequence) %>%
  width %>%
  summary
cat("\n Summary of quartiles in sequence lengths for BMP4 \n")
df_seqspec %>%
  filter(gene == "BMP4") %>%
  pull(sequence) %>%
  width %>%
  summary


## ---------------------------------------------------------------------------
# computing 1-mer, 2-mer and 3-mers with Biostrings

# 1-mers as proportions
df_1mer <- df_seqspec$sequence %>%
              Biostrings::DNAStringSet(.) %>%
              Biostrings::letterFrequency(letters = c("A", "T", "G", "C"), as.prob = TRUE) %>%
              as.data.frame %>%
              cbind(df_seqspec[, c("gene", "accession_id"), drop = T], .)

# 2-mers as proportions
df_2mer <- df_seqspec$sequence %>%
              Biostrings::DNAStringSet(.) %>%
              Biostrings::dinucleotideFrequency(., as.prob = TRUE) %>%
              as.data.frame %>%
              mutate_all(~ round(., 2)) %>% # round to two decimal places
              cbind(df_seqspec[, c("gene", "accession_id"), drop = T], .)

# 3-mers as proportions
df_3mer <- df_seqspec$sequence %>%
              Biostrings::DNAStringSet(.) %>%
              Biostrings::trinucleotideFrequency(., as.prob = TRUE) %>%
              as.data.frame %>%
              mutate_all(~ round(., 2)) %>% # round to two decimal places
              cbind(df_seqspec[, c("gene", "accession_id"), drop = T], .)


## ---------------------------------------------------------------------------
# add a column to each data frame to indicate k-mer size
df_1mer_hist <- df_1mer %>% mutate(kmer_size = "1-mer")
df_2mer_hist <- df_2mer %>% mutate(kmer_size = "2-mer")
df_3mer_hist <- df_3mer %>% mutate(kmer_size = "3-mer")

# combine data frames
df_combined <- bind_rows(df_1mer_hist, df_2mer_hist, df_3mer_hist)

# reshape to long format (each k-mer probability in one row)
df_long <- df_combined %>%
  pivot_longer(cols = -c(accession_id, gene, kmer_size), names_to = "kmer", values_to = "probability")

# facet by both gene and k-mer size for a 2x3 plot 
hist <- ggplot(df_long, aes(x = probability)) +
  geom_histogram(binwidth = 0.005, fill = "skyblue", color = "grey30") +
  facet_wrap(~ gene + kmer_size, scales = "free_y") +  # Independent y-axis for each combination of gene and k-mer size
  labs(x = "Probability of k-mer occurrence", y = "Number of k-mers",
       title = "Histogram of k-mer probabilities for differing k-mer sizes in AHSG and BMP4") +
  theme_minimal()

# Save the plot if it doesn't already exist
if (!file.exists("./figs/small_mer_histogram.png")) {
  ggsave(filename = "./figs/small_mer_histogram.png", plot = hist, width = 12, height = 8, dpi = 600)
}
# Show plot in R markdown (below)


## ----eval = FALSE-----------------------------------------------------------
## # I attempted to compute larger k-mers using Jellyfish (installed in the bash shell according to the instructions on gmarcais/jellyfish github)
## 
## # I commented out certain lines since Jellyfish is likely only installed on my computer
## 
## #system("jellyfish --version") # check that jellyfish is working in the current R proj. directory
## 
## # computing 7-, 14- and 21-mers
## #file_paths <- c(list.files("./doc/BMP4", full.names = T),
## #                list.files("./doc/AHSG", full.names = T))
## # jellyfish(file_paths) # takes a very long time to run
## 
## # load from tsv files
## if (file.exists("./doc/df_7mer.tsv")) {
##   df_7mer <- read_tsv("./doc/df_7mer.tsv", show_col_types = FALSE)
##   df_7mer[is.na(df_7mer)] <- 0 # replace NAs with 0
## } else {
##   df_7mer <- bind_rows( # load data for each k-mer size
##     kmer_counts_from_files("./doc/BMP4", 7),
##     kmer_counts_from_files("./doc/AHSG", 7)
##   )
##   write_tsv(df_7mer, "./doc/df_7mer.tsv")
##   rm(df_7mer) # must be removed from local memory to proceed, otherwise my computer crashes
## }
## # again for 14-mers. My computer was not able to run this, I had to run it on another computer that has 16GB RAM
## #if (file.exists("./doc/df_14mer.tsv")) {
## #  df_14mer <- read_tsv("./doc/df_14mer.tsv", show_col_types = FALSE)
## #  df_14mer[is.na(df_14mer)] <- 0 # ""
## #} else {# Load data for each k-mer size
## #  df_14mer <- bind_rows(
## #    kmer_counts_from_files("./doc/BMP4", 14),
## #    kmer_counts_from_files("./doc/AHSG", 14)
## #  )
## #  write_tsv(df_14mer, "./doc/df_14mer.tsv")
## #  rm(df_14mer) # must be removed from local memory to proceed, otherwise my computer crashes
## #}
## 
## # My computer was not able to compute 21-mers, neither was a 16GB RAM computer.
## #if (file.exists("./doc/df_21mer.tsv")) {
## #  df_21mer <- read_tsv("./doc/df_21mer.tsv", show_col_types = FALSE)
## #  df_21mer[is.na(df_21mer)] <- 0 # ""
## #} else {# Load data for each k-mer size * not run
## #  df_21mer <- bind_rows(
## #    kmer_counts_from_files("./doc/BMP4", 21),
## #    kmer_counts_from_files("./doc/AHSG", 21)
## #  )
## #  write_tsv(df_21mer, "./doc/df_21mer.tsv")
## #  rm(df_21mer) # must be removed from local memory to proceed, otherwise my computer crashes
## #}


## ---------------------------------------------------------------------------
# apply defined function to each set of kmer data
if(!file.exists("./doc/confusionMatrix_7mer.txt")) {
  map2(list(df_1mer, df_2mer, df_3mer, df_7mer), 
     c("1", "2", "3", "7"),
    RF_confMatrix)
}

# an example of the computed confusion matrices
cat("An example of the computed confusion matrices, for 1-mer data\n")
readRDS("./doc/confusionMatrix_1mer.txt")


## ----eval = FALSE, message=FALSE, warning=FALSE, sanitize=TRUE--------------
## # I included only the example for 1mers, as I struggled to vectorize the code and I would have been over the page limit. However, the only change I made was to "k".
## 
## # one-hot encode sequence data
## one_hot_encoded_data <- lapply(df_seqspec$sequence, one_hot_encode)
## 
## # first bind the feature matrix and labels
## feature_matrix <- t(sapply(df_seqspec$sequence, generate_kmer_counts, k = 1))
## training_data <- data.frame(feature_matrix, label = as.factor(df_seqspec$gene))
## 
## # split  data into training and test (70/30 split)
## trainIndex <- createDataPartition(training_data$label, p = 0.7, list = FALSE)
## testIndex <- setdiff(seq_len(nrow(training_data)), trainIndex)
## train_data <- training_data[trainIndex, ]
## test_data <- training_data[-trainIndex, ]
## 
## # need to define model input shape for tensorflow
## sequence_length <- ncol(one_hot_encoded_data[[1]])
## input_shape <- c(sequence_length, 4) # 4 for A, T, G and C
## 
## # define the layers and parameters in the model, according to reference [3]
## model <- keras_model_sequential()
## # add layers individually
## model$add(layer_conv_1d(filters = 64, kernel_size = 16, activation = "relu", input_shape = input_shape))
## model$add(layer_max_pooling_1d(pool_size = 13, strides = 13))
## model$add(bidirectional(layer = layer_lstm(units = 64, return_sequences = TRUE)))
## model$add(layer_dropout(rate = 0.5))
## model$add(layer_flatten())
## model$add(layer_dense(units = 512, activation = "relu"))
## model$add(layer_dropout(rate = 0.5))
## model$add(layer_dense(units = 1, activation = "sigmoid"))
## 
## # compile using TensorFlow
## tf$keras$models$Sequential$compile(
##   model,
##   loss = "binary_crossentropy",
##   optimizer = "adam",
##   metrics = list("accuracy")
## )
## 
## # convert list of matrices to array
## train_x <- array(unlist(one_hot_encoded_data[trainIndex]), dim = c(length(trainIndex), sequence_length, 4))
## test_x <- array(unlist(one_hot_encoded_data[-trainIndex]), dim = c(length(testIndex), sequence_length, 4))
## # add labels
## train_y <- as.numeric(train_data$label == "AHSG")  # Adjust labels as needed
## test_y <- as.numeric(test_data$label == "BMP4")
## 
## # convert train and test data to TensorFlow-compatible format
## train_x <- array_reshape(train_x, dim = c(dim(train_x)[1], dim(train_x)[2], dim(train_x)[3]))
## test_x <- array_reshape(test_x, dim = c(dim(test_x)[1], dim(test_x)[2], dim(test_x)[3]))
## train_y <- as.array(train_y)
## test_y <- as.array(test_y)
## 
## # I could only fit using TensorFlow's backend format
## history_1mer <- tf$keras$models$Sequential$fit(
##   model,
##   x = train_x,
##   y = train_y,
##   epochs = as.integer(30),
##   batch_size = as.integer(32),
##   validation_data = list(test_x, test_y)
## )
## 
## # save history object
## history_metrics <- list(
##   accuracy = history_1mer$history$accuracy,
##   val_accuracy = history_1mer$history$val_accuracy,
##   loss = history_1mer$history$loss,
##   val_loss = history_1mer$history$val_loss
## )
## saveRDS(history_metrics, file = "./doc/DNN_history_1mer.rds")# Create figure comparing accuracy per epoch for all the trained models
## 


## ---------------------------------------------------------------------------
# Create plot comparing accuracy, loss, val_accuracy and val_loss for the four (1, 2, 3, 7) models
p1 <- DNN_plot("./doc/DNN_history_1mer.rds", k = 1)
p2 <- DNN_plot("./doc/DNN_history_2mer.rds", k = 2)
p3 <- DNN_plot("./doc/DNN_history_3mer.rds", k = 3)
p4 <- DNN_plot("./doc/DNN_history_7mer.rds", k = 7)

# Output as one figure
combined_plot <- (p1 + p2 + p3 + p4) +
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = 'A')  # Labels each plot as A, B, C, D

# Save plot
if (!file.exists("./figs/DNN_performance.png")) {
  ggsave(filename = "./figs/DNN_performance.png", plot = combined_plot, width = 12, height = 8, dpi = 600)
}

