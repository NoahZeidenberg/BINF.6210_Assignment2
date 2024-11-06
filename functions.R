# get the species name, ID and annotation entry from a gene query list
extract_info <- function(entry) {
  # extract info using regex expressions
  gene_name <- str_extract(entry, "^[^\\n]+") # get all the characters before the first newline
  species <- str_extract(entry, "\\[([A-Za-z]+\\s[A-Za-z]+)") # extract based on the assumption species is in square brackets
  id <- str_extract(entry, "ID: \\d+") # regex expression to extract the numbers following "ID: "
  annotation <- str_extract(entry, "Annotation: .+") # regex expression to extract matching line
  
  # extract accession_id, bp_start, and bp_stop from the Annotation field
  accession_id <- str_extract(annotation, "[A-Z]+_\\d+\\.\\d")
  bp_start <- str_extract(annotation, "\\((\\d+)\\.\\.")
  bp_stop <- str_extract(annotation, "\\.\\.(\\d+)")
  
  # clean up extracted values
  species <- str_remove(species, "\\[") # Remove the leading square bracket
  id <- str_replace(id, "ID: ", "") # remove the leading "ID:" characters
  bp_start <- str_replace(bp_start, "\\(", "") %>% str_replace_all("\\.", "")  # remove bracket and periods
  bp_stop <- str_replace(bp_stop, "\\)", "") %>% str_replace_all("\\.", "") # ""
  
  return(data.frame(gene = gene_name, species_name = species, ID = id,
                    accession_id = accession_id,
                    bp_start = bp_start, bp_stop = bp_stop, 
                    stringsAsFactors = FALSE))
}

# remove header from fasta
clean_fasta <- function(sequence) { 
  
  # extract the sequence itself
  clean_seq <- sequence %>%
    strsplit(split = "\n") %>% # split by newline
    .[[1]] %>%
    .[-1] %>% # remove the first line (the header)
    paste(collapse = "") # collapse back into one string, without newline characters
  
  return(clean_seq)
}

one_hot_encode <- function(sequence) {
  dna_bases <- c("A", "C", "G", "T")
  seq_split <- strsplit(sequence, "")[[1]]
  one_hot_matrix <- sapply(dna_bases, function(base) as.integer(seq_split == base))
  t(one_hot_matrix)  # Transpose to get bases as columns
}

generate_kmer_counts <- function(sequence, k) {
  kmer_counts <- oligonucleotideFrequency(DNAString(sequence), width = k)
  as.numeric(kmer_counts)
}

# write fasta files to appropriate folder for later use with jellyfish and randomForest
write_fasta_files <- function(df_seqspec, output_dir = "./doc") {
  
  # helper function to write a single FASTA file
  write_fasta <- function(gene, accession_id, sequence) {
    
    # I already created a directory for each gene (BMP4 and AHSG)
    gene_dir <- file.path(output_dir, gene) # file.path automatically separates by '/'
    
    # define the file path
    fasta_file <- file.path(gene_dir, paste0(accession_id, ".fasta"))
    
    # write the sequence to a FASTA file
    cat(paste0(">", accession_id, "\n", sequence), file = fasta_file) # use cat's connection feature to print to the file
  }
  # apply write_fasta to each row of the dataframe using parallel (multi-argument) map
  purrr::pmap(df_seqspec, write_fasta)
  
  cat("FASTA files created.\n") # just to show success on output
}

jellyfish <- function(fasta_files) {
  
  # define a helper function to calculate k-mers for a single FASTA file
  calculate_kmers <- function(file_path) {
    
    # extract gene and accession id from the file path
    gene_dir <- dirname(file_path) # get parent folder (BMP4 or AHSG)
    accession_id <- tools::file_path_sans_ext(basename(file_path)) # remove file extension and path leading up to the file name
    
    # using k-mers of k-size 7, 14 and 21.
    kmer_sizes = c(7, 14, 21)
    
    # run Jellyfish count command + dump to text file
    purrr::walk(kmer_sizes, function(k) {
      
      # define paths for output files
      output_jf <- file.path(gene_dir, paste0(accession_id, "output_", k, "mer.jf"))
      output_txt <- file.path(gene_dir, paste0(accession_id, "kmer_counts_", k, "mer.txt"))
      
      # run Jellyfish count and dump commands through bash
      system2("jellyfish", args = c("count", "-m", k, "-s", "100M", "-o", output_jf, file_path))
      system2("jellyfish", args = c("dump", "-c", output_jf), stdout = output_txt)
      
      # No need for the .jf file
      file.remove(output_jf)
    })
  }
  
  # apply the function to each file in path
  walk(fasta_files, calculate_kmers)
  
  # output success to console
  cat("K-mer counting complete and results saved.\n")
}

load_jellyfish_kmers <- function(accession_id, gene, k) {
  # write path to variable using the same naming system as above
  kmer_file <- file.path("./doc", gene, paste0(accession_id, "kmer_counts_", k, "mer.txt"))
  
  # read counts into dataframe
  kmer_counts <- read.table(kmer_file, header = FALSE, col.names = c("kmer", "count"))
  return(kmer_counts %>%
           tidyr::pivot_wider(names_from = kmer, values_from = count, # transform s.t. columns are the k-mer name (e.g. 'AAAGTAA') and the row contains counts
                              values_fill = 0)) # any missing kmers are assigned a count of 0
}

# Combine all k-mers into a single dataframe ** could not compute with k > 7, my computer kept crashing ** 
combine_kmers <- function(df) {
  # load k-mers for the jellyfish files
  k_sizes <- c(7, 14, 21) 
  
  gene_data <- df %>%
    rowwise() %>% # using rowwise so that the x*y pairs are not automatically combined/flattened by mutate()
    mutate(kmers_7 = list(load_jellyfish_kmers(accession_id, gene, 7)), # generate kmers
           kmers_14 = list(load_jellyfish_kmers(accession_id, gene, 14)), # ""
           kmers_21 = list(load_jellyfish_kmers(accession_id, gene, 21))) %>% # ""
    unnest_wider(c(kmers_7, kmers_14, kmers_21)) %>% # expand kmers back to individual columns
    bind_cols(df_1mer, df_2mer, df_3mer, .)  # Add existing 1-mer, 2-mer, 3-mer counts
}

# load large-k k-mer counts from Jellyfish text files
kmer_counts_from_files <- function(gene_dir, k) {
  # list of files matching the predefined pattern
  files <- list.files(gene_dir, pattern = paste0("kmer_counts_", k, "mer.txt"), full.names = TRUE)
  
  # map over files and bind by row
  kmer_data <- purrr::map_dfr(files, function(file) {
    accession_id <- sub(paste0("kmer_counts_", k, "mer.txt"), "", basename(file))
    
    # read counts, convert to probabilities
    counts <- read.table(file, header = FALSE, col.names = c("kmer", "count"))
    counts[is.na(counts)] <- 0 # try removing NAs here instead of later
    total_count <- sum(counts$count)
    counts$probability <- counts$count / total_count
    
    # Pivot to make each k-mer its own column with probability values
    counts <- counts %>%
      select(-count) %>%  # Remove raw counts if only probabilities are needed
      pivot_wider(names_from = kmer, values_from = probability, values_fill = list(probability = 0))
    
    # Add accession_id and gene columns
    counts <- counts %>%
      mutate(accession_id = accession_id, gene = ifelse(grepl("BMP4", gene_dir), "BMP4", "AHSG")) %>%
      select(accession_id, gene, everything())  # Order columns
    
    return(counts)
  })
  
  return(kmer_data)
}

# train RF model
RF_confMatrix <- function(df_kmers, k) {
  data <- df_kmers %>% select(-accession_id) # remove accession_id
  data$gene <- as.factor(data$gene) # set gene as factor
  
  # partition the data into 70/30 split (train/test)
  training_index <- caret::createDataPartition(data$gene, p = 0.7, list = FALSE) # create indices for training
  train_data <- data[training_index, ] # 70%
  test_data <- data[-training_index, ] # 30%
  
  # Train the model given 500 trees (should be appropriate for k's 1-7)
  rf_model <- randomForest(gene ~ ., data = train_data, ntree = 500)
  predictions <- stats::predict(rf_model, test_data) # assess model performance
  
  # Output a confusion matrix 
  conf_matrix <- caret::confusionMatrix(predictions, test_data$gene)
  
  # Save the confusion matrix if running for the first time
  if (!file.exists(paste0("./doc/confusionMatrix_", k, "mer.txt"))) {
    readr::write_rds(conf_matrix, paste0("./doc/confusionMatrix_", k, "mer.txt"))
  }
  
  # show successful run
  cat(paste0("Successfully computed confusion matrix for rf model of ", k, "-mers\n"))
}

# load DNN history object into a dataframe, to build a line graph with shaded bands for train and test accuracy/loss values
load_history <- function(history_file) {
  history_df <- data.frame(readRDS(history_file), epoch = 1:30) %>%
    pivot_longer(cols = c(accuracy, val_accuracy, loss, val_loss), 
                 names_to = "metric",
                 values_to = "value") %>%
    mutate(type = ifelse(grepl("val", metric), "Validation", "Training"),
           metric = ifelse(grepl("accuracy", metric), "Accuracy", "Loss"))
  return(history_df)
}

# Generate ggplot
DNN_plot <- function(history_file, k) {
  
  # load history object in
  history_metrics <- load_history(history_file)
  
  # Separate accuracy and loss data
  accuracy_data <- history_metrics %>% filter(metric == "Accuracy")
  loss_data <- history_metrics %>% filter(metric == "Loss")
  
  # Calculate min and max for loss
  min_loss <- min(loss_data$value)
  max_loss <- max(loss_data$value)
  
  # Apply min-max scaling to loss and scale to 10%
  loss_data <- loss_data %>%
    mutate(scaled_loss = 0.1 * (value - min_loss) / (max_loss - min_loss))
  
  # Join scaled loss data with accuracy data for plotting
  plot_data <- accuracy_data %>%
    left_join(loss_data %>% select(epoch, type, scaled_loss), by = c("epoch", "type"))
  
  # Plot with ggplot2
  plot <- ggplot(plot_data, aes(x = epoch, y = value, color = type, fill = type)) +
    # Line plot for accuracy
    geom_line() +
    # Shaded error bands using scaled loss for margins around accuracy
    geom_ribbon(aes(ymin = value - scaled_loss, ymax = value + scaled_loss), alpha = 0.2) +
    labs(
      title = paste0("DNN performance for k-mer length ", k),
      subtitle = "Training and Testing Loss Values are Scaled to 10%",
      x = "Epoch",
      y = "Accuracy",
      color = "Dataset Type",
      fill = "Dataset Type"
    ) +
    theme_minimal() +
    theme(plot.subtitle=element_text(size=8, face="italic", color="grey20"),
          legend.position = "right",
          legend.key.size = unit(0.5, 'cm'))
  
  return(plot)
}
