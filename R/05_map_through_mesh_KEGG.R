# Load necessary libraries
library(RISmed)
library(tidyverse)
library(ReactomePA)
library(org.Hs.eg.db)
library(clusterProfiler)

# Read data
mesh <- XML::xmlToDataFrame("C:/Users/SaraAhsani-Nasab/OneDrive - Unit of Biostatistics Epidemiology and Public Health/PEA/MeSH/desc2023.xml")
terms <- mesh$DescriptorName
pubtator <- read.delim(("C:/Users/SaraAhsani-Nasab/OneDrive - Unit of Biostatistics Epidemiology and Public Health/PEA/Pubtator/gene2pubtatorcentral.gz"),
                       quote = "", header = TRUE,
                       col.names = c('PMID', 'Object', 'Gene', 'Gene_name', 
                                     'Dataset'))

# Seeds
seeds <- c(5, 10, 169, 3011, 19811)


# KEGG --------------------------------------------------------------------

# Function to query, get genes and pathways
query_to_pathwaysKEGG <- function(x) {
  somma <- NA
  matched_pathways <- character(0)
  
  tryCatch({
    # Query execution
    search_topic <- paste0('("', x, '"[Mesh]) AND ("2018/01/01"[Date - Entrez] : "2018/12/31"[Date - Entrez])')
    search_query <- EUtilsSummary(search_topic, retmax = 1000)
    
    # Gene retrieval
    filtered_data <- data.frame('PMID' = (search_query@PMID))
    common_als <- intersect(pubtator$PMID, as.numeric(filtered_data$PMID))
    annot_als <- pubtator[pubtator$PMID %in% common_als,]
    
    # Pathway analysis
    KEGG_res <- enrichKEGG(gene = annot_als$Gene, organism = "hsa",
                           keyType = "kegg", #keyType can be also "ncbi-geneid"
                           pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.01,
                           use_internal_data = FALSE)
    
    path <- dplyr::tibble('Pathway' = KEGG_res@result[["Description"]],
                          'q_value' = KEGG_res@result[["qvalue"]]) %>% 
      dplyr::filter(q_value < 0.01) %>% 
      rowid_to_column(var = "rowid")
    
    path$coviddi <- stringr::str_detect(path$Pathway, 
                                        stringr::regex("COVID|coronavirus|SARS-CoV-2|COVID-19", 
                                                       ignore_case = TRUE))
    
    somma <- sum(path$coviddi, na.rm = TRUE)
    matched_pathways <- path$Pathway[path$coviddi]
    pathway_sig <- path[path$Pathway %in% matched_pathways, ]
    
  }, error = function(e) {
    warning(paste("An error occurred for the term", x, ": ", e))
  })
  
  return(list("sum" = somma, "matched_pathways" = matched_pathways))
}

# Function to execute pipeline for a given seed
execute_pipelineKEGG <- function(seed) {
  set.seed(seed)
  terms_2 <- sample(terms, size = 2, replace = FALSE)
  
  results <- map(terms_2, query_to_pathwaysKEGG, .progress = "progress")
  
  somma_list <- map(results, 'sum')
  matched_pathways_list <- map(results, 'matched_pathways')
  
  pathways_df <- dplyr::tibble('term' = terms_2, 
                               'sum' = somma_list, 
                               'matched_pathways' = matched_pathways_list,
                               'seed' = seed)
  
  return(pathways_df)
}

# Execute pipeline for each seed and bind rows
final_results_dfKEGG <- map_dfr(seeds, execute_pipelineKEGG, .progress = "progress")

# If 'matched_pathways' is a list column in final_results_df
final_results_dfKEGG$matched_pathways <- map_chr(final_results_dfKEGG$matched_pathways, function(x) {
  paste(x, collapse = ",")
})

final_results_dfKEGG$sum <- map_chr(final_results_dfKEGG$sum, ~ {
  if (is.na(.x)) return("NA")
  return(as.character(.x))
})

write.csv(final_results_dfKEGG, "final_results_dfKEGG.csv")

save(final_results_dfKEGG, file = 'final_resutsKEGG.rda')


# Overlaps ----------------------------------------------------------------

## convert path$geneID to a numeric variable
split_strings <- strsplit(path$geneID, "/")
numeric_list <- lapply(split_strings, function(x) as.numeric(x))
path$geneID <- sapply(numeric_list, function(x) x)

## overlaps of genes among pathways (3 most significant/higher q-value)
row1_values <- path$geneID[[1]]
row2_values <- path$geneID[[2]]
row3_values <- path$geneID[[3]]

overlap_13 <- intersect(row1_values, row3_values)

if (length(overlap_13) > 0) {
  cat("Overlap between rows", 1, "and", 3, "occurs at values:", overlap_13, "\n")
} else {
  cat("There is no overlap between rows", 1, "and", 3, "\n")
}


overlap_13 <- intersect(row1_values, row3_values)

if (length(overlap_23) > 0) {
  cat("Overlap between rows", 2, "and", 3, "occurs at values:", overlap_23, "\n")
} else {
  cat("There is no overlap between rows", 2, "and", 3, "\n")
}