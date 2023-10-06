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


# Reactome ----------------------------------------------------------------

# Function to query, get genes and pathways
query_to_pathwaysREAC <- function(x) {
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
    Reactome_res.ALS <- ReactomePA::enrichPathway(gene = annot_als$Gene, pvalueCutoff = 0.01,
                                                  organism = "human", pAdjustMethod = "BH",
                                                  qvalueCutoff = 0.01, readable = T)
    
    path <- dplyr::tibble('Pathway' = Reactome_res.ALS@result[["Description"]],
                          'q_value' = Reactome_res.ALS@result[["qvalue"]]) %>% 
      dplyr::filter(q_value < 0.01)
    
    path$coviddi <- stringr::str_detect(path$Pathway, 
                                        stringr::regex("COVID|coronavirus|SARS-CoV-2|COVID-19", 
                                                       ignore_case = TRUE))
    
    somma <- sum(path$coviddi, na.rm = TRUE)
    matched_pathways <- path$Pathway[path$coviddi]
    
  }, error = function(e) {
    warning(paste("An error occurred for the term", x, ": ", e))
  })
  
  return(list("sum" = somma, "matched_pathways" = matched_pathways))
}

# Function to execute pipeline for a given seed
execute_pipelineREAC <- function(seed) {
  set.seed(seed)
  terms_2 <- sample(terms, size = 2, replace = FALSE)
  
  results <- map(terms_2, query_to_pathwaysREAC, .progress = "progress")
  
  somma_list <- map(results, 'sum')
  matched_pathways_list <- map(results, 'matched_pathways')
  
  pathways_df <- dplyr::tibble('term' = terms_2, 
                               'sum' = somma_list, 
                               'matched_pathways' = matched_pathways_list,
                               'seed' = seed)
  
  return(pathways_df)
}

# Execute pipeline for each seed and bind rows
final_results_dfREAC <- map_dfr(seeds, execute_pipelineREAC, .progress = "progress")

# If 'matched_pathways' is a list column in final_results_df
final_results_dfREAC$matched_pathways <- map_chr(final_results_dfREAC$matched_pathways, function(x) {
  paste(x, collapse = ",")
})

final_results_dfREAC$sum <- map_chr(final_results_dfREAC$sum, ~ {
  if (is.na(.x)) return("NA")
  return(as.character(.x))
})

write.csv(final_results_dfREAC, "final_results_dfREAC.csv")

save(final_results_dfREAC, file = 'final_resutsREAC.rda')
