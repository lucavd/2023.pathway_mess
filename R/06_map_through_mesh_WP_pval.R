# Load necessary libraries
library(RISmed)
library(tidyverse)
library(ReactomePA)
library(org.Hs.eg.db)
library(clusterProfiler)

# Read data
mesh <- XML::xmlToDataFrame(here::here("data/MeSH/desc2023.xml"))
diseases <- mesh |> filter(stringr::str_detect(TreeNumberList, stringr::regex("C", ignore_case = TRUE)))
terms <- diseases$DescriptorName
pubtator <- read.delim((here::here("data/Pubtator/gene2pubtatorcentral")),
                       quote = "", header = TRUE,
                       col.names = c('PMID', 'Object', 'Gene', 'Gene_name', 
                                     'Dataset'))

# Seeds
# seeds <- c(5, 10, 169, 3011, 19811)
seeds <- c(3011)


# WP ----------------------------------------------------------------------

# Function to query, get genes and pathways
query_to_pathwaysWP <- function(x) {
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
    WP_res <- enrichWP(gene = annot_als$Gene, organism = "Homo sapiens", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.01, qvalueCutoff = 0.01,)
    
    path <- dplyr::tibble('Pathway' = WP_res@result[["Description"]],
                          'q_value' = WP_res@result[["qvalue"]],
                          'genes'= WP_res@result[["geneID"]]) %>% 
      dplyr::filter(q_value < 0.01) %>% 
      rowid_to_column(var = "rowid") |> 
      mutate(term = paste0(x))
    
    path$coviddi <- stringr::str_detect(path$Pathway, 
                                        stringr::regex("COVID|coronavirus|SARS-CoV-2|COVID-19", 
                                                       ignore_case = TRUE))
    
  }, error = function(e) {
    warning(paste("An error occurred for the term", x, ": ", e))
  })
  
  return(path)
}

# Function to execute pipeline for a given seed
execute_pipelineWP <- function(seed) {
  set.seed(seed)
  terms_2 <- sample(terms, size = 5004, replace = FALSE)
  
  results <- map(terms_2, query_to_pathwaysWP, .progress = "progress")
  
  return(results)
}

# Execute pipeline for each seed and bind rows
final_results_dfWP <- map_dfr(seeds, execute_pipelineWP, .progress = "progress")

write.csv(final_results_dfWP, "final_results_dfWP_pval.csv")

save(final_results_dfWP, file = 'final_resutsWP_pval.rda')