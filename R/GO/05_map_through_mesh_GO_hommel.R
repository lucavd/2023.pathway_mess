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

# GO ----------------------------------------------------------------------

# Function to query, get genes and pathways
query_to_pathwaysGO <- function(x) {
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
    GO_res <- enrichGO(gene = annot_als$Gene, OrgDb = org.Hs.eg.db , 
                       keyType = "ENTREZID", ont = "ALL", 
                       pvalueCutoff = 0.01, pAdjustMethod = "hommel", qvalueCutoff = 0.01,
                       readable = TRUE, pool = TRUE)
    
    path <- dplyr::tibble('Pathway' = GO_res@result[["Description"]],
                          'q_value' = GO_res@result[["qvalue"]]) %>% 
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
execute_pipelineGO <- function(seed) {
  set.seed(seed)
  terms_2 <- sample(terms, size = 10, replace = FALSE)
  
  results <- map(terms_2, query_to_pathwaysGO, .progress = "progress")
  
  somma_list <- map(results, 'sum')
  matched_pathways_list <- map(results, 'matched_pathways')
  
  pathways_df <- dplyr::tibble('term' = terms_2, 
                               'sum' = somma_list, 
                               'matched_pathways' = matched_pathways_list,
                               'seed' = seed)
  
  return(pathways_df)
}

# Execute pipeline for each seed and bind rows
final_results_dfGO <- map_dfr(seeds, execute_pipelineGO, .progress = "progress")

# If 'matched_pathways' is a list column in final_results_df
final_results_dfGO$matched_pathways <- map_chr(final_results_dfGO$matched_pathways, function(x) {
  paste(x, collapse = ",")
})

final_results_dfGO$sum <- map_chr(final_results_dfGO$sum, ~ {
  if (is.na(.x)) return("NA")
  return(as.character(.x))
})

write.csv(final_results_dfGO, "output/final_results_dfGO_hommel.csv")

save(final_results_dfGO, file = 'output/final_resutsGO_hommel.rda')

