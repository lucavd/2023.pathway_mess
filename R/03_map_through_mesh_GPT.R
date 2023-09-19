# setup -------------------------------------------------------------------
#library(XML) # Uncomment if you use 'XML'
library(RISmed)
library(tidyverse)
library(ReactomePA)
library(org.Hs.eg.db)
#library(here) # Uncomment if you use 'here'

# databases ---------------------------------------------------------------
mesh <- XML::xmlToDataFrame(here::here("data/MeSH/desc2023.xml"))
terms <- mesh$DescriptorName
terms_10 <- sample(terms, size = 10, replace = FALSE)

pubtator <- read.delim(file = here::here("data/Pubtator/gene2pubtatorcentral"),
                       quote = "", header = TRUE,
                       col.names = c('PMID', 'Object', 'Gene', 'Gene_name', 
                                     'Dataset'))

# Function to Execute PubMed Query
execute_query <- function(x, start_date, end_date) {
  
  # Explicitly build the search query
  search_topic <- paste0('("', x, '"[Mesh]) AND ("', start_date, '"[Date - Entrez] : "', end_date, '"[Date - Entrez])')
  
  cat("Executing query: ", search_topic, "\n")
  
  # Debug: Drop into interactive mode
 # browser()
  
  # Run query and catch errors
  search_query <- tryCatch(
    {
      EUtilsSummary(search_topic, retmax = 1000)
    },
    error = function(e) {
      cat("Caught an error:", conditionMessage(e), "\n")
      return(NULL)
    },
    warning = function(w) {
      cat("Caught a warning:", conditionMessage(w), "\n")
    }
  )
  
  if (is.null(search_query) || !("ids" %in% slotNames(search_query))) {
    cat("No results or unexpected query object. Query was: ", search_topic, "\n")
    return(data.frame())
  }
  
  filtered_data <- data.frame('PMID' = as.numeric(PMID(EUtilsGet(search_query))))
  cat("Found ", nrow(filtered_data), " PMIDs for query: ", search_topic, "\n")
  
  return(filtered_data)
}




# Function for Pathway Enrichment Analysis
perform_pathway_analysis <- function(common_pmids, pubtator) {
  annot_common <- pubtator[pubtator$PMID %in% common_pmids,]
  
  enrich_result <- ReactomePA::enrichPathway(gene = annot_common$Gene,
                                             pvalueCutoff = 0.01,
                                             organism = "human",
                                             pAdjustMethod = "BH",
                                             qvalueCutoff = 0.01,
                                             readable = TRUE)
  
  if (is.null(enrich_result)) {
    return(tibble())
  }
  
  pathway_data <- as_tibble(enrich_result@result[, c("Description", "qvalue")])
  return(pathway_data)
}

# Main Function
main_function <- function(term, pubtator_data, start_date = "2018/01/01", end_date = "2018/12/31") {
  filtered_data <- execute_query(term, start_date, end_date)
  
  if (nrow(filtered_data) == 0) {
    return(tibble())
  }
  
  common_pmids <- intersect(pubtator_data$PMID, filtered_data$PMID)
  
  if (length(common_pmids) == 0) {
    return(tibble())
  }
  
  pathway_data <- perform_pathway_analysis(common_pmids, pubtator_data)
  
  if (nrow(pathway_data) == 0) {
    return(tibble())
  }
  
  pathway_data <- pathway_data %>% filter(qvalue < 0.01)
  pathway_data$covid_related <- str_detect(pathway_data$Description, regex("COVID|coronavirus|SARS-CoV-2|COVID-19", ignore_case = TRUE))
  
  sum_covid_related <- sum(pathway_data$covid_related)
  
  return(tibble('Mesh_Term' = term,
                'Num_COVID_Related_Pathways' = sum_covid_related))
}

# Execute the main function for each term in 'terms_10'
final_results <- map_dfr(terms_10, function(term) {
  
  cat("Processing term:", term, "\n")
  
  filtered_data <- execute_query(term, "2018/01/01", "2018/12/31")
  
  if (nrow(filtered_data) == 0) {
    cat("No PMIDs found for term:", term, "\n")
    return(tibble())
  }
  
  common_pmids <- intersect(pubtator$PMID, filtered_data$PMID)
  
  if (length(common_pmids) == 0) {
    cat("No overlapping PMIDs for term:", term, "\n")
    return(tibble())
  }
  
  pathway_data <- perform_pathway_analysis(common_pmids, pubtator)
  
  if (nrow(pathway_data) == 0) {
    cat("No significant pathway enrichment for term:", term, "\n")
    return(tibble())
  }
  
  pathway_data <- pathway_data %>% filter(qvalue < 0.01)
  pathway_data$covid_related <- str_detect(pathway_data$Description, regex("COVID|coronavirus|SARS-CoV-2|COVID-19", ignore_case = TRUE))
  
  sum_covid_related <- sum(pathway_data$covid_related)
  
  result_tibble <- tibble('Mesh_Term' = term,
                          'Num_COVID_Related_Pathways' = sum_covid_related)
  
  return(result_tibble)
})

cat("Final Results:\n")
print(final_results)

