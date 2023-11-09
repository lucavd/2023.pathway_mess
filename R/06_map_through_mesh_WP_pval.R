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
  # Initialize path as an empty tibble
  path <- tibble(
    PathwayID = character(),
    Pathway = character(),
    q_value = numeric(),
    geneCOUNT = integer(),
    genes = character(),
    coviddi = logical(),
    terms = character()
  )
  
  tryCatch({
    # Query execution
    search_topic <- paste0('("', x, '"[Mesh]) AND ("2018/01/01"[Date - Entrez] : "2018/12/31"[Date - Entrez])')
    search_query <- EUtilsSummary(search_topic, retmax = 1000)
    
    # Gene retrieval
    filtered_data <- data.frame('PMID' = (search_query@PMID))
    common_als <- intersect(pubtator$PMID, as.numeric(filtered_data$PMID))
    annot_als <- pubtator[pubtator$PMID %in% common_als,]
    
    # Pathway analysis
    if (length(common_als) > 0) {
      WP_res <- enrichWP(gene = annot_als$Gene, organism = "Homo sapiens", 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 0.01, qvalueCutoff = 0.01)
      if (!is.null(WP_res) && !is.null(WP_res@result) && nrow(WP_res@result) > 0) {
        path <- tibble(
          PathwayID = WP_res@result[["ID"]],
          Pathway = WP_res@result[["Description"]],
          q_value = WP_res@result[["qvalue"]],
          geneCOUNT = WP_res@result[["Count"]],
          genes = WP_res@result[["geneID"]]
        ) %>% 
          filter(q_value < 0.01) %>% 
          rowid_to_column(var = "rowid") %>% 
          mutate(terms = x)
        
        path$coviddi <- str_detect(path$Pathway, 
                                   regex("COVID|coronavirus|SARS-CoV-2|COVID-19", 
                                         ignore_case = TRUE))
      }
    }
    
    # Return the resulting path tibble (could be empty if no results or errors were encountered)
    return(path)
    
  }, error = function(e) {
    warning(paste("An error occurred for the term", x, ": ", e))
    # Return the pre-initialized empty `path` to ensure the return value is a tibble
    return(path)
  })
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