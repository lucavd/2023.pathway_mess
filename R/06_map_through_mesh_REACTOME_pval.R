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


# Reactome ----------------------------------------------------------------------

# Function to query, get genes and pathways
query_to_pathwaysREAC <- function(x) {
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
      Reactome_res.ALS <- ReactomePA::enrichPathway(gene = annot_als$Gene, pvalueCutoff = 0.01,
                                                    organism = "human", pAdjustMethod = "BH",
                                                    qvalueCutoff = 0.01, readable = TRUE)
      
      # Check if the results are non-empty before constructing the tibble
      if (!is.null(Reactome_res.ALS) && !is.null(Reactome_res.ALS@result) && nrow(Reactome_res.ALS@result) > 0) {
        path <- tibble(
          PathwayID = Reactome_res.ALS@result[["ID"]],
          Pathway = Reactome_res.ALS@result[["Description"]],
          q_value = Reactome_res.ALS@result[["qvalue"]],
          geneCOUNT = Reactome_res.ALS@result[["Count"]],
          genes = Reactome_res.ALS@result[["geneID"]]
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
    # Return the pre-initialized empty `path` to ensure the return value is always a tibble
    return(path)
  })
}


# Function to execute pipeline for a given seed
execute_pipelineREAC <- function(seed) {
  set.seed(seed)
  terms_2 <- sample(terms, size = 5004, replace = FALSE)
  
  results <- map(terms_2, query_to_pathwaysREAC, .progress = "progress")
  
  return(results)
}

# Execute pipeline for each seed and bind rows
final_results_dfREAC <- map_dfr(seeds, execute_pipelineREAC, .progress = "progress")

write.csv(final_results_dfREAC, "final_results_dfREAC_pval.csv")

save(final_results_dfREAC, file = 'final_resutsREAC_pval.rda')