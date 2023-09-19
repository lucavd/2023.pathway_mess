
# setup -------------------------------------------------------------------
#library(XML)
library(RISmed)
library(tidyverse)
library(ReactomePA)
library(org.Hs.eg.db)
#library(here)


# databases ---------------------------------------------------------------

mesh <- XML::xmlToDataFrame(here::here("data/MeSH/desc2023.xml"))

diseases = mesh |> filter(stringr::str_detect(TreeNumberList, stringr::regex("C", ignore_case = TRUE)))

terms <- mesh$DescriptorName

terms_100 <- sample(terms, size = 100, replace = FALSE)

pubtator <- read.delim(file = here::here("data/Pubtator/gene2pubtatorcentral"),
                       quote = "", header = TRUE,
                       col.names = c('PMID', 'Object', 'Gene', 'Gene_name', 
                                     'Dataset'))

# query -------------------------------------------------------------------

execute_query <- function(x) {
 
   search_topic <- paste0('("',x,'"[Mesh]) AND ("2018/01/01"[Date - Entrez] : "2018/12/31"[Date - Entrez])')
  
   search_query <- EUtilsSummary(search_topic, retmax = 1000)
}
  
get_genes <- function(x) {
  
  filtered_data <- data.frame('PMID'= (x@PMID))
  
  common_als <- (intersect(pubtator$PMID, as.numeric(filtered_data$PMID)))
  
  annot_als <- pubtator[pubtator$PMID %in% common_als,]
  
} 

get_pathways <- function(x) {
  # Initialize somma as NA to handle cases where no sum can be computed
  somma <- NA 
  
  # Exception handling
  tryCatch({
    Reactome_res.ALS <- ReactomePA::enrichPathway(gene= x$Gene, pvalueCutoff= 0.01,
                                                  organism = "human", pAdjustMethod = "BH", 
                                                  qvalueCutoff = 0.01, readable=T)
    
    path <- dplyr::tibble('Pathway' = Reactome_res.ALS@result[["Description"]],
                          'q_value' = Reactome_res.ALS@result[["qvalue"]]) %>% 
      dplyr::filter(q_value < 0.01)
    
    path$coviddi <- stringr::str_detect(path$Pathway, 
                                        stringr::regex("COVID|coronavirus|SARS-CoV-2|COVID-19", 
                                                       ignore_case = TRUE))
    
    somma <- sum(path$coviddi, na.rm = TRUE)
    
  }, error = function(e) {
    warning("An error occurred: ", e)
  })
  
  return(somma)
}



results <- map(terms_100, execute_query, .progress = "progress")
genes <- map(results, get_genes, .progress = "progress")
pathways <- map(genes, get_pathways, .progress = "progress")

names(pathways) <- terms_100

pathways_df <- dplyr::tibble('term' = terms_100, 'somma' = purrr::flatten_dbl(pathways))
