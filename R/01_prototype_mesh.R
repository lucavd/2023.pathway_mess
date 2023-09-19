library(xml2)
mesh <- read_xml(here::here("data/MeSH/desc2023.xml"))

mesh <- XML::xmlToDataFrame(here::here("data/MeSH/desc2023.xml"))

terms <- mesh$DescriptorName

paste0(terms[1])

query_string <- paste0('("', terms[1], '"[Mesh]) AND ("2018/01/01"[Date - Entrez] : "2018/12/31"[Date - Entrez])')


# no more of 3 requests per second or during weekends or nights (see RISmed help)

search_topic <- '("Amyotrophic Lateral Sclerosis"[Mesh]) AND ("2018/01/01"[Date - Entrez] : "2018/12/31"[Date - Entrez])'

search_query <- EUtilsSummary(search_topic, retmax = 10000)
#search_query <- EUtilsSummary(query_string, retmax = 10000)

summary(search_query) # check if more than 10000

pubmed_data <- data.frame('Title' = ArticleTitle(EUtilsGet(search_query)),
                          'Abstract' = AbstractText(EUtilsGet(search_query)),
                          'PMID'= PMID(EUtilsGet(search_query)))

# pubmed_data$Title[1] <- 'COVID' # test for the following code

filtered_data <- pubmed_data %>%
  filter(
    !str_detect(Title, regex("COVID|coronavirus|SARS-CoV-2|COVID-19", 
                             ignore_case = TRUE)),
    !str_detect(Abstract, regex("COVID|coronavirus|SARS-CoV-2|COVID-19", 
                                ignore_case = TRUE))
  )

# could be run once at the beginning
pubtator <- read.delim(file = here::here("data/Pubtator/gene2pubtatorcentral"),
                       quote = "", header = TRUE,
                       col.names = c('PMID', 'Object', 'Gene', 'Gene_name', 
                                     'Dataset'))


common_als <- (intersect(pubtator$PMID, as.numeric(filtered_data$PMID)))
annot_als <- pubtator[pubtator$PMID %in% common_als,]

Reactome_res.ALS <- enrichPathway(gene= annot_als$Gene,pvalueCutoff= 0.01,
                                  organism = "human", pAdjustMethod = "BH", 
                                  qvalueCutoff = 0.05, readable=T)

path <- tibble('Pathway' = Reactome_res.ALS@result[["Description"]],
               'q_value' = Reactome_res.ALS@result[["qvalue"]])

path$coviddi <- str_detect(path$Pathway, regex("COVID|coronavirus|SARS-CoV-2|COVID-19", 
                                               ignore_case = TRUE))

sum(path$coviddi)
