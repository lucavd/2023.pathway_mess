library(RISmed)
library(ReactomePA)
library(org.Hs.eg.db)
library(tidyverse)

devtools::install_github('skoval/RISmed')

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install('DOSE')
# BiocManager::install('ReactomePA')
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("ReactomePA")

########################################################################
###############' LITERATURE RETRIEVAL
# modify queries, add parenthesis, change date and correct "

# ALS -----------------------------------------------------------

search_topic <- '("amyotrophic lateral sclerosis" [TIAB]) AND ("biomarkers"[MeSH Terms] OR "biomarkers"[TIAB] OR "biomarker"[TIAB]) AND ("2015/01/01"[Date - Entrez] : "2023/06/13"[Date - Entrez])'

search_query <- EUtilsSummary(search_topic, retmax = 10000)
summary(search_query)
pubmed_data_als <- search_query@PMID



# CAD -------------------------------------------------------------

search_topic <- '("CAD"[TIAB] OR "coronary artery disease"[TIAB]) AND ("biomarkers"[MeSH Terms]OR "biomarkers"[TIAB] OR "biomarker"[TIAB]) AND ("2015/01/01"[Date - Entrez] : "2023/06/13"[Date - Entrez])'

search_query <- EUtilsSummary(search_topic, retmax = 10000)
summary(search_query)
pubmed_data_cad <- search_query@PMID



# cardiomiopathy ----------------------------------------------------

search_topic <- '("cardiomyopathy"[TIAB]) AND ("biomarkers"[MeSH Terms] OR "biomarkers"[TIAB] OR "biomarker"[TIAB]) AND ("2015/01/01"[Date - Entrez] : "2023/06/13"[Date - Entrez])'

search_query <- EUtilsSummary(search_topic, retmax = 10000)
summary(search_query)
pubmed_data_cardiomio <- search_query@PMID


# ARRHYTHMIA ----------------------------------------------

search_topic <- '("arrhythmia"[TIAB]) AND ("biomarkers"[MeSH Terms] OR"biomarkers"[TIAB] OR "biomarker"[TIAB]) AND ("2015/01/01"[Date - Entrez] : "2023/06/13"[Date - Entrez])'

search_query <- EUtilsSummary(search_topic, retmax = 10000)
summary(search_query)
pubmed_data_arr <- search_query@PMID


# HEART FAILURE ----------------------------------------------

search_topic <- '("heart failure"[TIAB]) AND ("biomarkers"[MeSH Terms] OR "biomarkers"[TIAB] OR "biomarker"[TIAB]) AND ("2015/01/01"[Date - Entrez] : " 2023/06/13"[Date - Entrez])'
# max retrievable are 9999 
search_query <- EUtilsSummary(search_topic, retmax = 10000)
summary(search_query)
pubmed_data_hf <- search_query@PMID


# RESPIRATORY FAILURE ----------------------------------------------

search_topic <- '("respiratory failure"[TIAB]) AND ("biomarkers"[MeSH Terms] OR "biomarkers"[TIAB] OR "biomarker"[TIAB]) AND ("2015/01/01"[Date - Entrez] : "2023/06/13"[Date - Entrez])'

search_query <- EUtilsSummary(search_topic, retmax = 10000)
summary(search_query)
pubmed_data_rf <- search_query@PMID


# PNEUMONIA ----------------------------------------------

search_topic <- '("pneumonia"[TIAB]) AND ("biomarkers"[MeSH Terms] OR "biomarkers"[TIAB] OR "biomarker"[TIAB]) AND ("2015/01/01"[Date - Entrez] : "2023/06/13"[Date - Entrez])'

search_query <- EUtilsSummary(search_topic, retmax = 10000)
summary(search_query)
pubmed_data_pneumo <- search_query@PMID


# METABOLIC SYNDROME ----------------------------------------------

search_topic <- '("metabolic syndrome"[TIAB]) AND ("biomarkers"[MeSH Terms] OR "biomarkers"[TIAB] OR "biomarker"[TIAB]) AND ("2015/01/01"[Date - Entrez] : "2023/06/13"[Date - Entrez])'

search_query <- EUtilsSummary(search_topic, retmax = 10000)
summary(search_query)
pubmed_data_met <- search_query@PMID


# DYSTROPHY ----------------------------------------------

search_topic <- '("dystrophy"[TIAB]) AND ("biomarkers"[MeSH Terms] OR "biomarkers"[TIAB] OR "biomarker"[TIAB]) AND ("2015/01/01"[Date - Entrez] : "2023/06/13"[Date - Entrez])'

search_query <- EUtilsSummary(search_topic, retmax = 10000)
summary(search_query)
pubmed_data_dys <- search_query@PMID

