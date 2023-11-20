library(tidyverse)


# Reactome ----------------------------------------------------------------

load("final_resutsREAC_pval.rda")
df <- final_results_dfREAC

covid_related <- df[df$coviddi == TRUE, ]

random_terms <- sample(covid_related$terms, 10, replace = FALSE)
print(random_terms)

df_random_terms <- df[df$terms %in% random_terms, ]

df_random_terms$terms <- as.factor(df_random_terms$terms)


### covid pathways' POSITIONS for 10 MeSH terms (based on row ids)

for (i in levels(df_random_terms$terms)) {
  
  term <- df_random_terms[df_random_terms$terms == i, ]
  
  p <- term[term$coviddi == TRUE, "rowid"]  
  
  print(paste("Position of COVID-related pathways for the MeSH term", i, ":", p))
}

### covid pathway's significance (summary of the q-values)

## for all pathways
for (i in levels(df_random_terms$terms)) {
  
  term <- df_random_terms[df_random_terms$terms == i, ]
  
  q <- summary(term$q_value)
  
  print(paste("Summary of q-values for the MeSH term", i, ":", q))
}

## for covid pathways
for (i in levels(df_random_terms$terms)) {
  
  term <- df_random_terms[df_random_terms$terms == i, ]
  
  covid_in_term <- term[term$coviddi == TRUE, ]
  
  q <- summary(covid_in_term$q_value)
  
  print(paste("Summary of COVID-related q-values for the MeSH term", i, ":", q))
}


# KEGG --------------------------------------------------------------------

load("final_resutsKEGG_pval.rda")
df <- final_results_dfKEGG

covid_related <- df[df$coviddi == TRUE, ]

random_terms <- sample(covid_related$terms, 10, replace = FALSE)
print(random_terms)

df_random_terms <- df[df$terms %in% random_terms, ]

df_random_terms$terms <- as.factor(df_random_terms$terms)


### covid pathways' POSITIONS for 10 MeSH terms (based on row ids)

for (i in levels(df_random_terms$terms)) {
  
  term <- df_random_terms[df_random_terms$terms == i, ]
  
  p <- term[term$coviddi == TRUE, "rowid"]  
  
  print(paste("Position of COVID-related pathways for the MeSH term", i, ":", p))
}

### covid pathway's significance (summary of the q-values)

## for all pathways
for (i in levels(df_random_terms$terms)) {
  
  term <- df_random_terms[df_random_terms$terms == i, ]
  
  q <- summary(term$q_value)
  
  print(paste("Summary of q-values for the MeSH term", i, ":", q))
}

## for covid pathways
for (i in levels(df_random_terms$terms)) {
  
  term <- df_random_terms[df_random_terms$terms == i, ]
  
  covid_in_term <- term[term$coviddi == TRUE, ]
  
  q <- summary(covid_in_term$q_value)
  
  print(paste("Summary of COVID-related q-values for the MeSH term", i, ":", q))
}


# WP ----------------------------------------------------------------------

load("final_resutsWP_pval.rda")
df <- final_results_dfWP

covid_related <- df[df$coviddi == TRUE, ]

random_terms <- sample(covid_related$terms, 10, replace = FALSE)
print(random_terms)

df_random_terms <- df[df$terms %in% random_terms, ]

df_random_terms$terms <- as.factor(df_random_terms$terms)


### covid pathways' POSITIONS for 10 MeSH terms (based on row ids)

for (i in levels(df_random_terms$terms)) {
  
  term <- df_random_terms[df_random_terms$terms == i, ]
  
  p <- term[term$coviddi == TRUE, "rowid"]  
  
  print(paste("Position of COVID-related pathways for the MeSH term", i, ":", p))
}

### covid pathway's significance (summary of the q-values)

## for all pathways
for (i in levels(df_random_terms$terms)) {
  
  term <- df_random_terms[df_random_terms$terms == i, ]
  
  q <- summary(term$q_value)
  
  print(paste("Summary of q-values for the MeSH term", i, ":", q))
}

## for covid pathways
for (i in levels(df_random_terms$terms)) {
  
  term <- df_random_terms[df_random_terms$terms == i, ]
  
  covid_in_term <- term[term$coviddi == TRUE, ]
  
  q <- summary(covid_in_term$q_value)
  
  print(paste("Summary of COVID-related q-values for the MeSH term", i, ":", q))
}




# EXAMPLE -----------------------------------------------------------------

## for position: flu
flu <- df_random_terms[df_random_terms$terms == "Influenza, Human", ]

p <- flu[flu$coviddi == TRUE, "rowid"]


## for summary of q-values: colitis
Colitis <- df_random_terms[df_random_terms$terms == "Colitis", ]

summary(Colitis$q_value)

Colitis_covid <- Colitis[Colitis$coviddi == TRUE, ]

summary(Colitis_covid$q_value)
