library(tidyverse)

# Reactome ----------------------------------------------------------------

load("final_resutsREAC_pval.rda")
df <- final_results_dfREAC

covid_related <- df %>% filter(coviddi)

df_terms <- df %>% filter(terms %in% covid_related$terms) %>% 
  mutate(terms = as.factor(terms))

# Function to get positions of COVID-related pathways for each term
get_positions <- function(term) {
  term_data <- df_terms %>% filter(terms == term)
  positions <- term_data %>% filter(coviddi) %>% pull(rowid)
  tibble(term = term, type = "Position", value = list(positions))
}

positions_df_REAC <- map_df(levels(df_terms$terms), get_positions)

# whether each positive match term has covid-related pathways 
# within its top 20 and 100 
positions_df_REAC <- positions_df_REAC %>%
  rowwise() %>%
  mutate(within_top_20 = ifelse(any(value < 20), 1, 0))

positions_df_REAC <- positions_df_REAC %>%
  rowwise() %>%
  mutate(within_top_100 = ifelse(any(value < 100), 1, 0))

### repeat previous command only for those terms that have >20 pathways

max_pathways <- df %>% 
  filter(terms %in% covid_related$terms) %>%
  group_by(terms) %>%
  summarise(max_id = max(rowid))

terms_morequal_20 <- max_pathways %>% 
  filter(max_id >= 20)

df_terms <- df %>% filter(terms %in% terms_morequal_20$terms) %>% 
  mutate(terms = as.factor(terms))

get_positions20 <- function(term) {
  term_data <- df_terms %>% filter(terms == term)
  positions <- term_data %>% filter(coviddi) %>% pull(rowid)
  tibble(term = term, type = "Position", value = list(positions))
}

positions20_df_REAC <- map_df(levels(df_terms$terms), get_positions20)

positions20_df_REAC <- positions20_df_REAC %>%
  rowwise() %>%
  mutate(within_top_20 = ifelse(any(value < 20), 1, 0))

### repeat previous command only for those terms that have >100 pathways

terms_morequal_100 <- max_pathways %>% 
  filter(max_id >= 100)

df_terms <- df %>% filter(terms %in% terms_morequal_100$terms) %>% 
  mutate(terms = as.factor(terms))

get_positions100 <- function(term) {
  term_data <- df_terms %>% filter(terms == term)
  positions <- term_data %>% filter(coviddi) %>% pull(rowid)
  tibble(term = term, type = "Position", value = list(positions))
}

positions100_df_REAC <- map_df(levels(df_terms$terms), get_positions100)

positions100_df_REAC <- positions100_df_REAC %>%
  rowwise() %>%
  mutate(within_top_100 = ifelse(any(value < 100), 1, 0))


# KEGG --------------------------------------------------------------------

load("final_resutsKEGG_pval.rda")
df <- final_results_dfKEGG

covid_related <- df %>% filter(coviddi)

df_terms <- df %>% filter(terms %in% covid_related$terms) %>% 
  mutate(terms = as.factor(terms))

# Function to get positions of COVID-related pathways for each term
get_positions <- function(term) {
  term_data <- df_terms %>% filter(terms == term)
  positions <- term_data %>% filter(coviddi) %>% pull(rowid)
  tibble(term = term, type = "Position", value = list(positions))
}

positions_df_KEGG <- map_df(levels(df_terms$terms), get_positions)

# whether each positive match term has covid-related pathways 
# within its top 20 and 100 
positions_df_KEGG <- positions_df_KEGG %>%
  rowwise() %>%
  mutate(within_top_20 = ifelse(any(value < 20), 1, 0))

positions_df_KEGG <- positions_df_KEGG %>%
  rowwise() %>%
  mutate(within_top_100 = ifelse(any(value < 100), 1, 0))

### repeat previous command only for those terms that have >20 pathways

max_pathways <- df %>% 
  filter(terms %in% covid_related$terms) %>%
  group_by(terms) %>%
  summarise(max_id = max(rowid))

terms_morequal_20 <- max_pathways %>% 
  filter(max_id >= 20)

df_terms <- df %>% filter(terms %in% terms_morequal_20$terms) %>% 
  mutate(terms = as.factor(terms))

get_positions20 <- function(term) {
  term_data <- df_terms %>% filter(terms == term)
  positions <- term_data %>% filter(coviddi) %>% pull(rowid)
  tibble(term = term, type = "Position", value = list(positions))
}

positions20_df_KEGG <- map_df(levels(df_terms$terms), get_positions20)

positions20_df_KEGG <- positions20_df_KEGG %>%
  rowwise() %>%
  mutate(within_top_20 = ifelse(any(value < 20), 1, 0))

### repeat previous command only for those terms that have >100 pathways

terms_morequal_100 <- max_pathways %>% 
  filter(max_id >= 100)

df_terms <- df %>% filter(terms %in% terms_morequal_100$terms) %>% 
  mutate(terms = as.factor(terms))

get_positions100 <- function(term) {
  term_data <- df_terms %>% filter(terms == term)
  positions <- term_data %>% filter(coviddi) %>% pull(rowid)
  tibble(term = term, type = "Position", value = list(positions))
}

positions100_df_KEGG <- map_df(levels(df_terms$terms), get_positions100)

positions100_df_KEGG <- positions100_df_KEGG %>%
  rowwise() %>%
  mutate(within_top_100 = ifelse(any(value < 100), 1, 0))

# WP ----------------------------------------------------------------------

load("final_resutsWP_pval.rda")
df <- final_results_dfWP

covid_related <- df %>% filter(coviddi)

df_terms <- df %>% filter(terms %in% covid_related$terms) %>% 
  mutate(terms = as.factor(terms))

# Function to get positions of COVID-related pathways for each term
get_positions <- function(term) {
  term_data <- df_terms %>% filter(terms == term)
  positions <- term_data %>% filter(coviddi) %>% pull(rowid)
  tibble(term = term, type = "Position", value = list(positions))
}

positions_df_WP <- map_df(levels(df_terms$terms), get_positions)

# to check each positive match term has covid-related pathways 
# within its top 20 and 100 
positions_df_WP <- positions_df_WP %>%
  rowwise() %>%
  mutate(within_top_20 = ifelse(any(value < 20), 1, 0))

positions_df_WP <- positions_df_WP %>%
  rowwise() %>%
  mutate(within_top_100 = ifelse(any(value < 100), 1, 0))

### repeat previous command only for those terms that have >20 pathways

max_pathways <- df %>% 
  filter(terms %in% covid_related$terms) %>%
  group_by(terms) %>%
  summarise(max_id = max(rowid))

terms_morequal_20 <- max_pathways %>% 
  filter(max_id >= 20)

df_terms <- df %>% filter(terms %in% terms_morequal_20$terms) %>% 
  mutate(terms = as.factor(terms))

get_positions20 <- function(term) {
  term_data <- df_terms %>% filter(terms == term)
  positions <- term_data %>% filter(coviddi) %>% pull(rowid)
  tibble(term = term, type = "Position", value = list(positions))
}

positions20_df_WP <- map_df(levels(df_terms$terms), get_positions20)

positions20_df_WP <- positions20_df_WP %>%
  rowwise() %>%
  mutate(within_top_20 = ifelse(any(value < 20), 1, 0))

### repeat previous command only for those terms that have >100 pathways

terms_morequal_100 <- max_pathways %>% 
  filter(max_id >= 100)

df_terms <- df %>% filter(terms %in% terms_morequal_100$terms) %>% 
  mutate(terms = as.factor(terms))

get_positions100 <- function(term) {
  term_data <- df_terms %>% filter(terms == term)
  positions <- term_data %>% filter(coviddi) %>% pull(rowid)
  tibble(term = term, type = "Position", value = list(positions))
}

positions100_df_WP <- map_df(levels(df_terms$terms), get_positions100)

positions100_df_WP <- positions100_df_WP %>%
  rowwise() %>%
  mutate(within_top_100 = ifelse(any(value < 100), 1, 0))


