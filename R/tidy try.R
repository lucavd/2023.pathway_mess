library(tidyverse)

# Assuming final_results_dfREAC is your initial dataframe
df <- final_results_dfREAC

# Filter rows where coviddi is TRUE
covid_related <- df %>% filter(coviddi)

# Sample 994 unique terms
random_terms <- sample(covid_related$terms, 994, replace = FALSE)

# Filter rows with selected random terms
df_random_terms <- df %>% filter(terms %in% random_terms) %>% 
  mutate(terms = as.factor(terms))

# Function to get positions of COVID-related pathways for each term
get_positions <- function(term) {
  term_data <- df_random_terms %>% filter(terms == term)
  positions <- term_data %>% filter(coviddi) %>% pull(rowid)
  tibble(term = term, type = "Position", value = list(positions))
}

# Function to get summary of q-values for each term
get_q_values <- function(term) {
  term_data <- df_random_terms %>% filter(terms == term)
  q_values <- summary(term_data$q_value)
  tibble(term = term, type = "Q Value Summary", value = list(q_values))
}

# Function to get summary of COVID-related q-values for each term
get_covid_q_values <- function(term) {
  term_data <- df_random_terms %>% filter(terms == term, coviddi)
  q_values <- summary(term_data$q_value)
  tibble(term = term, type = "COVID Q Value Summary", value = list(q_values))
}

# Create separate dataframes
positions_df <- map_df(levels(df_random_terms$terms), get_positions)
q_values_df <- map_df(levels(df_random_terms$terms), get_q_values)
covid_q_values_df <- map_df(levels(df_random_terms$terms), get_covid_q_values)

# Combine the dataframes
combined_results <- reduce(list(positions_df, q_values_df, covid_q_values_df), full_join, by = "term")

# View the combined results
print(combined_results)
