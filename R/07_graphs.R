library(tidyverse)

# Assuming your data is in a dataframe called df
load('R/final_resutsREAC_pval.rda')

df <- final_results_dfREAC

# Calculate means and standard errors for q_value and geneCOUNT
summary_df <- df %>%
  group_by(coviddi) %>%
  summarise(
    avg_q_value = -log2(mean(q_value)),
    avg_gene_count = mean(geneCOUNT),
    se_q_value = sd(q_value) / sqrt(n()),
    se_gene_count = sd(geneCOUNT) / sqrt(n())
  )

# Plot for Average Q-value with error bars
q_value_plot <- ggplot(summary_df, aes(x = as.factor(coviddi), y = avg_q_value, fill = as.factor(coviddi))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = avg_q_value - se_q_value, ymax = avg_q_value + se_q_value),
                width = 0.2, position = position_dodge(.9)) +
  labs(x = "Related to COVID-19", y = "-log2(Average Q-value)", fill = "COVID-19 Related") +
  theme_minimal()

# Plot for Average Gene Count with error bars
gene_count_plot <- ggplot(summary_df, aes(x = as.factor(coviddi), y = avg_gene_count, fill = as.factor(coviddi))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = avg_gene_count - se_gene_count, ymax = avg_gene_count + se_gene_count),
                width = 0.2, position = position_dodge(.9)) +
  labs(x = "Related to COVID-19", y = "Average Gene Count", fill = "COVID-19 Related") +
  theme_minimal()

# Print the plots
print(q_value_plot)
print(gene_count_plot)

# If you want to save the plots
ggsave("q_value_plot.png", q_value_plot, width = 8, height = 4, dpi = 600)
ggsave("gene_count_plot.png", gene_count_plot, width = 8, height = 4, dpi = 600)


