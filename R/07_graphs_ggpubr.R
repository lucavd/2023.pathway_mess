library(tidyverse)
library(ggpubr)

# BAR PLOT-- ggpubr -------------------------------------------------------

# REACTOME
load('final_resutsREAC_pval.rda')
dfREAC <- final_results_dfREAC

summary_dfREAC <- dfREAC %>%
  group_by(coviddi) %>%
  summarise(
    avg_q_value = -log2(mean(q_value)),
    avg_gene_count = mean(geneCOUNT),
    se_q_value = sd(q_value) / sqrt(n()),
    se_gene_count = sd(geneCOUNT) / sqrt(n())
  )
summary_dfREAC$coviddi <-  as.factor(summary_dfREAC$coviddi)
summary_dfREAC$database <- as.factor(c("Reactome", "Reactome"))

q_value_barplot <- ggbarplot(summary_dfREAC, x = "coviddi", y = "avg_q_value",
                             fill = "coviddi",
                             color = "white",
                             palette = "plos",
                             width = 0.3,
                             format.scale = TRUE,
                             xlab = "\nPathways Related to COVID-19",
                             ylab ="-log2(Average q-value)\n",
                             font.x = c(12, "bold", "black"),
                             font.y = c(12, "bold", "black"),
                             x.text.angle = 90,
                             legend = "none"
)
q_value_barplot <- q_value_barplot + 
  geom_errorbar(data = NULL,
                aes(ymin = avg_q_value - se_q_value, ymax = avg_q_value + se_q_value),
                width = 0.1,
                stat = "identity",
                position = "identity",
                inherit.aes = TRUE)

print(q_value_barplot)
ggsave("q_value_barplot_REAC.png", q_value_barplot, width = 6, height = 6, dpi = 600)


gene_count_barplot <- ggbarplot(summary_dfREAC, x = "coviddi", y = "avg_gene_count",
                                fill = "coviddi",
                                color = "white",
                                palette = "plos",
                                width = 0.3,
                                format.scale = TRUE,
                                xlab = "\nPathways Related to COVID-19",
                                ylab ="Average Gene Count\n",
                                font.x = c(12, "bold", "black"),
                                font.y = c(12, "bold", "black"),
                                x.text.angle = 90,
                                legend = "none"
)
gene_count_barplot <- gene_count_barplot +
  geom_errorbar(data = NULL,
                aes(ymin = avg_gene_count - se_gene_count, ymax = avg_gene_count + se_gene_count),
                width = 0.2,
                stat = "identity",
                position = "identity",
                inherit.aes = TRUE)

print(gene_count_barplot)
ggsave("gene_count_barplot_REAC.png", gene_count_barplot, width = 6, height = 6, dpi = 600)


# KEGG
load('final_resutsKEGG_pval.rda')
dfKEGG <- final_results_dfKEGG

summary_dfKEGG <- dfKEGG %>%
  group_by(coviddi) %>%
  summarise(
    avg_q_value = -log2(mean(q_value)),
    avg_gene_count = mean(geneCOUNT),
    se_q_value = sd(q_value) / sqrt(n()),
    se_gene_count = sd(geneCOUNT) / sqrt(n())
  )
summary_dfKEGG$coviddi <-  as.factor(summary_dfKEGG$coviddi)
summary_dfKEGG$database <- as.factor(c("KEGG", "KEGG"))

q_value_barplot <- ggbarplot(summary_dfKEGG, x = "coviddi", y = "avg_q_value",
                             fill = "coviddi",
                             color = "white",
                             palette = "plos",
                             width = 0.3,
                             format.scale = TRUE,
                             xlab = "\nPathways Related to COVID-19",
                             ylab ="-log2(Average q-value)\n",
                             font.x = c(12, "bold", "black"),
                             font.y = c(12, "bold", "black"),
                             x.text.angle = 90,
                             legend = "none"
)
q_value_barplot <- q_value_barplot + 
  geom_errorbar(data = NULL,
                aes(ymin = avg_q_value - se_q_value, ymax = avg_q_value + se_q_value),
                width = 0.1,
                stat = "identity",
                position = "identity",
                inherit.aes = TRUE)

print(q_value_barplot)
ggsave("q_value_barplot_KEGG.png", q_value_barplot, width = 6, height = 6, dpi = 600)


gene_count_barplot <- ggbarplot(summary_dfKEGG, x = "coviddi", y = "avg_gene_count",
                                fill = "coviddi",
                                color = "white",
                                palette = "plos",
                                width = 0.3,
                                format.scale = TRUE,
                                xlab = "\nPathways Related to COVID-19",
                                ylab ="Average Gene Count\n",
                                font.x = c(12, "bold", "black"),
                                font.y = c(12, "bold", "black"),
                                x.text.angle = 90,
                                legend = "none"
)
gene_count_barplot <- gene_count_barplot +
  geom_errorbar(data = NULL,
                aes(ymin = avg_gene_count - se_gene_count, ymax = avg_gene_count + se_gene_count),
                width = 0.2,
                stat = "identity",
                position = "identity",
                inherit.aes = TRUE)

print(gene_count_barplot)
ggsave("gene_count_barplot_KEGG.png", gene_count_barplot, width = 6, height = 6, dpi = 600)


# WP
load('final_resutsWP_pval.rda')
dfWP <- final_results_dfWP

summary_dfWP <- dfWP %>%
  group_by(coviddi) %>%
  summarise(
    avg_q_value = -log2(mean(q_value)),
    avg_gene_count = mean(geneCOUNT),
    se_q_value = sd(q_value) / sqrt(n()),
    se_gene_count = sd(geneCOUNT) / sqrt(n())
  )
summary_dfWP$coviddi <-  as.factor(summary_dfWP$coviddi)
summary_dfWP$database <- as.factor(c("WP", "WP"))

q_value_barplot <- ggbarplot(summary_dfWP, x = "coviddi", y = "avg_q_value",
                             fill = "coviddi",
                             color = "white",
                             palette = "plos",
                             width = 0.3,
                             format.scale = TRUE,
                             xlab = "\nPathways Related to COVID-19",
                             ylab ="-log2(Average q-value)\n",
                             font.x = c(12, "bold", "black"),
                             font.y = c(12, "bold", "black"),
                             x.text.angle = 90,
                             legend = "none"
)
q_value_barplot <- q_value_barplot + 
  geom_errorbar(data = NULL,
                aes(ymin = avg_q_value - se_q_value, ymax = avg_q_value + se_q_value),
                width = 0.1,
                stat = "identity",
                position = "identity",
                inherit.aes = TRUE)

print(q_value_barplot)
ggsave("q_value_barplot_WP.png", q_value_barplot, width = 6, height = 6, dpi = 600)

gene_count_barplot <- ggbarplot(summary_dfWP, x = "coviddi", y = "avg_gene_count",
                                fill = "coviddi",
                                color = "white",
                                palette = "plos",
                                width = 0.3,
                                format.scale = TRUE,
                                xlab = "\nPathways Related to COVID-19",
                                ylab ="Average Gene Count\n",
                                font.x = c(12, "bold", "black"),
                                font.y = c(12, "bold", "black"),
                                x.text.angle = 90,
                                legend = "none"
)
gene_count_barplot <- gene_count_barplot +
  geom_errorbar(data = NULL,
                aes(ymin = avg_gene_count - se_gene_count, ymax = avg_gene_count + se_gene_count),
                width = 0.2,
                stat = "identity",
                position = "identity",
                inherit.aes = TRUE)

print(gene_count_barplot)
ggsave("gene_count_barplot_WP.png", gene_count_barplot, width = 6, height = 6, dpi = 600)


# MERGED ------------------------------------------------------------------

load('final_resutsREAC_pval.rda')
dfREAC <- final_results_dfREAC

summary_dfREAC <- dfREAC %>%
  group_by(coviddi) %>%
  summarise(
    avg_q_value = -log2(mean(q_value)),
    avg_gene_count = mean(geneCOUNT),
    se_q_value = sd(q_value) / sqrt(n()),
    se_gene_count = sd(geneCOUNT) / sqrt(n())
  )
summary_dfREAC$coviddi <-  as.factor(summary_dfREAC$coviddi)
summary_dfREAC$database <- as.factor(c("Reactome", "Reactome"))


load('final_resutsKEGG_pval.rda')
dfKEGG <- final_results_dfKEGG

summary_dfKEGG <- dfKEGG %>%
  group_by(coviddi) %>%
  summarise(
    avg_q_value = -log2(mean(q_value)),
    avg_gene_count = mean(geneCOUNT),
    se_q_value = sd(q_value) / sqrt(n()),
    se_gene_count = sd(geneCOUNT) / sqrt(n())
  )
summary_dfKEGG$coviddi <-  as.factor(summary_dfKEGG$coviddi)
summary_dfKEGG$database <- as.factor(c("KEGG", "KEGG"))


load('final_resutsWP_pval.rda')
dfWP <- final_results_dfWP

summary_dfWP <- dfWP %>%
  group_by(coviddi) %>%
  summarise(
    avg_q_value = -log2(mean(q_value)),
    avg_gene_count = mean(geneCOUNT),
    se_q_value = sd(q_value) / sqrt(n()),
    se_gene_count = sd(geneCOUNT) / sqrt(n())
  )
summary_dfWP$coviddi <-  as.factor(summary_dfWP$coviddi)
summary_dfWP$database <- as.factor(c("WP", "WP"))

## barplot
summary_databases <- rbind(summary_dfREAC, summary_dfKEGG, summary_dfWP)

q_value_barplot <- ggbarplot(summary_databases, x = "coviddi", y = "avg_q_value",
                             fill = "database",
                             color = "white",
                             palette = "plos",
                             width = 0.4,
                             position = position_dodge(0.5),
                             format.scale = TRUE,
                             xlab = "\nPathways Related to COVID-19",
                             ylab ="-log2(Average q-value)\n",
                             font.x = c(12, "bold", "black"),
                             font.y = c(12, "bold", "black"),
                             x.text.angle = 90,
                             legend = "none",
                             #add = "mean_se",
                             #error.plot = "errorbar"
)
q_value_barplot <- q_value_barplot +
  geom_errorbar(data = NULL,
                aes(ymin = avg_q_value - se_q_value, ymax = avg_q_value + se_q_value),
                width = 0.1,
                stat = "identity",
                position = position_dodge2(width = 0.1, padding = 0.1),
                inherit.aes = TRUE)

print(q_value_barplot)
