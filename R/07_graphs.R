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
library(tidyverse)
library(patchwork)

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

## merged summary databases
summary_databases <- rbind(summary_dfREAC, summary_dfKEGG, summary_dfWP)


#### barplot

## q-vqlue
q_value_barplot <- ggplot(summary_databases,
                          aes(x = coviddi, y= avg_q_value, fill = database)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_dodge(0.6)) +
  geom_errorbar(aes(ymin = avg_q_value - se_q_value, ymax = avg_q_value + se_q_value),
                width = 0.15,
                position = position_dodge(0.6)) +
  labs(#x = "\nPathways Related to COVID-19",
    x=" ",y = "-log2(Average q-value)\n",
       fill = "Database") +
  theme_classic() +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right"
        ) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    #axis.title.x = element_text(size = 14, face = "bold"), 
    axis.title.y = element_text(size = 14, face = "bold")
  ) 

print(q_value_barplot)
ggsave("q_value_barplot.png", q_value_barplot, width = 6, height = 6, dpi = 600)

## gene count
gene_count_barplot <- ggplot(summary_databases,
                          aes(x = coviddi, y= avg_gene_count, fill = database)) +
  geom_bar(stat = 'identity', width = 0.5, position = position_dodge(0.6)) +
  geom_errorbar(aes(ymin = avg_gene_count - se_gene_count,
                    ymax = avg_gene_count + se_gene_count),
                width = 0.15,
                position = position_dodge(0.6)) +
  labs(#x = "\nPathways Related to COVID-19",
       x=" ", y = "Average Gene Count\n",
       fill = "Database") +
  theme_classic() +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right"
        ) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    #axis.title.x = element_text(size = 14, face = "bold"), 
    axis.title.y = element_text(size = 14, face = "bold")
  ) 

print(gene_count_barplot)
ggsave("gene_count_barplot.png", gene_count_barplot, width = 6, height = 6, dpi = 600)


## combined plot
grouped_plot <- q_value_barplot + gene_count_barplot +
  plot_layout(guides = "collect") +
  theme(plot.margin = unit(c(1, 0, 0, 1), "cm")) + 
  labs(x ="\nPathways Related to COVID-19") +
  theme(
    axis.title.x = element_text(size = 14, face = "bold", hjust = -3.5)
  )

print(grouped_plot)
ggsave("grouped_barplot.png", grouped_plot, width = 10, height = 5, dpi = 600)


# Statistical test --------------------------------------------------------

library(coin)

# Reactome

wilcox_test_qval <- wilcox.test(q_value ~ coviddi, data = final_results_dfREAC)
wilcox_test_gc <- wilcox.test(geneCOUNT ~ coviddi, data = final_results_dfREAC)

effect_size_qval_REAC <- final_results_dfREAC |> wilcox_effsize(q_value ~ coviddi)
effect_size_gc_REAC <- final_results_dfREAC |> wilcox_effsize(geneCOUNT ~ coviddi)

effect_size_qval_REAC$effsize
effect_size_qval_REAC$magnitude

effect_size_gc_REAC$effsize
effect_size_gc_REAC$magnitude

# KEGG

wilcox_test_qval <- wilcox.test(q_value ~ coviddi, data = final_results_dfKEGG)
wilcox_test_gc <- wilcox.test(geneCOUNT ~ coviddi, data = final_results_dfKEGG)

effect_size_qval_KEGG <- final_results_dfKEGG |> wilcox_effsize(q_value ~ coviddi)
effect_size_gc_KEGG <- final_results_dfKEGG |> wilcox_effsize(geneCOUNT ~ coviddi)

effect_size_qval_KEGG$effsize
effect_size_qval_KEGG$magnitude

effect_size_gc_KEGG$effsize
effect_size_gc_KEGG$magnitude
# WP

wilcox_test_qval <- wilcox.test(q_value ~ coviddi, data = final_results_dfWP)
wilcox_test_gc <- wilcox.test(geneCOUNT ~ coviddi, data = final_results_dfWP)

effect_size_qval_WP <- final_results_dfWP |> wilcox_effsize(q_value ~ coviddi)
effect_size_gc_WP <- final_results_dfWP |> wilcox_effsize(geneCOUNT ~ coviddi)

effect_size_qval_WP$effsize
effect_size_qval_WP$magnitude

effect_size_gc_WP$effsize
effect_size_gc_WP$magnitude

