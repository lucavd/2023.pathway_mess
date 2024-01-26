# Query execution
search_topic <- paste0('("', 'Pulmonary Arterial Hypertension', '"[Mesh]) AND ("2018/01/01"[Date - Entrez] : "2018/12/31"[Date - Entrez])')
search_query <- EUtilsSummary(search_topic, retmax = 1000)

# Gene retrieval
filtered_data <- data.frame('PMID' = (search_query@PMID))
common_als <- intersect(pubtator$PMID, as.numeric(filtered_data$PMID))
annot_als <- pubtator[pubtator$PMID %in% common_als,]

# Pathway analysis
GO_res <- enrichGO(gene = annot_als$Gene, OrgDb = org.Hs.eg.db , 
                   keyType = "ENTREZID", ont = "ALL", 
                   pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.01,
                   readable = TRUE, pool = TRUE)

path <- dplyr::tibble('Pathway' = GO_res@result[["Description"]],
                      'q_value' = GO_res@result[["qvalue"]]) %>% 
  dplyr::filter(q_value < 0.01) %>% 
  rowid_to_column(var = "rowid")

path$coviddi <- stringr::str_detect(path$Pathway, 
                                    stringr::regex("COVID|coronavirus|SARS-CoV-2|COVID-19", 
                                                   ignore_case = TRUE))

somma <- sum(path$coviddi, na.rm = TRUE)
matched_pathways <- path$Pathway[path$coviddi]
pathway_sig <- path[path$Pathway %in% matched_pathways, ]
