# jaccard index ---------------------------------------------------

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

#### common genes ######
pubtator <- read.delim(file = here::here("gene2pubtator"),quote = "", header = TRUE,
                       col.names = c('PMID', 'Object', 'Gene', 'Gene_name', 'Dataset'))


common_als <- (intersect(pubtator$PMID, as.numeric(pubmed_data_als) )  )
annot_als <- pubtator[pubtator$PMID %in% common_als,]
annot_als$TYPE <- "ALS"

common_cad <- (intersect(pubtator$PMID, as.numeric(pubmed_data_cad) )  )
annot_cad <- pubtator[pubtator$PMID %in% common_cad,]
annot_cad$TYPE <- "CAD"

common_cardiomio <- (intersect(pubtator$PMID, as.numeric(pubmed_data_cardiomio) )  )
annot_cardiomio <- pubtator[pubtator$PMID %in% common_cardiomio,]
annot_cardiomio$TYPE <- "Cardiomyopathy"

common_arr <- (intersect(pubtator$PMID, as.numeric(pubmed_data_arr) )  )
annot_arr <- pubtator[pubtator$PMID %in% common_arr,]
annot_arr$TYPE <- "Arrhythmia"

common_hf <- (intersect(pubtator$PMID, as.numeric(pubmed_data_hf) )  )
annot_hf <- pubtator[pubtator$PMID %in% common_hf,]
annot_hf$TYPE <- "Heart Failure"

common_rf <- (intersect(pubtator$PMID, as.numeric(pubmed_data_rf) )  )
annot_rf <- pubtator[pubtator$PMID %in% common_rf,]
annot_rf$TYPE <- "Respiratory Failure"

common_pneumo <- (intersect(pubtator$PMID, as.numeric(pubmed_data_pneumo) )  )
annot_pneumo <- pubtator[pubtator$PMID %in% common_pneumo,]
annot_pneumo$TYPE <- "Pneumonia"

common_met <- (intersect(pubtator$PMID, as.numeric(pubmed_data_met) )  )
annot_met <- pubtator[pubtator$PMID %in% common_met,]
annot_met$TYPE <- "Metabolic Syndrome"

common_dys <- (intersect(pubtator$PMID, as.numeric(pubmed_data_dys) )  )
annot_dys <- pubtator[pubtator$PMID %in% common_dys,]
annot_dys$TYPE <- "Dystrophy"


# All_annot <- rbind(annot_als, 
#                    annot_cad, 
#                    annot_cardiomio, 
#                    annot_arr,
#                    annot_hf,
#                    annot_rf,
#                    annot_pneumo,
#                    annot_met,
#                    annot_dys
#                    )



# jaccard -----------------------------------------------------------------


jaccard(annot_als$Gene, annot_cad$Gene)
1 - jaccard(annot_als$Gene, annot_cad$Gene)

jaccard(annot_als$Gene, annot_cardiomio$Gene)
1 - jaccard(annot_als$Gene, annot_cardiomio$Gene)

jaccard(annot_als$Gene, annot_arr$Gene)
1 - jaccard(annot_als$Gene, annot_arr$Gene)

jaccard(annot_als$Gene, annot_hf$Gene)
1 - jaccard(annot_als$Gene, annot_hf$Gene)

jaccard(annot_als$Gene, annot_rf$Gene)
1 - jaccard(annot_als$Gene, annot_rf$Gene)

jaccard(annot_als$Gene, annot_pneumo$Gene)
1 - jaccard(annot_als$Gene, annot_pneumo$Gene)

jaccard(annot_als$Gene, annot_met$Gene)
1 - jaccard(annot_als$Gene, annot_met$Gene)

jaccard(annot_als$Gene, annot_dys$Gene)
1 - jaccard(annot_als$Gene, annot_dys$Gene)




# diagrammi---------------------------------------------------------

VennDiagram::venn.diagram(list(ALS=annot_als$Gene,
                               CAD=annot_cad$Gene),
                          here::here("ALS/venn_als_cad.tiff"))

Reactome_res.ALS <- enrichPathway(gene= annot_als$Gene,pvalueCutoff= 0.01,organism = "human", pAdjustMethod = "BH", readable=T)
Reactome_res.cad <- enrichPathway(gene= annot_cad$Gene,pvalueCutoff= 0.01,organism = "human", pAdjustMethod = "BH", readable=T)
Reactome_res.cardiomio <- enrichPathway(gene= annot_cardiomio$Gene,pvalueCutoff= 0.01,organism = "human", pAdjustMethod = "BH", readable=T)
Reactome_res.arr <- enrichPathway(gene= annot_arr$Gene,pvalueCutoff= 0.01,organism = "human", pAdjustMethod = "BH", readable=T)
Reactome_res.dys <- enrichPathway(gene= annot_dys$Gene,pvalueCutoff= 0.01,organism = "human", pAdjustMethod = "BH", readable=T)
Reactome_res.hf <- enrichPathway(gene= annot_hf$Gene,pvalueCutoff= 0.01,organism = "human", pAdjustMethod = "BH", readable=T)
Reactome_res.met <- enrichPathway(gene= annot_met$Gene,pvalueCutoff= 0.01,organism = "human", pAdjustMethod = "BH", readable=T)
Reactome_res.pneumo <- enrichPathway(gene= annot_pneumo$Gene,pvalueCutoff= 0.01,organism = "human", pAdjustMethod = "BH", readable=T)
Reactome_res.rf <- enrichPathway(gene= annot_rf$Gene,pvalueCutoff= 0.01,organism = "human", pAdjustMethod = "BH", readable=T)

VennDiagram::venn.diagram(list(ALS=Reactome_res.ALS@result[Reactome_res.ALS@result$p.adjust<0.01,]$ID,
                               CAD = Reactome_res.cad@result[Reactome_res.cad@result$p.adjust<0.01,]$ID),
                          here::here("ALS/venn_als_cad_path.tiff"))  

db_list <- list(Reactome_res.ALS@result,
                Reactome_res.cad@result,
                Reactome_res.cardiomio@result,
                Reactome_res.arr@result,
                Reactome_res.dys@result,
                Reactome_res.hf@result,
                Reactome_res.met@result,
                Reactome_res.pneumo@result,
                Reactome_res.rf@result)

#common_path <- reduce(db_list, inner_join(by = join_by(ID)))

common_path <- plyr::join_all(db_list, by = 'ID', type = 'inner', match = 'first')


common_path[1:9] %>% 
  filter(p.adjust < 0.01) |> 
  filter(Count > 150) %>% 
  ggplot2::ggplot(aes(y = fct_reorder(as.factor(Description), Count),  x = Count, fill = p.adjust)) + 
  geom_col() + xlab("") + ylab("")

ggsave('common_path.png', device = 'png', dpi = 600)
