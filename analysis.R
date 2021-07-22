# Sript written in R 3.4.4

########################
# Loading need package #
########################
require(edgeR)
require(preprocessCore)
require(magrittr)
require(gplots)
require(RColorBrewer)
require(tidyverse)
require(stringi)

require(devtools)
install_github("wjawaid/enrichR")
install.packages("enrichR")
require(enrichR)




#####################################################
# 1. Loading and prepare necessary data to analysis #
#####################################################
samples <- data.frame(read.table('/home/ifpan/projects/ifpan-nalepa-smog/samples.table',header = TRUE, sep = "\t"))[,c(1,2)]
samples$file <- substr(samples$file,8,9)

groups <- data.frame(read.table('/home/ifpan/projects/ifpan-nalepa-smog/analysis/groups.csv',header = TRUE, sep = ","))

samples$group <- groups$group[match(samples$file, groups$ID)]
samples$treatment <- groups$treatment[match(samples$file, groups$ID)]

# write sample info
write.csv(samples, "full-sample-info.csv", row.names = FALSE)

fpkm <- data.frame(read.table('/home/ifpan/projects/ifpan-nalepa-smog/genes.fpkm_table',header = TRUE, sep = "\t"))

mm.genes <- data.frame(read.table('/home/ifpan/projects/ifpan-nalepa-smog/analysis/mart_export.txt', header = TRUE, sep = "\t"))

fpkm$gene.name <- mm.genes$Gene.name[match(fpkm$tracking_id, mm.genes$Gene.stable.ID)]

rownames(fpkm) <- fpkm$tracking_id

fpkms.normalised <- data.matrix(fpkm[,c(-1, -26)])
fpkms.normalised <- normalize.quantiles(fpkms.normalised,copy=FALSE)
fpkms.log <- log2(fpkms.normalised + 1)

rm(fpkms.normalised)

# write fpkmslog
write.csv(fpkms.log, "fpkms-log.csv")

#remove fpkms.log that have rowMeans < 1
fpkms.log[rowMeans(fpkms.log) < 1,] <- NA


###########################################
# 2. Statistical analysis - two-way anova #
###########################################
stat <- function(counts,
                 group,
                 treatment) {
  if (is.na(counts[1])) {
    c(rep(NA, 3))
  } else {
    unlist(summary(aov(counts ~ 
                         group
                       * treatment
                       )))[c(17,18,19)]
  }
}

results <- data.frame(fpkm$gene.name, fpkm$tracking_id, row.names = rownames(fpkm))

apply(fpkms.log,
      1,
      stat,
      group = samples$group,
      treatment = samples$treatment
) %>% t %>%
  data.frame() %>%
  bind_cols(results,.) -> results


colnames(results) <- c('gene.name', 'gene.id', 'group.p', 'treatment.p', 'group.treatment.p')

results$group.p %>% p.adjust(., method = 'fdr') -> results$group.fdr
results$treatment.p %>% p.adjust(., method = 'fdr') -> results$treatment.fdr
results$group.treatment.p %>% p.adjust(., method = 'fdr') -> results$group.treatment.fdr

results <- cbind(results, fpkms.log)

write.csv(results, "two-way-anova-results.csv")


#######################################
# 3. Preparation of object to heatmap #
#######################################
mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)

cut.threshold <- function(x, threshold = 2.5) {
  x[x > threshold] <- threshold
  x[x < -threshold] <- -threshold
  x
}

column.color.marker <- samples %>% 
  mutate(group.teratment = paste(group, treatment, sep = "_")) %>% 
  mutate(color = case_when(group.teratment == "eae_lap" ~ "red", 
                           group.teratment == "eae_water" ~ "green", 
                           group.teratment == "ctrl_lap" ~ "blue", 
                           group.teratment == "ctrl_water" ~ "orange")) %>% 
  .[,6]
  
# Function to create heatmap
heatmap.genes <- function(x, title = "", sizeRow = 0.8, lhei = c(1,5)) {
  x %>% 
    apply(1, scale) %>%
    t %>%
    apply(1, cut.threshold, threshold = 3) %>%
    t %>%
    `colnames<-`(colnames(x)) %>%
    heatmap.2(
      distfun = function(x) as.dist(1-cor(t(x))),
      col=rev(morecols(50)),trace="none",
      Colv = FALSE,
      main=title,
      scale="row",
      sepwidth = c(0.3,0.3),
      labRow=fpkm$gene.name[match(rownames(.), rownames(fpkm))],
      srtCol = 90,
      cexRow = sizeRow,
      offsetCol = 0.1,
      ColSideColors = column.color.marker,
      symkey=FALSE,
      cexCol = 1.1,
      lhei = lhei,
      margins = c(4, 7)
    )
  
  legend("top", title = "Legend",legend=c("ctrl_water","ctrl_lap", "eae_water", "eae_lap"),
         fill=c("orange","blue", "green", "red"), cex=1.0, box.lty=0, ncol = 2)
  
}


#######################################
# 4. Visualisation results on heatmap #
#######################################
# interaction
png("/home/ifpan/projects/ifpan-nalepa-smog/analysis/heatmap-two-way-ANOVA-interaction.png",
    width = 4960,
    height = 7016,
    res = 600)

fpkms.log %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "id") %>% 
  filter(id %in% {filter(results, group.treatment.p < 0.01) %>% .[,2] %>% as.character()}) %>% 
  column_to_rownames(., var = "id") %>% 
  heatmap.genes(., title = "interaction, pvalue = 0.01, n = 105")

dev.off()

# List of genes from interaction
fpkm$gene.name[match(rownames(filter(results, group.treatment.p < 0.01) ), rownames(fpkm))] %>% 
  as.vector() %>% 
  t %>% 
  t %>% 
  as.data.frame() %>% 
  write_tsv(., "/home/ifpan/projects/ifpan-nalepa-smog/analysis/gene-files/anova-interaction-genename.tsv", col_names = FALSE)


# group effect
png("/home/ifpan/projects/ifpan-nalepa-smog/analysis/heatmap-two-way-ANOVA-group.effect.png",
    width = 4960,
    height = 7016,
    res = 600)

fpkms.log %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "id") %>% 
  filter(id %in% {filter(results, group.fdr < 0.01) %>% .[,2] %>% as.character()}) %>% 
  column_to_rownames(., var = "id") %>% 
  heatmap.genes(., title = "group effect, fdr =  0.01, n = 113")

dev.off()

# List of genes from group effect
fpkm$gene.name[match(rownames(filter(results, group.fdr < 0.01) ), rownames(fpkm))] %>% 
  as.vector() %>% 
  t %>% 
  t %>% 
  as.data.frame() %>% 
  write_tsv(., "/home/ifpan/projects/ifpan-nalepa-smog/analysis/gene-files/anova-group-genename.tsv", col_names = FALSE)


# treatment effect from treatment effect
png("/home/ifpan/projects/ifpan-nalepa-smog/analysis/heatmap-two-way-ANOVA-treatment.effect.png",
    width = 4960,
    height = 7016,
    res = 600)

fpkms.log %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "id") %>% 
  filter(id %in% {filter(results, treatment.fdr < 0.01) %>% .[,2] %>% as.character()}) %>% 
  column_to_rownames(., var = "id") %>%
  heatmap.genes(., title = "treatment effect, fdr = 0.01, n = 140", sizeRow = 0.65)

dev.off()

# List of genes from treatment effect
fpkm$gene.name[match(rownames(filter(results, treatment.fdr < 0.01) ), rownames(fpkm))] %>% 
  as.vector() %>% 
  t %>% 
  t %>% 
  as.data.frame() %>% 
  write_tsv(., "/home/ifpan/projects/ifpan-nalepa-smog/analysis/gene-files/anova-treatment-genename.tsv", col_names = FALSE)




###################
# second-stage ANOVA #
###################
fpkms.log %>% 
  as.data.frame %>% 
  rownames_to_column(., var = "tracking_id") %>% 
  filter(tracking_id %in% fpkm$tracking_id[match(rownames(filter(results, group.fdr < 0.05) ), rownames(fpkm))]) %>% 
  column_to_rownames(., var = "tracking_id") %>% 
  as.matrix() -> fpkms.log.twostage

filter(results, group.fdr < 0.05) %>% .[, 1:5] %>% 
  mutate(treatment.fdr = p.adjust(treatment.p, method = 'fdr')) -> results.twostage


# treatment effect
png("/home/ifpan/projects/ifpan-nalepa-smog/analysis/heatmap-second-stage-ANOVA-treatment.effect.png",
    width = 4960,
    height = 3508,
    res = 600)

fpkms.log.twostage %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "id") %>% 
  filter(id %in% {filter(results.twostage, treatment.p < 0.05) %>% .[,2] %>% as.character()}) %>% 
  column_to_rownames(., var = "id") %>% 
  heatmap.genes(., title = "treatment, pvalue = 0.01, n = 27", sizeRow = 1.1,  lhei = c(3,5))

dev.off()

filter(results.twostage, treatment.p < 0.05)  %>% 
  select(gene.name) %>%
  write_tsv(., "/home/ifpan/projects/ifpan-nalepa-smog/analysis/gene-files/anova-treatment-secondstage.tsv", col_names = FALSE)



##############################
# Download data from enrichR #
##############################
dbs <-  c("WikiPathway_2021_Human",
          "KEGG_2021_Human", 
          "GO_Biological_Process_2021",
          "GO_Molecular_Function_2018",
          "Descartes_Cell_Types_and_Tissue_2021")

# enrichR; ANOVA
# interaction
enrichr({fpkm$gene.name[match(rownames(filter(results, group.treatment.p < 0.01) ), rownames(fpkm))] %>% 
    as.vector()}, dbs) %>% 
  bind_rows(., .id = "database.name") %>% 
  select(Term, Genes, P.value, Adjusted.P.value) %>%
  write_tsv(., "/home/ifpan/projects/ifpan-nalepa-smog/analysis/enrichR-gene-interaction.tsv")
  

# group effect
enrichr({fpkm$gene.name[match(rownames(filter(results, group.fdr < 0.01) ), rownames(fpkm))] %>% 
    as.vector()}, dbs) %>%
  bind_rows(., .id = "database.name") %>% 
  select(Term, Genes, P.value, Adjusted.P.value) %>%
  write_tsv(., "/home/ifpan/projects/ifpan-nalepa-smog/analysis/enrichR-gene-group-effect.tsv")


# treatment
enrichr({fpkm$gene.name[match(rownames(filter(results, treatment.fdr < 0.01) ), rownames(fpkm))] %>% 
    as.vector()}, dbs) %>%
  bind_rows(., .id = "database.name") %>% 
  select(Term, Genes, P.value, Adjusted.P.value) %>%
  write_tsv(., "/home/ifpan/projects/ifpan-nalepa-smog/analysis/enrichR-gene-treatment-effect.tsv")

# enrichR; secondstage
enrichr({filter(results.twostage, treatment.p < 0.05)  %>% .[, 1] %>% as.character()}, dbs) %>%
  bind_rows(., .id = "database.name") %>% 
  select(Term, Genes, P.value, Adjusted.P.value) %>%
  write_tsv(., "/home/ifpan/projects/ifpan-nalepa-smog/analysis/enrichR-gene-group-effect-secondstage.tsv")


# # enrichR; ANOVA; gene.markers
# # interaction
# enrichr({results.markers %>% 
#     filter(group.treatment.p < 0.05) %>%
#     .[, 1] %>% as.vector()}, dbs)
# 
# # group effect
# enrichr({results.markers %>% 
#     filter(group.fdr < 0.05) %>%
#     .[, 1] %>% as.vector()}, dbs)
# 
# # treatment
# enrichr({results.markers %>% 
#     filter(treatment.p < 0.05) %>%
#     .[, 1] %>% as.vector()}, dbs)


############################
# 5. Analysis marker genes #
############################
marker.genes <- read.table('/home/ifpan/projects/ifpan-nalepa-smog/analysis/gene-markers.tsv', header = TRUE, sep = "\t") %>%
  rename(gene.name = location) %>%
  rename(location = marker)


fpkm %>% 
  filter(gene.name %in% {
    marker.genes %>% 
      filter(region %in% c("Cortex", "Hippocampus,Cortex", "CNS")) %>% 
      filter(!str_detect(location, "possible batch effect")) %>% .[, 5]}) -> fpkm.markers


fpkms.normalised.markers <- data.matrix(fpkm.markers[, c(-1, -26)])
fpkms.normalised.markers <- normalize.quantiles(fpkms.normalised.markers, copy = FALSE)
fpkms.log.markers <- log2(fpkms.normalised.markers + 1)

rm(fpkms.normalised.markers)


results.markers <- data.frame(fpkm.markers$gene.name, fpkm.markers$tracking_id, row.names = rownames(fpkm.markers))

apply(fpkms.log.markers,
      1,
      stat,
      group = samples$group,
      treatment = samples$treatment
) %>% t %>%
  data.frame() %>%
  bind_cols(results.markers,.) -> results.markers


colnames(results.markers) <- c('gene.name', 'gene.id', 'group.p', 'treatment.p', 'group.treatment.p')

results.markers$group.p %>% p.adjust(., method = 'fdr') -> results.markers$group.fdr
results.markers$treatment.p %>% p.adjust(., method = 'fdr') -> results.markers$treatment.fdr
results.markers$group.treatment.p %>% p.adjust(., method = 'fdr') -> results.markers$group.treatment.fdr

results.markers <- cbind(results.markers, fpkms.log.markers)

# Visualisation results two-way ANOVA for gene.markers
# interaction
png("/home/ifpan/projects/ifpan-nalepa-smog/analysis/heatmap-ANOVA-interaction-gene.markers.png",
    width = 4960,
    height = 3508,
    res = 600)

results.markers %>% 
  filter(group.treatment.p < 0.05) %>% 
  .[, c(9:32)] %>%
  heatmap.genes(., title = "group.treatment, p = 0.05, n = 6", sizeRow = 1.3, lhei = c(3,5))

dev.off()

# marker genes, interaction
marker.genes %>% 
  filter(gene.name %in% {results.markers %>% 
      filter(group.treatment.p < 0.05) %>% 
      .[, c(9:32)] %>% 
      mutate(gene.name = fpkm$gene.name[match(rownames(.), rownames(fpkm))]) %>% 
      .[, 25] %>% 
      as.vector()}) %>%
  select(-location) %>%
  select(gene.name) %>% 
  write_tsv(., "/home/ifpan/projects/ifpan-nalepa-smog/analysis/gene-files/anova-interaction-genemarker.tsv", col_names = FALSE)


# group.effect
png("/home/ifpan/projects/ifpan-nalepa-smog/analysis/heatmap-ANOVA-group.effect-gene.markers.png",
    width = 4960,
    height = 3508,
    res = 600)

results.markers %>% 
  filter(group.fdr < 0.05) %>% 
  .[, c(9:32)] %>% 
  heatmap.genes(., title = "group.effect, fdr = 0.05, n = 6", sizeRow = 1.2, lhei = c(3,5))

dev.off()

# marker genes, group.effect
marker.genes %>% 
  filter(gene.name %in% {results.markers %>% 
      filter(group.fdr < 0.05) %>% 
      .[, c(9:32)] %>% 
      mutate(gene.name = fpkm$gene.name[match(rownames(.), rownames(fpkm))]) %>% 
      .[, 25] %>% 
      as.vector()}) %>%
  select(-location) %>%
  select(gene.name) %>%
  write_tsv(., "/home/ifpan/projects/ifpan-nalepa-smog/analysis/gene-files/anova-group-genemarker.tsv", col_names = FALSE)


# treatment.effect
png("/home/ifpan/projects/ifpan-nalepa-smog/analysis/heatmap-ANOVA-treatment.effect-gene.markers.png",
    width = 4960,
    height = 3508,
    res = 600)

results.markers %>% 
  filter(treatment.p < 0.05) %>% 
  .[, c(9:32)] %>% 
  heatmap.genes(., title = "treatment, p = 0.05, n = 13", sizeRow = 1.2, lhei = c(3,5))

dev.off()

# marker genes treatment.effect
marker.genes %>% 
  filter(gene.name %in% {results.markers %>% 
      filter(treatment.p < 0.05) %>% 
      .[, c(9:32)] %>% 
      mutate(gene.name = fpkm$gene.name[match(rownames(.), rownames(fpkm))]) %>% 
      .[, 25] %>% 
      as.vector()}) %>%
  select(-location) %>%
  select(gene.name) %>%
  write_tsv(., "/home/ifpan/projects/ifpan-nalepa-smog/analysis/gene-files/anova-treatment-genemarker.tsv", col_names = FALSE)



