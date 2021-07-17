# Instillation and loading of packages
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager") 
BiocManager::install(version = "3.9") 
BiocManager::install("edgeR")
BiocManager::install("preprocessCore")

install.packages("magrittr")
install.packages("dplyr")
install.packages("rebus")
install.packages("gplots")
install.packages("readr")
install.packages("RColorBrewer")
install.packages("GlobalEnv")
install.packages("stringi")
install.packages("statmod")
install.packages("gridExtra")
install.packages("tidyr")

require(edgeR)
require(magrittr)
require(dplyr)
require(tibble)
require(preprocessCore)
require(rebus)
require(gplots)
require(readr)
require(RColorBrewer)
require(stringi)
require(statmod)
require(grid)
require(gridExtra)
require(tidyr)


#####################################################
# 1. Loading and prepare necessary data to analysis #
#####################################################
# create sample info
samples <- data.frame(read.table('/home/rstudio/data/analysis/samples.table',header = TRUE, sep = "\t"))[,c(1,2)]
samples$file <- substr(samples$file,8,9)

groups <- data.frame(read.table('/home/rstudio/data/analysis/groups.csv',header = TRUE, sep = ","))
samples$group <- groups$group[match(samples$file, groups$ID)]
samples$treatment <- groups$treatment[match(samples$file, groups$ID)]

fpkm <- data.frame(read.table('/home/rstudio/data/genes.fpkm_table',header = TRUE, sep = "\t"))

mm.genes <- data.frame(read.table('/home/rstudio/data/analysis/mart_export.txt', header = TRUE, sep = "\t"))
rownames(fpkm) <- fpkm$tracking_id

fpkm$gene.name <- mm.genes$Gene.name[match(fpkm$tracking_id, mm.genes$Gene.stable.ID)]

# create object to edgeR
fpkm.toEdgeR <- fpkm %>% .[, -1] %>% .[, -25]
group <- samples %>% mutate(gt = paste(group, treatment, sep = "_")) %>% .[, 5] %>% factor()

# normalized fpkm
fpkms.normalised <- data.matrix(fpkm[,c(-1, -26)])
fpkms.normalised <- normalize.quantiles(fpkms.normalised,copy=FALSE)
fpkms.log <- log2(fpkms.normalised + 1)

rm(fpkms.normalised)

fpkms.log[rowMeans(fpkms.log) < 1,] <- NA

############################
# Prepare data for heatmap #
############################
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

heatmap.genes <- function(x, title = "", sizeRow = 0.8, lhei = c(1,5)) {
  fpkms.log %>% 
    as.data.frame() %>% 
    rownames_to_column(., var = "esemblid") %>% 
    filter(esemblid %in% x$esemblid) %>%
    column_to_rownames(., var = "esemblid") %>%
    apply(1, scale) %>%
    t %>%
    apply(1, cut.threshold, threshold = 3) %>%
    t %>%
    heatmap.2(
      distfun = function(x) as.dist(1-cor(t(x))),
      col=rev(morecols(50)),trace="none",
      Colv = FALSE,
      main=title,
      scale="row",
      sepwidth = c(0.3,0.3),
      labRow=fpkm$gene.name[match(rownames(.), rownames(fpkm))],
      srtCol = 0,
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

# Function to extract gene.name
extract.gene.name <- function(df){
  df %>%
  .[, 1] %>% 
    as.data.frame() %>% 
    left_join(., fpkm[, c(1, 26)], by = c("." = "tracking_id")) %>% 
    .[, 2] %>% 
    t %>% 
    t %>% 
    as.data.frame()
}

#####################
# analysis in edgeR #
#####################
########################################
# 2. prepare data to analysis in edgeR #
########################################

# Creating DGEList object
y <- DGEList(counts = fpkm.toEdgeR, group=group)

# Filter and normalization
keep <- filterByExpr(y)
summary(keep)

y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples

# Estimating dispersion
y <- estimateDisp(y)


###########################################
# 3 Testing and choosing gene by contrast #
###########################################

design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- glmQLFit(y, design)


my.contrasts <- makeContrasts(
  ctr.lapvswater = ctrl_lap - ctrl_water,
  eae.lapvswater = eae_lap - eae_water,
  water.ctrvseae =  ctrl_water - eae_water,
  lap.ctrvseae = ctrl_lap - eae_lap,
  group.effect = ((ctrl_water + ctrl_lap)/2) - ((eae_water + eae_lap)/2),
  treatment.effect = (ctrl_water + eae_water) - (ctrl_lap + eae_lap),
  interaction = (ctrl_lap - ctrl_water) - (eae_lap - eae_water),
  levels = design
)


# interaction
png("/home/rstudio/data/analysis/heatmap-edgeR-interaction.png",
    width = 4960,
    height = 7016,
    res = 600)

glmQLFTest(fit, contrast = my.contrasts[, "interaction"]) %>%
  topTags(sort.by = "PValue", n = 1000) %>%  
  as.data.frame() %>% 
  rownames_to_column(., var = "esemblid") %>%
  filter(PValue < 0.05) %>% 
  heatmap.genes(., title = "interaction edger, pvalue = 0.05, n = 77",
                sizeRow = 1.1) 

dev.off()

# # List of genes from interaction
# glmQLFTest(fit, contrast = my.contrasts[, "interaction"]) %>%
#   topTags(sort.by = "PValue", n = 1000) %>%
#   as.data.frame() %>%
#   rownames_to_column(., var = "esemblid") %>%
#   filter(PValue < 0.05) %>% 
#   extract.gene.name(.) %>%
#   write.table(., file = "/home/rstudio/data/analysis/edger-interaction-genename.tsv", sep = "\t", col.names = FALSE,
#               row.names = FALSE, quote = FALSE)
  

# group effect
png("/home/rstudio/data/analysis/heatmap-edgeR-group-effect.png",
    width = 4960,
    height = 7016,
    res = 600)

glmQLFTest(fit, contrast = my.contrasts[, "group.effect"]) %>%
  topTags(sort.by = "PValue", n = 10000) %>%  
  as.data.frame() %>% 
  rownames_to_column(., var = "esemblid") %>%
  filter(PValue < 0.01) %>%
  heatmap.genes(., title = "group edger, pvalue = 0.01, n = 143", sizeRow = 0.6)

dev.off()

# # List of genes from group.effect
# glmQLFTest(fit, contrast = my.contrasts[, "group.effect"]) %>%
#   topTags(sort.by = "PValue", n = 10000) %>%  
#   as.data.frame() %>% 
#   rownames_to_column(., var = "esemblid") %>%
#   filter(PValue < 0.01) %>% 
#   extract.gene.name(.) %>%
#   write.table(., file = "/home/rstudio/data/analysis/edger-group-genename.tsv", 
#               sep = "\t", 
#               col.names = FALSE,
#               row.names = FALSE, quote = FALSE)

# Treatment effect
png("/home/rstudio/data/analysis/heatmap-edgeR-treatment-effect.png",
    width = 4960,
    height = 7016,
    res = 600)

glmQLFTest(fit, contrast = my.contrasts[, "treatment.effect"]) %>%
  topTags(sort.by = "PValue", n = 10000) %>%  
  as.data.frame() %>% 
  rownames_to_column(., var = "esemblid") %>%
  filter(PValue < 0.01) %>%
  heatmap.genes(., title = "treatment edger, pvalue = 0.01, n = 105", sizeRow = 0.8)

dev.off()

# # List of genes from treatment.effect
# glmQLFTest(fit, contrast = my.contrasts[, "treatment.effect"]) %>%
#   topTags(sort.by = "PValue", n = 10000) %>%  
#   as.data.frame() %>% 
#   rownames_to_column(., var = "esemblid") %>%
#   filter(PValue < 0.01) %>% 
#   extract.gene.name(.) %>%
#   write.table(., file = "/home/rstudio/data/analysis/edger-treatment-genename.tsv", 
#               sep = "\t", 
#               col.names = FALSE,
#               row.names = FALSE, quote = FALSE)
