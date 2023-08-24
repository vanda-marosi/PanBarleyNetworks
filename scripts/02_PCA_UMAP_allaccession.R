# script for pca & umap after varstab for each accession separately
library("tidyverse")
library("cowplot")
library("FactoMineR")
library("factoextra")
library("RColorBrewer")
library("umap")
library("patchwork")
set.seed(42)
# import filtered count tables and joined metatable
setwd("/home/vanda/Documents/PanBarley_transcriptome/01_vsdTPMs_peraccession/")
fnames <- list.files(pattern = "*_vsdTPM.csv")
vsdlist = list()
collist = list()
# first create an overlap of common hogs across all accessions
for (i in fnames) {
    vsd_df <- read.table(i, row.names = 1, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    vsd_t <- as.data.frame(t(vsd_df), stringsAsFactors = TRUE) %>% rownames_to_column("ID")
    cols <- colnames(vsd_t)
    collist[[i]] <- cols
}
str(collist)
common_hogs <- Reduce(intersect, collist)
str(common_hogs)
# there are 11,816 commonly highly expressed hogs
for (i in fnames) {
    vsd_df <- read.table(i, row.names = 1, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    vsd_t <- as.data.frame(t(vsd_df), stringsAsFactors = TRUE) %>% 
                rownames_to_column("ID") %>%
                select(matches(common_hogs))
    vsdlist[[i]] <- vsd_t
}
all_vsd <- do.call(rbind, vsdlist)
dim(all_vsd)
setwd("/home/vanda/Documents/PanBarley_transcriptome/00_meta_peraccession/")
meta <- read.table("PanBaRT20_geneTPM_meta.csv", sep = ",", header = TRUE, stringsAsFactors = TRUE)
dim(meta)
#####
# PCA
#####
cnt_pca <- inner_join(all_vsd, meta, by = "ID") %>%
                mutate(id = ID) %>% 
                column_to_rownames("ID") %>% 
                rename(ID = id) %>% 
                select(ID, Accession, Tissue, Batch, starts_with("chr"))
# convert to factors
cols <- c("ID", "Accession", "Tissue", "Batch")
cnt_pca[cols] <- lapply(cnt_pca[cols], factor)
# drop extra factor levels
cnt_pca <- droplevels(cnt_pca)
# run pca
dim(cnt_pca)
pca <- PCA(X = cnt_pca, scale.unit = FALSE, ncp = 4, quali.sup = c(1:4), graph = F)
# visualize
scree <- fviz_eig(pca, addlabels = TRUE)
pca1 <- fviz_pca_ind(X = pca, label = "none", legend.title = "Tissue", title = "PanBarley variance stabilized TPM counts", addEllipses = FALSE, axes = c(1, 2), habillage = cnt_pca$Tissue, pointsize = 6) +
                scale_color_manual(values = c("#D91E36", "#611C35", "#2E5077", "#48A9A6", "#C4A69D"),
                                   limits = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root"),
                                   labels = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root")) +
                scale_shape_manual(limits = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root"),
                                   labels = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root"),
                                      values = c(15,19,17,18,20)) +
                scale_fill_manual(values = c("#D91E36", "#611C35", "#2E5077", "#48A9A6", "#C4A69D"),
                                  limits = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root"),
                                  labels = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root")) +
                theme(
                    plot.title = element_text(face = "bold", size = 28),
                    axis.text=element_text(size=25, color = "black"), #change font size of axis text
                    axis.title=element_text(size=32), #change font size of axis titles
                    legend.text=element_text(size=25),
                    legend.title=element_text(size=28),
                    legend.key.size = unit(1.3, 'cm'),
                    #legend.position=c(0.85, 0.7),
                    legend.background = element_rect(fill="white", size=0.5, colour ="transparent"),
                    legend.key = element_blank(),
                    plot.background = element_rect(fill = "white", colour = "transparent"),
                    panel.background = element_rect(fill = "transparent", colour = "grey"),
                    panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
                    panel.grid.major = element_line(colour = "lightgrey", linetype = "solid", size = 0.75,),
                    strip.background = element_rect(fill = "transparent", colour = "transparent"),
                    axis.ticks=element_line(colour="black"),
                    axis.ticks.length = unit(3, "pt")) + 
                guides(color = guide_legend(override.aes = list(size = 6))) +
                xlab("PC1") + 
                ylab("PC2")
pca2 <- fviz_pca_ind(X = pca, label = "none", legend.title = "Tissue", title = "PanBarley variance stabilized TPM counts", addEllipses = FALSE, axes = c(3, 4), habillage = cnt_pca$Tissue, pointsize = 6) +
                scale_color_manual(values = c("#D91E36", "#611C35", "#2E5077", "#48A9A6", "#C4A69D"),
                                   limits = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root"),
                                   labels = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root")) +
                scale_shape_manual(limits = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root"),
                                   labels = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root"),
                                      values = c(15,19,17,18,20)) +
                scale_fill_manual(values = c("#D91E36", "#611C35", "#2E5077", "#48A9A6", "#C4A69D"),
                                  limits = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root"),
                                  labels = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root")) +
                theme(
                    plot.title = element_text(face = "bold", size = 28),
                    axis.text=element_text(size=25, color = "black"), #change font size of axis text
                    axis.title=element_text(size=32), #change font size of axis titles
                    legend.text=element_text(size=25),
                    legend.title=element_text(size=28),
                    legend.key.size = unit(1.3, 'cm'),
                    #legend.position=c(0.85, 0.7),
                    legend.background = element_rect(fill="white", size=0.5, colour ="transparent"),
                    legend.key = element_blank(),
                    plot.background = element_rect(fill = "white", colour = "transparent"),
                    panel.background = element_rect(fill = "transparent", colour = "grey"),
                    panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
                    panel.grid.major = element_line(colour = "lightgrey", linetype = "solid", size = 0.75,),
                    strip.background = element_rect(fill = "transparent", colour = "transparent"),
                    axis.ticks=element_line(colour="black"),
                    axis.ticks.length = unit(3, "pt")) + 
                guides(color = guide_legend(override.aes = list(size = 6))) +
                xlab("PC3") + 
                ylab("PC4")
pca3 <- fviz_pca_ind(X = pca, label = "none", legend.title = "Accession", title = "PanBarley variance stabilized TPM counts", addEllipses = FALSE, axes = c(1, 2), habillage = cnt_pca$Accession, pointsize = 6) +
                scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,23,25)) +
                theme(
                    plot.title = element_text(face = "bold", size = 28),
                    axis.text=element_text(size=25, color = "black"), #change font size of axis text
                    axis.title=element_text(size=32), #change font size of axis titles
                    legend.text=element_text(size=25),
                    legend.title=element_text(size=28),
                    legend.key.size = unit(1.3, 'cm'),
                    #legend.position=c(0.85, 0.7),
                    legend.background = element_rect(fill="white", size=0.5, colour ="transparent"),
                    legend.key = element_blank(),
                    plot.background = element_rect(fill = "white", colour = "transparent"),
                    panel.background = element_rect(fill = "transparent", colour = "grey"),
                    panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
                    panel.grid.major = element_line(colour = "lightgrey", linetype = "solid", size = 0.75,),
                    strip.background = element_rect(fill = "transparent", colour = "transparent"),
                    axis.ticks=element_line(colour="black"),
                    axis.ticks.length = unit(3, "pt")) + 
                guides(color = guide_legend(override.aes = list(size = 6))) +
                xlab("PC1") + 
                ylab("PC2")
######
# UMAP
######
umap = cnt_pca[, grep("chr", colnames(cnt_pca))]
umap_labels = cnt_pca[1:8]
# create umap table
u <- umap(umap)
head(u)
# inspect table
head(u$layout, 3)
udf <- as.data.frame(u$layout)
glimpse(udf)
# visualize
u1 <- ggplot(udf, aes(x = V1, y = V2, color = umap_labels$Tissue, shape = umap_labels$Accession)) + 
                labs(title= "PanBarley variance stabilized TPM counts", x ="UMAP1", y = "UMAP2") +
                geom_point(size = 5) +
                scale_color_manual(name = "Tissue", values = c("#D91E36", "#611C35", "#2E5077", "#48A9A6", "#C4A69D"), 
                                    limits = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root"),
                                     labels = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root")) +
                scale_shape_manual(name = "Accession",
                                   values = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,23,25)) +
                theme(
                    plot.title = element_text(face = "bold", size = 28),
                    axis.text=element_text(size=25, color = "black"), #change font size of axis text
                    axis.title=element_text(size=32), #change font size of axis titles
                    legend.text=element_text(size=25),
                    legend.title=element_text(size=28),
                    legend.key.size = unit(1.3, 'cm'),
                    legend.background = element_rect(fill="white", size=0.5, colour ="transparent"),
                    legend.key = element_blank(),
                    plot.background = element_rect(fill = "white", colour = "transparent"),
                    panel.background = element_rect(fill = "transparent", colour = "grey"),
                    panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
                    panel.grid.major = element_line(colour = "lightgrey", linetype = "solid", size = 0.75,),
                    strip.background = element_rect(fill = "transparent", colour = "transparent"),
                    axis.ticks=element_line(colour="black"),
                    axis.ticks.length = unit(3, "pt")) + 
                guides(color = guide_legend(override.aes = list(size = 6)))
# combine figures and save as composites
comp <- (pca1 + pca2) / (pca3 + scree)
setwd("/lustre/groups/pgsb/projects/panbarley_networks/02_PCA_UMAP_peraccession/")
ggsave("Panbarley_pca_varstab.png"), plot = comp, device = "png", width = 32, height = 20)
ggsave("Panbarley_umap_varstab.png"), plot = u1, device = "png", width = 8, height = 5)
sessionInfo()
