# script for pca & umap after varstab for each accession separately
library("tidyverse")
library("cowplot")
library("FactoMineR")
library("factoextra")
library("RColorBrewer")
library("umap")
library("patchwork")
set.seed(42)
# import filtered count tables and metatables
getwd()
fnames <- list.files(pattern = "*_vsdTPM.csv*")
cntlist = list()
for (i in fnames) {
    counts <- read.table(i, row.names = 1, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    cnt <- counts %>% mutate_if(is.numeric, as.integer)
    cntlist[[i]] <- cnt
}
getwd()
wanted  <- list.files(pattern = ("_meta.csv"))
unwanted <- list.files(pattern = ("PanBaRT20_geneTPM_meta.csv"))
fnames <- base::setdiff(wanted, unwanted)
#str(fnames)
metalist = list()
for (i in fnames) {
    m <- read.table(i, row.names = 1, header = TRUE, sep = ",", stringsAsFactors = TRUE)
    #cnt <- cnt %>% mutate_if(is.numeric, as.integer)
    metalist[[i]] <- m
}
#
getwd()
#
cols <- c("ID", "Accession", "Tissue", "Batch")
#
for (i in 1:20) {
    #if (metalist[[i]]$Accession[1] == "Hockett") next
    if (all(colnames(cntlist[[i]]) == rownames(metalist[[i]]))) {
        #cnt_cols <- colnames(cntlist[[i]])
        cnt <- as.data.frame(t(cntlist[[i]]), 
                             stringsAsFactors = TRUE) %>% rownames_to_column("ID")
        #cnt$ID <- cnt_cols
        #print(head(cnt, 2))
        meta <- metalist[[i]] %>% rownames_to_column("ID")
        cnt_pca <- inner_join(cnt, meta, by = "ID") %>%
                            mutate(id = ID) %>% 
                            column_to_rownames("ID") %>% 
                            rename(ID = id) %>% 
                            select(ID, Accession, Tissue, Batch, starts_with("chr"))
        # convert to factors
        cnt_pca[cols] <- lapply(cnt_pca[cols], factor)
        # drop extra factor levels
        cnt_pca <- droplevels(cnt_pca)
        # run pca
        print(dim(cnt_pca))
#####
# PCA
#####
        pca <- PCA(X = cnt_pca, scale.unit = FALSE, ncp = 4, quali.sup = c(1:4), graph = F)
        # visualize
        scree <- fviz_eig(pca, addlabels = TRUE)
        pca1 <- fviz_pca_ind(X = pca, 
                             label = "none", 
                             legend.title = "Tissue", 
                             title = paste0(meta$Accession[1], " variance stabilized TPM counts"), 
                             addEllipses = FALSE, axes = c(1, 2), habillage = cnt_pca$Tissue, pointsize = 6) +
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
        pca2 <- fviz_pca_ind(X = pca, 
                             label = "none", 
                             legend.title = "Tissue", 
                             title = paste0(meta$Accession[1], " variance stabilized TPM counts"), 
                             addEllipses = FALSE, axes = c(3, 4), habillage = cnt_pca$Tissue, pointsize = 6) +
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
######
# UMAP
######
        umap = cnt_pca[, grep("chr", colnames(cnt_pca))]
        umap_labels = cnt_pca[1:8]
        # create umap table
        u <- umap(umap, n_neighbors=13)
        print(head(u))
        # inspect table
        print(head(u$layout, 3))
        udf <- as.data.frame(u$layout)
        print(glimpse(udf))
        # visualize
        u1 <- ggplot(udf, aes(x = V1, y = V2, color = umap_labels$Tissue, shape = umap_labels$Tissue)) + 
                labs(title= paste0(meta$Accession[1], " variance stabilized TPM counts"), 
                     x ="UMAP1", y = "UMAP2") +
                geom_point(size = 5) +
                scale_color_manual(name = "Tissue", 
                                   values = c("#D91E36", "#611C35", "#2E5077", "#48A9A6", "#C4A69D"), 
                                   limits = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root"),
                                   labels = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root")) +
                scale_shape_manual(name = "Tissue",
                                   limits = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root"),
                                   labels = c("Caryopsis", "Inflorescence", "Coleoptiles", "Shoot", "Root"),
                                   values = c(15,19,17,18,20)) +
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
        comp <- (u1 + scree) / (pca1 + pca2)
        ggsave(paste0(metalist[[i]]$Accession[1], "_umap_pca_varstab.png"), 
               plot = comp, device = "png", 
               width = 32, 
               height = 20)
        }
}
sessionInfo()
