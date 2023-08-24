# WGCNA
# for Pan-Barley project
library("WGCNA")
library("tidyverse")
options(stringsAsFactors = FALSE)
enableWGCNAThreads(10)
set.seed(42) # For reproducibility of results
# import vsd tables
setwd("/lustre/groups/pgsb/projects/panbarley_networks/01_vsdTPMs_peraccession/")
fnames  <- list.files(pattern = ("_vsdTPM.csv"))
str(fnames)
vsdlist = list()
# vsd table should be : GeneIDs as rows, sample names as cols
for (i in fnames) {
        df_vsd <- read.table(i, sep=",", row.names = 1, header = TRUE, stringsAsFactors = FALSE)
        df_ID <- tibble::rownames_to_column(df_vsd, "ID")
        vsd <- as.data.frame(t(df_ID[-c(0:1)]))
        colnames(vsd) <- df_ID$ID
        print(dim(vsd))
        # Get the first part of the file name
        new_name <- unlist(strsplit(i, "_"))[1]
        vsdlist[[new_name]] <- vsd
}
str(vsdlist[1])
# import metatables
setwd("/lustre/groups/pgsb/projects/panbarley_networks/00_meta_peraccession/")
wanted  <- list.files(pattern = ("_meta.csv"))
unwanted <- list.files(pattern = ("PanBaRT20_geneTPM_meta.csv"))
meta_fnames <- base::setdiff(wanted, unwanted)
str(meta_fnames)
metalist = list()
for (i in meta_fnames) {
        df <- read.table(i, sep=",", header = TRUE, stringsAsFactors = TRUE)
        print(dim(df))
        # Get the first part of the file name
        new_name <- unlist(strsplit(i, "_"))[1]
        metalist[[new_name]] <- df
}
str(metalist[1])

# Store softpower estimates for later
softPowerVector <- vector("numeric", length = length(vsdlist))

# for bubbleplot
mid <- 0

# Loop through vsdlist
setwd("/lustre/groups/pgsb/projects/panbarley_networks/03_wgcna_results/")
for (i in names(vsdlist)) {
    print(i)
    vsd <- vsdlist[[i]] 
    meta <- metalist[[i]]
    # Check for quality
    check <- goodSamplesGenes(vsd, verbose = 3)
    # Remove problematic samples&genes
    if (!check$allOK) {
        # Print the gene and sample names that were removed:
        if (sum(!check$goodGenes) > 0) {
            printFlush(paste("Removing genes:", paste(names(vsd)[!check$goodGenes], collapse = ", ")));
        }
        if (sum(!check$goodSamples) > 0) {
            printFlush(paste("Removing samples:", paste(rownames(vsd)[!check$goodSamples], collapse = ", ")));
        }
        # Remove the offending genes and samples from the data:
        vsd <- vsd[check$goodSamples, check$goodGenes]
    }
    # Choose a set of soft-thresholding powers
    powers <- c(c(1:14), seq(from = 15, to = 21, by = 2))
    sft_w <- pickSoftThreshold(vsd, powerVector = powers, verbose = 5, networkType = "unsigned")

    # Plot the results
    pdf(file = paste(i, "_NetworkTopology.pdf", sep = ""), width = 8, height = 5)
    par(mfrow = c(1, 2))
    cex1 = 0.9
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft_w$fitIndices[,1], -sign(sft_w$fitIndices[,3])*sft_w$fitIndices[,2],
    xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,unsigned R^2",type="n",
    main = paste("Scale independence"))
    text(sft_w$fitIndices[,1], -sign(sft_w$fitIndices[,3])*sft_w$fitIndices[,2],
    labels=powers,cex=cex1,col="red")
    # this line corresponds to using an R^2 cut-off of h
    abline(h=0.87,col="red") # here cuts at 6
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft_w$fitIndices[,1], sft_w$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
    text(sft_w$fitIndices[,1], sft_w$fitIndices[,5], labels=powers, cex=cex1,col="red")
    dev.off()

    # Get the softPower value
    softPowerVector[i] <- sft_w$powerEstimate
    softPower <- sft_w$powerEstimate
    print(softPower)
    # Replacing NA with 9
    if (is.na(sft_w$powerEstimate)) {
        softPowerVector[i] <- 9
        softPower <- 9
    } else {
        softPowerVector[i] <- sft_w$powerEstimate
        softPower <- sft_w$powerEstimate
    }
    
    # Calculate the adjacency matrix
    # this is a metric to determine which genes j and i are connected
    adjacency <- adjacency(vsd,
                        type = "unsigned",
                        power = softPower)
    # calculate the topological overlap matrix (TOM) - more robust than adjancency
    # this estimates how many connections gene j and i have in common
    TOM <- TOMsimilarity(adjacency)
    colnames(TOM) <- colnames(vsd)
    rownames(TOM) <- colnames(vsd) # maybe here rownames?
    # save rds recommended here
    collectGarbage() # to clean R
    
    # Clustering
    geneTree = hclust(as.dist(1-TOM), method = "average");
    # We like large modules, so we set the minimum module size relatively high:
    minModuleSize = 30;
    # Module identification using dynamic tree cut:
    unmergedLabels = cutreeDynamic(dendro = geneTree,
                               distM = 1-TOM,
                               deepSplit = 2,
                               cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE );
    # The  label  0  isreserved for genes not assigned to any of the modules.
    unmergedColors = labels2colors(unmergedLabels)
    print(table(unmergedColors))
    print(table(unmergedLabels))
    sizeGrWindow(8,6);
    pdf(file = paste(i, "_Dendrogram_unmerged.pdf", sep = ""), wi = 8, he = 6)
    plotDendroAndColors(geneTree,
                        unmergedColors, 
                        "Dynamic Tree Cut",
                        dendroLabels = FALSE, 
                        hang = 0.03, 
                        addGuide = TRUE, 
                        guideHang = 0.05)
    dev.off()

    # Calculate Module Eigengenes (ME)
    MElist <- moduleEigengenes(vsd, colors = unmergedColors)
    MEs <- MElist$eigengenes
    # calculate dissimilarity of eigengenes
    MEdiss <- 1-cor(MEs)
    MEtree <- hclust(as.dist(MEdiss), method = "average")
    print(dim(MEs))
    # clustering of TOM Module Eigengenes
    sizeGrWindow(14,8);
    pdf(file = paste(i, "_MEclustering_unmerged.pdf", sep = ""), wi = 8, he = 6)
    par(mfrow = c(1,1))
    plot(MEtree, main = paste0(i, " Clustering of Module Eigengenes"), xlab = "", sub = "")
    abline(h=0.1, col = "red")
    #abline(h=0.05, col = "green")
    dev.off()
    
    # I decided to merge modules for which the eigengenes have a distance smaller than 0.1
    # this roughly corresponds to a correlation of 0.9
    merge <- mergeCloseModules(vsd, unmergedColors, cutHeight = 0.1)
    mergedColors <- merge$colors
    # Eigengenes of the new merged modules:
    mergedMEs = merge$newMEs
    table(mergedColors)
    #
    sizeGrWindow(12, 9)
    pdf(file = paste(i, "_Dendogram_merged.pdf", sep = ""), wi = 8, he = 6)
    plotDendroAndColors(geneTree, 
                        cbind(unmergedColors, mergedColors),
                        c("Dynamic Tree Cut", "Merged at 0.05"),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
    dev.off()
    
    # Rename to moduleColors
    moduleColors = mergedColors
    # Construct numerical labels corresponding to the colors
    colorOrder = c("grey", standardColors(50));
    moduleLabels = match(moduleColors, colorOrder)-1;
    MEs = mergedMEs;
    
    # Create heatmap displaying the relationship of the identifyed modules
    ## using the topological overlap dissimilarity matrix (TOM)
    #### Lighter shades of yellow indicate closer expression neighborhood
    disTOM <- 1-TOM
    sizeGrWindow(8,5);
    pdf(file = paste(i, "_TOM.pdf", sep = ""), wi = 8, he = 5)
    TOMplot(disTOM^7, geneTree, as.character(mergedColors))
    dev.off()
    
    # Relating consensus modules to traits
    # Tissues
    # use wgcna func to binarise, subset tissue & sample cols
    tissue <- meta$Tissue
    sample <- meta$ID
    bin_t = binarizeCategoricalVariable(tissue, # same as binarizeCategoricalColumns()
                                   includePairwise = FALSE,
                                   includeLevelVsAll = TRUE,
                                   minCount = 1,
                                   dropFirstLevelVsAll = FALSE);
    bin_t <- data.frame(sample, bin_t)
    rownames(bin_t) <- bin_t[,1]
    bin_t[,1] <- NULL
    sorted_tissue <- unique(sort(tissue))
    colnames(bin_t) <- sorted_tissue
    print(head(bin_t, 3))
    print(dim(bin_t))
    moduleTraitCor = cor(merge$newMEs, bin_t, use = "p");
    moduleTraitPvalue = corPvalueFisher(moduleTraitCor, 544)
    MEColors = substring(names(merge$newMEs), 3);
    MEColorNames = paste(MEColors, sep="");
    
    # Plot the module-trait relationship table as in WGCNA tutorial
    sizeGrWindow(8,14)
    pdf(file = paste(i, "_ModuleTissueCorrelation_WGCNAtutorialplot.pdf", sep = ""), wi = 8, he = 14);
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)
    par(mar = c(6, 13, 3, 2.2)); #bottom, left, top, and right
    labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(bin_t),
               yLabels = MEColorNames,
               xLabelsAngle = 30,
               ySymbols = MEColorNames,
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5, #0.3,
               cex.lab.x = 1,
               cex.lab.y = 1.5,
               cex.legendLabel = 1,
               legendLabel = "Pearson's correlation",
               zlim = c(-1,1),
               main = paste0(i, " Module Eigengene - tissue relationship"))
    dev.off()

    # Extract module-ortholog lists for comparison plots
    ortlist <- row.names(TOM)
    ort_mod <- as.data.frame(cbind(ortlist, moduleLabels, moduleColors))
    colnames(ort_mod) <- c("ID", "module", "colors")
    print(str(ort_mod))

    # Create cultivar-specific module names
    module_colors <- unique(ort_mod$colors)
    # Based on the presence of grey module, start numbering from 0 or 1
    if ("grey" %in% module_colors) {
      number_of_modules <- as.character(0:(length(module_colors))-1)
    } else {
      number_of_modules <- as.character(1:length(module_colors))
    }
    cultivar_abreviation <- rep(c(paste0(i, "_")), times=length(number_of_modules))
    cult_num <- as.data.frame(paste(cultivar_abreviation, number_of_modules, sep=""))
    print(str(cult_num))
    sorted_ort_mod <- ort_mod %>%
                   group_by(colors) %>%
                   count() %>%
                   rename(size = n) %>%
                   arrange(desc(size)) %>%
                   ungroup()
    new_module_names <- cbind(sorted_ort_mod, cult_num)
    colnames(new_module_names) <- c("colors", "size", "module_name")
    
    # Intersect ID-mod table with new module name table and save
    converted_ort_mod <- inner_join(ort_mod, new_module_names, by ="colors") %>% select(ID, module_name, size, colors)
    write.table(converted_ort_mod, 
                file = paste(i, "_ortmodnames.csv", sep = ""), 
                append = FALSE, quote = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)
    
    # Create better looking bubble plot
    print(str(moduleTraitPvalue))
    fisher <- as.data.frame(moduleTraitPvalue)
    fisher <- fisher %>% 
                rownames_to_column("colors") %>%                     
                mutate(colors = map_chr(colors, ~str_replace(.x, "ME", "")))
    fisher_colors <- fisher$colors
    fisher <- fisher %>% gather(key = "tissue", value = "fisher", Caryopsis:Shoot)
    print(str(fisher))
    pearson <- as.data.frame(moduleTraitCor)
    pearson <- cbind(pearson, fisher_colors) # we can do this because order is the same
    colnames(pearson) <- c("Caryopsis", "Coleoptiles", "Inflorescence", "Root", "Shoot", "colors")
    pearson <- pearson %>% gather(key = "tissue", value = "pearson", Caryopsis:Shoot)
    print(str(pearson))
    pf <- inner_join(pearson, fisher, by = c("colors", "tissue"))
    print(str(pf))
    # intersect with converted names
    pf_conv <- inner_join(pf, converted_ort_mod, by = "colors") %>%
                mutate(Order = as.integer(str_extract(module_name, "_[0-9]+"))) %>%
                arrange(Order)
    print(str(pf_conv))
    
   # Save correlation tables
    write.table(pf_conv, 
                file = paste(i, "_modtraitcor_tissue.csv", sep = ""), 
                append = FALSE, quote = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE)    
    
    # create bubble plot
    p1 <- pf_conv %>%  
        ggplot(aes(x=fct_inorder(module_name), y=tissue, size=fisher, color=pearson)) +
        geom_point() +
        labs(title = paste(i, " module-tissue correlation", sep=""), x = "", y = "") +
        scale_size("Fishers exact p-value",range = c(7, .1)) +
        scale_color_gradient2("Pearson correlation", midpoint = mid, low="blue", mid = "white", high="red",
                         breaks=c(-0.75,-0.25,0,0.5,0.9),labels=c(-0.75,-0.25,0,0.5,0.9)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust =1, size=16, color = "black"),
                axis.text.y=element_text(size=20, color = "black", face="bold"),
                axis.title.y=element_text(size=20, color="black"),
                title=element_text(size=20),
                legend.position="right", 
                legend.key.width=unit(1, "cm"), 
                legend.text=element_text(size=14), 
                legend.title=element_text(size=16),
                plot.background = element_rect(fill = "transparent", colour = "transparent"),
                panel.background = element_rect(fill = "transparent", colour = "grey"),
                panel.border = element_rect(fill = NA, colour = "black", size = 1),
                panel.grid.major = element_line(colour = "lightgrey", linetype = "solid", size = 0.5,),
                strip.background = element_rect(fill = "transparent", colour = "transparent"),
                axis.ticks=element_line(colour="black"),
                axis.ticks.length = unit(8, "pt"))
    # decide plot length based on number of modules
    if (length(module_colors <= 25)) { 
      15 <- width_of_plot
    } else if (length(module_colors <= 30 {
      17 <- width_of_plot
    } else if (length(module_colors <= 40 {
      20 <- width_of_plot
    } else if (length(module_colors <= 50 {
      22 <- width_of_plot
    } else if (length(module_colors <= 60 {
      24 <- width_of_plot
    } else if (length(module_colors <= 70 {
      26 <- width_of_plot
    } else {28 <- width_of_plot}
    ggsave(plot = p1, filename = paste(i, "_ModuleTissueCorrelation.pdf", sep=""), width = width_of_plot, height = 7)
} 

# Create a data frame to store the collected softpower results
softPowerResults <- data.frame(VsdName = names(vsdlist), SoftPower = softPowerVector)
write.csv(softPowerResults, file = "wgcna_soft_power_results.csv", row.names = FALSE)

sessionInfo()










