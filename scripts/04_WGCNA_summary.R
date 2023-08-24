library("tidyverse")
set.seed(42) # For reproducibility of results
# import wgcna results tables
setwd("/lustre/groups/pgsb/projects/panbarley_networks/03_wgcna_results/")
fnames  <- list.files(pattern = ("_ortmodnames.csv"))
str(fnames)
wgcnalist = list()
# vsd table should be : GeneIDs as rows, sample names as cols
for (i in fnames) {
        df <- read.table(i, sep=",", header = TRUE, stringsAsFactors = FALSE)
        # Get the first part of the file name
        new_name <- unlist(strsplit(i, "_"))[1]
        wgcnalist[[new_name]] <- df
}
str(wgcnalist[1])
wgcna_all <- do.call(rbind, wgcnalist)
rownames(wgcna_all) <- NULL
dim(wgcna_all)
