# script for deseq2 varstab
library("tidyverse")
library("DESeq2")
set.seed(42)
# import filtered count tables and metatables
getwd()
fnames <- list.files(pattern = "PanBaRT20_geneTPM_ort_filt_*")
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
    m <- read.table(i, row.names = 1, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    #cnt <- cnt %>% mutate_if(is.numeric, as.integer)
    metalist[[i]] <- m
}
#
getwd()
# 
for (i in 1:20) {
    if (all(colnames(cntlist[[i]]) == rownames(metalist[[i]]))) {
        # create design
        dds <- DESeqDataSetFromMatrix(countData = cntlist[[i]],
                                   colData = metalist[[i]],
                                   design = ~ Tissue)
        dds <- estimateSizeFactors(dds)
        print(summary(sizeFactors(dds)))
        dds <- estimateDispersions(dds)
        vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
        vsdMat <- assay(vsd)
        print(dim(vsdMat))
        #print(head(cntlist[[i]], 1))
        #print(head(metalist[[i]], 3))
        write.table(vsdMat, 
                    file = paste0(metalist[[i]]$Accession[1], "_vsdTPM.csv"), 
                    append = FALSE, 
                    quote = FALSE, 
                    sep = ",", 
                    dec = ".", 
                    row.names = TRUE, 
                    col.names = TRUE)
        }
}
sessionInfo()