library(tidyverse)
library(clusterProfiler)
set.seed(42)
#
# import modules2gene table
getwd()
ortmod_names <- list.files(pattern = ("_ortmodnames.csv"))
str(ortmod_names)
datalist_ortmod = list()
for (i in ortmod_names) {
    ortmod_df <- read.table(i, sep=",", header = TRUE, stringsAsFactors = FALSE)
    ortmod_df <- ortmod_df %>%
                    select(ID, module_name) #%>%
                    #mutate(Accession = str_split_n(module_name, "_", 1))
    datalist_ortmod[[i]] <- ortmod_df
}    
mod2gene <- do.call(rbind, datalist_ortmod)
rownames(mod2gene) <- NULL
dim(mod2gene)
head(mod2gene)
# import term2gene tables
# and fuse them into single df per species
getwd()
unwanted  <- list.files(pattern = ("all_"))
wanted <- list.files(pattern = ("_panbarley_mercator.csv"))
mercator_names <- base::setdiff(wanted, unwanted)
str(mercator_names)
datalist_term2gene = list()
for (i in mercator_names) {
        mercator_df <- read.table(i, sep=",", stringsAsFactors = FALSE)
        mercator_df <- mercator_df %>% 
                    mutate(Filename = paste0(i)) %>% 
                    rename(Term = V1, Gene = V2) %>% 
                    mutate(Level = str_extract(Filename, regex("level[1-8]"))) %>% 
                    select(-Filename)
        datalist_term2gene[[i]] <- mercator_df
}
term2gene <- do.call(rbind, datalist_term2gene)
rownames(term2gene) <- NULL
dim(term2gene)
head(term2gene)
##
print("Imports were successful")
# 
getwd()
# create levels
levels = c("level1", "level2", "level3", "level4", "level5", "level6", "level7")
level = as.list(levels)
names(level) = levels
print("Assignments were successful")
#
for(i in unique(mod2gene$module_name)) {
		temp <- mod2gene %>% 
                	filter(module_name == i)
		print(head(temp))
    		genelist <- temp$ID
		#print(head(genelist, 2))
		print("Input was created, enrichments starts")
    		# run enrichment
        for (g in 1:7) {
            print(level[g])
            tempset <- term2gene %>% 
                        filter(Level == paste0(level[g])) %>%
                        select(-Level)                        
            names(tempset) <- NULL
            print(head(tempset, 3))
    		enrich <- enricher(gene = genelist, 
                       	pvalueCutoff = 0.05,
                       	pAdjustMethod = "BH",
                       	#universe = get(paste0("universe_", subg[j])),
                       	minGSSize = 10,
                       	maxGSSize = 500, 
                       	qvalueCutoff = 0.2,
                       	TERM2GENE=tempset)
            print(str(enrich))
    		res <- data.frame(enrich)
            print(head(res))
            if (length(res$ID) >= 1) {
    		write.csv(res, file =paste(level[g], "_module_",i,"_Mercator_clusterprof.csv", sep=""))
    		ggsave(filename = paste(level[g], "_module_",i,"_Mercator_plot.pdf", sep=""), 
			plot = dotplot(enrich),
			width = 15, 
			height = 10)
            }
            else 
            next
        }
}
#
sessionInfo()
