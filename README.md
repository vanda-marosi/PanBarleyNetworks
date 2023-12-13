# PanBarley Co-expression Networks Analysis pipeline

This repository contains the analysis scripts for the pan-barley transcriptome co-expression pipeline. 

### Table of Contents

- [1. Preparation of count & metadata tables](#preparation-of-count---metadata-tables)
- [2. Variance Stabilizing Transformation](#variance-stabilizing-transformation)
- [3. Principal Component Analysis](#principal-component-analysis)
- [4. Weighted Gene Correlation Network Analysis](#weighted-gene-correlation-network-analysis)
- [5. Mercator Functional Enrichment Analysis](#mercator-functional-enrichment-analysis)
- [6. Comparison of Co-expression Networks with Jaccard-distance and Cosine-similarity](#comparison-of-co-expression-networks-with-jaccard-distance-and-cosine-similarity)
- [7. Louvain Method for Community Detection](#louvain-method-for-community-detection)
- [8. OrthoFinder Analysis](#orthofinder-analysis)
- [References](#references)
---

## Preparation of count & metadata tables
R version 4.2.3 (Team, R. D. C. (2010). R: A language and environment for statistical computing. (No Title)) with package set tidyverse (v2.0.0) (Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L. D. A., François, R., ... & Yutani, H. (2019). Welcome to the Tidyverse. Journal of open source software, 4(43), 1686.) revealed oOut of the 13,680 single copy genes, that 13,652 genes have passed the filter of lowly expressing genes (TPM > 0.5 in minimum 2 biological replicates) in all accessions and were taking part of the following analysis.

The [00_DataWrangling.ipynb](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/00_DataWrangling.ipynb) Jupyter Notebook <sup>1</sup> demonstrates how to divide the overall count table (`PanBaRT20_tpm_genes.csv`) from all 297 samples into accession specific sample tables and how to subset the single-copy orthologs from all genes, using the "core;single" keywords on the `GsRTD_and_PanBaRT20_match.tsv` table. The generated `"PanBaRT20_geneTPM_ort_filt_*.csv"` and `"*_meta.csv"` tables serve as input for the next step of the pipeline.

---

## Variance Stabilizing Transformation
R package DESeq2 (v1.34.0) is used to perform variance stabilizing transformation for normalization applying the design “~Tissue” with the argument “blind = FALSE” in the [01_varstabDESeq2.R](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/01_varstabDESeq2.R) script <sup>3</sup>.

---

## Principal Component Analysis 
From the original 297 samples, principal component analysis (PCA) revealed that one sample was a severe outlier with R scripts [02_PCA_UMAP_all.ipynb](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/02_PCA_UMAP_all.ipynb) and [02_PCA_UMAP_peraccession.R](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/02_PCA_UMAP_peraccession.R) for all samples & for each accession separately with FactoMineR (v2.8) (Sebastien Le, Julie Josse, Francois Husson (2008). FactoMineR: An R Package for Multivariate Analysis. Journal of Statistical Software, 25(1), 1-18. 10.18637/jss.v025.i01). Therefore, this sample (ZDM01467_In2) was removed from further analysis. 

---

## Weighted Gene Correlation Network Analysis
Co-expression networks were built for all accessions’ single copy orthologs separately using the WGCNA package (v1.69) (https://doi.org/10.1186/1471-2105-9-559), as follows in [03_WGCNA_peraccession_1to1orthologs.R](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/03_WGCNA_peraccession_1to1orthologs.R). The soft power threshold used for scale-free topology was determined using the automatic detection function “pickSoftThreshold”. Unsigned networks were calculated with Pearson correlation to obtain the adjacency matrix that accounted for both positive and negative correlations among gene-pairs by applying methods dynamicTreeCut and TOMtype with a minimum module size = 30, “unsigned” network type and otherwise default settings. This resulted in twenty networks, with different module numbers for each accession.

This approach was followed by merging closely clustered modules by a cut height 0.1 that resulted in the final set of modules for each accession. Merged modules have been ordered by their size (number of genes) and renamed with accession specific names followed by the order of modules such as “Akashinriki_1” or “B1K_3”. This resulted in 738 modules across all accessions in total. As WGCNA excluded a few genes with insufficient variation from each of the genotype-specific networks using the “goodSamplesGenes” function, 91 genes have been filtered out from the next steps of analysis, leaving 13,561 single copy genes to analyze further. 

Module eigengenes were calculated using the “moduleEigengenes” function of the WGCNA package, which represents the expression profile of a given module by taking the first principal component of the module. Pearson correlation with Fisher’s-exact significance test was calculated among a binarized tissue table and the eigengenes.

---

## Mercator Functional Enrichment Analysis
Mercator4 (v5.0) (PMID: 34448161) functional annotation was performed and enrichment with over-representation analysis (ORA) of sets of module-genes was done using R package clusterProfiler (v4.6) (doi:10.1016/j.xinn.2021.100141) with Benjamini-Hochberg FDR correction p-value cut-off 0.05 as in [04_Mercator_FuncEnrich.R](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/04_Mercator_FuncEnrich.R).

---

## Comparison of Co-expression Networks with Jaccard-distance and Cosine-similarity 
To compare the shared number of orthologs across accession-specific network modules, Jaccard distance and Cosine-similarity was calculated for each module-pairs, and as Cosine-similarity proved to be more sensitive against module size dependent comparisons, it was used in the next step of the analysis as in [05_Network_comparison_Jaccard_Cosine.py](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/05_Network_comparison_Jaccard_Cosine.py) using Python (v3.10.12) (Van Rossum, G., & Drake Jr, F. L. (1995). Python reference manual. Centrum voor Wiskunde en Informatica Amsterdam).

---

## Louvain Method for Community Detection 
Across all accessions, a meta-network was established by computing edges using the Cosine similarity values and nodes as modules using networkx (v3.1) (cite Aric A. Hagberg, Daniel A. Schult and Pieter J. Swart, “Exploring network structure, dynamics, and function using NetworkX”, in Proceedings of the 7th Python in Science Conference (SciPy2008), Gäel Varoquaux, Travis Vaught, and Jarrod Millman (Eds), (Pasadena, CA USA), pp. 11–15, Aug 2008), and clustered further into 6 distinct communities via Louvain-community detection (De Meo, P., Ferrara, E., Fiumara, G., & Provetti, A. (2011, November). Generalized louvain method for community detection in large networks. In 2011 11th international conference on intelligent systems design and applications (pp. 88-93). IEEE.). Final network layout was calculated using the Netgraph (v4.12.11) (cite Brodersen, P. J. N., (2023). Netgraph: Publication-quality Network Visualisations in Python. Journal of Open Source Software, 8(87), 5372, https://doi.org/10.21105/joss.05372) package with edge bundle option as in [06_Networkx_louvain_define_communities.py](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/06_Networkx_louvain_define_communities.py).

---

## OrthoFinder Analysis
Ortholog detection among Oryza sativa (IRGSP-1.0) (PMID:23299411), Arabidopsis thaliana (TAIR10.57) (https://doi.org/10.1111/tpj.13415) and selected barley accessions Morex, Akashinriki and B1K representative proteomes was performed using OrthoFinder (2.5.5) (Emms, D. M., & Kelly, S. (2019). OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome biology, 20, 1-14.).  
---

## References

* 1 Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8.
*
*
*
*

