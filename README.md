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
R version 4.2.3<sup>1</sup> with package set tidyverse (v2.0.0)<sup>2</sup> revealed out of the 13,680 single copy genes, that 13,652 genes have passed the filter of lowly expressing genes (TPM > 0.5 in minimum 2 biological replicates) in all accessions and were taking part of the following analysis.

The [00_DataWrangling.ipynb](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/00_DataWrangling.ipynb) JupyterLab (v3.6.2)<sup>3</sup> notebook demonstrates how to divide the overall count table (`PanBaRT20_tpm_genes.csv`) from all 297 samples into accession specific sample tables and how to subset the single-copy orthologs from all genes, using the `core;single` keywords on the `GsRTD_and_PanBaRT20_match.tsv` table. The generated `PanBaRT20_geneTPM_ort_filt_*.csv` and `*_meta.csv` tables serve as input for the next step of the pipeline.

---

## Variance Stabilizing Transformation
R package DESeq2 (v1.34.0)<sup>4</sup> was used to perform variance stabilizing transformation for normalization applying the design `~Tissue` with the argument `blind = FALSE` as in [01_varstabDESeq2.R](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/01_varstabDESeq2.R).

---

## Principal Component Analysis 
From the original 297 samples, principal component analysis (PCA) revealed that one sample was a severe outlier as in [02_PCA_UMAP_all.ipynb](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/02_PCA_UMAP_all.ipynb) and [02_PCA_UMAP_peraccession.R](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/02_PCA_UMAP_peraccession.R) for all samples & for each accession separately with FactoMineR (v2.8)<sup>5</sup>. Therefore, this sample (ZDM01467_In2) was removed from further analysis. 

---

## Weighted Gene Correlation Network Analysis (WGCNA)
Co-expression networks were built for all accessions’ single copy orthologs separately using the WGCNA package (v1.69)<sup>6</sup>, as follows in [03_WGCNA_peraccession_1to1orthologs.R](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/03_WGCNA_peraccession_1to1orthologs.R). The soft power threshold used for scale-free topology was determined using the automatic detection function `pickSoftThreshold`. Unsigned networks were calculated with Pearson correlation to obtain the adjacency matrix that accounted for both positive and negative correlations among gene-pairs by applying methods `dynamicTreeCut` and `TOMtype` with a `minimum module size = 30`, `unsigned` network type and otherwise default settings. This resulted in twenty networks, with different module numbers for each accession.

This approach was followed by merging closely clustered modules by a `cut height = 0.1` that resulted in the final set of modules for each accession. Merged modules have been ordered by their size (number of genes) and renamed with accession specific names followed by the order of modules such as “`Akashinriki_1`” or “`B1K_3`”. This resulted in 738 modules across all accessions in total. As WGCNA excluded a few genes with insufficient variation from each of the genotype-specific networks using the `goodSamplesGenes` function, 91 genes have been filtered out from the next steps of analysis, leaving 13,561 single copy genes to analyze further. 

Module eigengenes were calculated using the `moduleEigengenes` function of the WGCNA package, which represents the expression profile of a given module by taking the first principal component of the module. Pearson correlation with Fisher’s-exact significance test was calculated among a binarized tissue table and the eigengenes.

---

## Mercator Functional Enrichment Analysis
Mercator4 (v5.0)<sup>7</sup> functional annotation was performed and enrichment with over-representation analysis of sets of module-genes was done using R package clusterProfiler (v4.6)<sup>8</sup> with Benjamini-Hochberg FDR correction p-value cut-off 0.05 as in [04_Mercator_FuncEnrich.R](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/04_Mercator_FuncEnrich.R).

---

## Comparison of Co-expression Networks with Jaccard-distance and Cosine-similarity 
To compare the shared number of orthologs across accession-specific network modules, Jaccard distance and Cosine-similarity was calculated for each module-pairs, and as Cosine-similarity proved to be more sensitive against module size dependent comparisons, it was used in the next step of the analysis as in [05_Network_comparison_Jaccard_Cosine.py](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/05_Network_comparison_Jaccard_Cosine.py) using Python (v3.10.12)<sup>9</sup>.

---

## Louvain Method for Community Detection 
Across all accessions, a meta-network was established by computing edges using the Cosine similarity values and nodes as modules using networkx (v3.1)<sup>10</sup>, and clustered further into 6 distinct communities via Louvain-community detection<sup>11</sup>. Final network layout was calculated using the Netgraph (v4.12.11)<sup>12</sup> package with edge bundle option as in [06_Networkx_louvain_define_communities.py](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/06_Networkx_louvain_define_communities.py).

---

## OrthoFinder Analysis

Ortholog detection among Oryza sativa (IRGSP-1.0)<sup>13</sup>, Arabidopsis thaliana (TAIR10.57)<sup>14</sup> and selected barley accessions Morex, Akashinriki and B1K representative proteomes was performed using OrthoFinder (2.5.5)<sup>15</sup> with default parameters.  

---

## References

* (1) [R Core Team (2022). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.](https://www.R-project.org/)
* (2) [Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L. D. A., François, R., ... & Yutani, H. 2019 Welcome to the Tidyverse. Journal of open source software, 4(43), 1686.](https://joss.theoj.org/papers/10.21105/joss.01686)
* (3) Kluyver, T., Ragan-Kelley, B., Fernando Perez, Granger, B., Bussonnier, M., Frederic, J., … Willing, C. (2016). Jupyter Notebooks – a publishing format for reproducible computational workflows. In F. Loizides & B. Schmidt (Eds.), Positioning and Power in Academic Publishing: Players, Agents and Agendas (pp. 87–90).
* (4) [Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology, 15(12), 1-21.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8?ref=https://githubhelp.com)
* (5) [Lê, S., Josse, J., & Husson, F. (2008). FactoMineR: an R package for multivariate analysis. Journal of statistical software, 25, 1-18.](https://www.jstatsoft.org/article/view/v025i01)
* (6) [Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. BMC bioinformatics, 9(1), 1-13.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559?ref=https://githubhelp.com)
* (7) [Lohse, M., Nagel, A., Herter, T., May, P., Schroda, M., Zrenner, R., ... & Usadel, B. (2014). M ercator: a fast and simple web server for genome scale functional annotation of plant sequence data (Vol. 37, No. 5, pp. 1250-1258).](https://onlinelibrary.wiley.com/doi/full/10.1111/pce.12231)
* (8) [Wu, T., Hu, E., Xu, S., Chen, M., Guo, P., Dai, Z., ... & Yu, G. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The innovation, 2(3).](https://www.cell.com/the-innovation/pdf/S2666-6758(21)00066-7.pdf)
* (9) [Van Rossum, G., & Drake, F. L. (1995). Python reference manual (Vol. 111, pp. 1-52). Amsterdam: Centrum voor Wiskunde en Informatica.](https://www.cse.unsw.edu.au/~en1811/python-docs/python-3.8.14-docs-pdf/tutorial.pdf)
* (10) [Hagberg, A., Swart, P., & S Chult, D. (2008). Exploring network structure, dynamics, and function using NetworkX (No. LA-UR-08-05495; LA-UR-08-5495). Los Alamos National Lab.(LANL), Los Alamos, NM (United States).](https://www.osti.gov/biblio/960616)
* (11) [De Meo, P., Ferrara, E., Fiumara, G., & Provetti, A. (2011, November). Generalized louvain method for community detection in large networks. In 2011 11th international conference on intelligent systems design and applications (pp. 88-93). IEEE.](https://ieeexplore.ieee.org/abstract/document/6121636?casa_token=y6ZSooOQw78AAAAA:OmLzOZjY9IleAl4cqaKRoLrvHz2VJYbK-JnBrDPHcc5v3DCHemsvqd0-LAk0QU1eEo8_r8mg0Q)
* (12) [Brodersen, P. J. (2023). Netgraph: Publication-quality Network Visualisations in Python. The Journal of Open Source Software, 8(87), 5372.](https://joss.theoj.org/papers/10.21105/joss.05372.pdf)
* (13) [Sakai, H., Lee, S. S., Tanaka, T., Numa, H., Kim, J., Kawahara, Y., ... & Itoh, T. (2013). Rice Annotation Project Database (RAP-DB): an integrative and interactive database for rice genomics. Plant and Cell Physiology, 54(2), e6-e6.](https://academic.oup.com/pcp/article/54/2/e6/1876125)
* (14) [Cheng, C. Y., Krishnakumar, V., Chan, A. P., Thibaud‐Nissen, F., Schobel, S., & Town, C. D. (2017). Araport11: a complete reannotation of the Arabidopsis thaliana reference genome. The Plant Journal, 89(4), 789-804.](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.13415)
* (15) [Emms, D. M., & Kelly, S. (2019). OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome biology, 20, 1-14.](https://link.springer.com/article/10.1186/s13059-019-1832-y)

