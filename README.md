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

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>
[GitHub Pages](https://pages.github.com/)
---

## Preparation of count & metadata tables
The R script [00_DataWrangling.ipynb](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/00_DataWrangling.ipynb) demonstrates how to divide the overall count table (`PanBaRT20_tpm_genes.csv`) from all 297 samples into accession specific sample tables and how to subset the single-copy orthologs from all genes, using the "core;single" keywords on the `GsRTD_and_PanBaRT20_match.tsv` table. The generated `"PanBaRT20_geneTPM_ort_filt_*.csv"` and `"*_meta.csv"` tables serve as input tables for the next step of the analysis.

Out of the total number of 13,680 single copy genes, 13,652 genes have passed the filter of lowly expressing genes (TPM > 0.5 in minimum 2 biological replicates) in all accessions and were taking part of the following analysis. 

---

## Variance Stabilizing Transformation
R package DESeq2 (v1.34.0) is used to perform variance stabilizing transformation for normalization applying the design “~Tissue” with the argument “blind = FALSE” in the [01_varstabDESeq2.R](https://github.com/vanda-marosi/PanBarleyNetworks/blob/main/scripts/01_varstabDESeq2.R) script <sup>1</sup>.

---

## Principal Component Analysis 
From the original 297 samples, principal component analysis (PCA) revealed that one sample was a severe outlier with R scripts 02_PCA_UMAP_all.ipynb and 02_PCA_UMAP_peraccession.R for all samples & for each accession separately. Therefore, this sample (ZDM01467_In2) was removed from further analysis. 

---

## Weighted Gene Correlation Network Analysis
Co-expression networks were built for all accessions’ single copy orthologs separately using the WGCNA package (v1.69) (https://doi.org/10.1186/1471-2105-9-559), as follows. The soft power threshold used for scale-free topology was determined using the automatic detection function “pickSoftThreshold”. Unsigned networks were calculated with Pearson correlation to obtain the adjacency matrix that accounted for both positive and negative correlations among gene-pairs by applying methods dynamicTreeCut and TOMtype with a minimum module size = 30, “unsigned” network type and otherwise default settings. This resulted in twenty networks, with different module numbers for each accession.

This approach was followed by merging closely clustered modules by a cut height 0.1 that resulted in the final set of modules for each accession. Merged modules have been ordered by their size (number of genes) and renamed with accession specific names followed by the order of modules such as “Akashinriki_1” or “B1K_3”. This resulted in 738 modules across all accessions in total. As WGCNA excluded a few genes with insufficient variation from each of the genotype-specific networks using the “goodSamplesGenes” function, 91 genes have been filtered out from the next steps of analysis, leaving 13,561 single copy genes to analyze further. 

Module eigengenes were calculated using the “moduleEigengenes” function of the WGCNA package, which represents the expression profile of a given module by taking the first principal component of the module. Pearson correlation with Fisher’s-exact significance test was calculated among a binarized tissue table and the eigengenes.

---

## Mercator Functional Enrichment Analysis
Mercator4 (v5.0) (PMID: 34448161) functional annotation was performed and enrichment with over-representation analysis (ORA) of sets of module-genes was done using R package clusterProfiler (v4.6) (doi:10.1016/j.xinn.2021.100141) with Benjamini-Hochberg FDR correction p-value cut-off 0.05

---

## Comparison of Co-expression Networks with Jaccard-distance and Cosine-similarity 
To compare the shared number of orthologs across accession-specific network modules, Jaccard distance and Cosine-similarity was calculated for each module-pairs, and as Cosine-similarity proved to be more sensitive against small-to-big module comparisons, it was used in the next step of the analysis.

---

## Louvain Method for Community Detection 
Across all accessions, a meta-network was established by computing edges using the Cosine values and nodes as modules using networkx (v) (Hagberg, A., Swart, P., & S Chult, D. (2008). Exploring network structure, dynamics, and function using NetworkX (No. LA-UR-08-05495; LA-UR-08-5495). Los Alamos National Lab. (LANL), Los Alamos, NM (United States).), and clustered further into 6 distinct communities via the Louvain method for community detection (https://iopscience.iop.org/article/10.1088/1742-5468/2008/10/P10008). Final network layout was calculated using the Netgraph (v) (https://doi.org/10.21105/joss.05372) package with edge bundle option.

---

## OrthoFinder Analysis
In order to identify well-known rice orthologs

---

## References

* 1 Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8.
*
*
*
*

