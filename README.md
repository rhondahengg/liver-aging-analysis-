# liver-aging-analysis
This project analyses single-cell RNA sequencing (scRNA-seq) data from human liver samples to investigate age-associated cellular and molecular differences between Pediatrics and Adult donors.

# dataset
Source: Edgar et al. (2025) —
“Single-Cell Atlas Of Human Pediatric Liver Reveals Age-Related Hepatic Gene Signatures”
(CZ CELLxGENE)

The dataset was subsampled to approximately 300 cells per donor to support efficient exploratory analysis and balanced age-group comparisons.

# objectives
This project aims to answer the following biological questions:
1. Do liver cell-type proportions differ between Pediatric and Adult donors?
2. Which genes are differentially expressed in Adult vs Pediatric hepatocytes?
3. Which biological pathways are enriched in age-associated liver gene expression programs?

# overall workflow
1. Preprocessing and Quality Control\
   did the standard preprocessing scRNA-seq using Scanpy
2. Cell-Type Composition Analysis\
   key liver cell populations analyzed: hepatocytes, cholangiocytes, Kupffer cells, hepatic sinusoidal endothelial cells, stellate cells
3. Hepatocyte Differential Expression Analysis\
   cells were subsetted from the full dataset and aggregated into donor-level pseudo-bulk samples
   differential expression analysis was performed between: peds vs adults
4. Functional and Pathway Enrichment Analysis\
   Gene set enrichment analysis (GSEA) and enrichment testing were performed using GSEApy with: Reactome pathways, KEGG pathways, MSigDB      Hallmark pathways, Gene Ontology Biological Processes

### note
liver.py is supposed to be in a jupyter notebook format to be able to load the results and graphs/tables but github says there is an error in .ipynb format 
