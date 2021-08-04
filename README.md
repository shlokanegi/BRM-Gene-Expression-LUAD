# BRM-Gene-Expression-LUAD
# Project Overview
## Gene Expression Profiles, Biological Pathways, and Disease Characteristics associated with BRM Expression Levels in Lung Adenocarcinoma Tumors

### Overview
BRM (encoded by _SMARCA2_) is one of the catalytic units of the SWI/SNF chromatin remodeler. This protein and its gene variations have been found to be associated with the risk of susceptibility and prognosis in multiple human cancers. We aimed to investigate whether gene expression levels of BRM in lung adenocarcinoma (LUAD) tumor tissues were associated with expression profiles of other genes, which pathways these genes clustered in, and whether they were associated with the survival times of LUAD patients in the TCGA dataset. Our findings revealed genes whose expression levels correlated with that of BRM. GO and pathway enrichment analysis under medium and high classification stringency respectively, unanimously revealed mitochondrial translation as most enriched pathway. Finally, survival analysis showed that _NFIX, TUBA1C, RASSF5, DTYMK_ and _CKS1B_ were significantly associated with survival outcomes in LUAD patients.

<img width="1433" alt="MITACS-poster-final" src="https://user-images.githubusercontent.com/66521525/128172502-d0f4c14c-12c0-4134-9a4e-fb0d9f9f144e.png">


### IMPLEMENTATION is divided into the following parts :-
#### 1) Data Collection and Pre-processing
RSEM Normalized RNA-Seq data from The Cancer Genome Atlas (TCGA) Lung Adenocarcinoma (LUAD) cohort, was downloaded from the GDAC Firehose Website (Broad GDAC Firehose). 

#### 2) Data Categorization
LUAD tumor samples were categorized into three groups, based on their _SMARCA2_ gene expression: the lowest 25 percentile samples having the least _SMARCA2_ gene expression values were classified as BRM-Low tumor samples, the middle 26 to 75 percentile samples with medium _SMARCA2_ gene expression values were classified as BRM-Intermediate tumor samples and top 75 to 100 percentile samples having the highest _SMARCA2_ gene expression values were classified as BRM-High tumor samples. Intermediate Samples were removed from further analysis 

#### 3) Gene Expression Analysis
The _SMARCA2_ gene expression was analysed in LUAD tumor samples. Histograms were created using the hist() function from the base package in R. We compared the expression of all genes between BRM-High and BRM-Low groups to identify genes associated with BRM gene expression. For this, we performed the Mann-Whitney U Test, using the wilcox.test() function from the stats package in R, and adjusted the p-values using Bonferroni Correction. Threshold for significance was set to 2.82E-06 after correction. Most significant genes were studied further for their role in cancer. Boxplots for top 15 significant genes were created using the ggboxplot() function requiring the ggplot and ggpubr package in R.

#### 4) Gene Ontology and Reactome Pathway Enrichment Analysis
We analysed the top 3000 genes using the Functional Annotation Clustering tool of Database for Annotation, Visualization and Integrated Discovery (DAVID) online tool (DAVID Functional Annotation Bioinformatics Microarray Analysis). GO clusters were obtained under medium classification stringency, whereas the Reactome Pathway clusters were obtained under high classification stringency. Further, we performed Overrepresentation Analysis using the Reactome online tool (Reactome Pathway Database) to visualize the coverage of our enriched pathways by the clustered genes.

#### 4) Survival Analysis
The OncoLnc online tool (OncoLnc) was utilized to study the overall survival (OS) of the top 15 significant genes. Genes with FDR corrected Cox regression p-value < 0.05 were considered to be significantly related to survival outcomes.  


### RESULTS and CONCLUSION
* BRM expression levels are variable across LUAD tumors
* 6034 genes were associated with BRM expression levels
* Mitochondrial Translation was found to be the most enriched pathway
* _NFIX_, _TUBA1C_, _RASSF5_, _CKS1B_ and _DTYMK_ were associated with survival outcomes.








