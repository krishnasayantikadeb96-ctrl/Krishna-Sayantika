install.packages(c("msigdbr", "biomaRt", "dplyr", "readr"))
library(msigdbr)
library(biomaRt)
library(dplyr)
library(readr)


#Codes for Transcription factors
gene_df <- read_csv("my_genes.csv")
gene_list <- unique(gene_df$genes)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
tf_query <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "go_id", "name_1006"),
  filters = "go", 
  values = "GO:0003700", 
  mart = ensembl
)
tf_genes <- unique(tf_query$hgnc_symbol)
head(tf_genes)
gene_df <- read.csv("my_genes.csv", stringsAsFactors = FALSE)
gene_list <- unique(gene_df$genes)
tf_result <- data.frame(
  Gene = gene_list,
  Is_Transcription_Factor = gene_list %in% tf_genes
)
print(tf_result)

#Codes for Signalling kinases
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
signaling_kinase_gos <- c(
  "GO:0004672",  
  "GO:0004709",  
  "GO:0004708",  
  "GO:0004707",  
  "GO:0004713",  
  "GO:0004714",  
  "GO:0000165"   
)
sig_kinase_query <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "go_id", "name_1006"),
  filters = "go",
  values = signaling_kinase_gos,
  mart = ensembl
)
signaling_kinase_genes <- unique(sig_kinase_query$hgnc_symbol)
head(signaling_kinase_genes)
length(signaling_kinase_genes)
gene_df <- read.csv("my_genes.csv", stringsAsFactors = FALSE)
gene_list <- unique(gene_df$genes)

signaling_kinase_result <- data.frame(
  Gene = gene_list,
  Is_Signaling_Kinase = gene_list %in% signaling_kinase_genes
)
print(signaling_kinase_result)

# Codes for GPCRs
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gpcr_go_terms <- c("GO:0007186", "GO:0004930", "GO:0038032", "GO:0007218")
gpcr_query <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "go_id", "name_1006"),
  filters = "go",
  values = gpcr_go_terms,
  mart = ensembl
)
gpcr_genes <- unique(gpcr_query$hgnc_symbol)
gene_df <- read.csv("my_genes.csv", stringsAsFactors = FALSE)
gene_list <- unique(gene_df$genes)
gpcr_result <- data.frame(
  Gene = gene_list,
  Is_GPCR_Signaling = gene_list %in% gpcr_genes
)
print(gpcr_result)

#Codes for Growth factors
growth_factor_go <- c("GO:0008083")
growth_query <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "go_id", "name_1006"),
  filters = "go",
  values = growth_factor_go,
  mart = ensembl
)
growth_factor_genes <- unique(growth_query$hgnc_symbol)
gene_df <- read.csv("my_genes.csv", stringsAsFactors = FALSE)
gene_list <- unique(gene_df$genes)
growth_factor_result <- data.frame(
  Gene = gene_list,
  Is_Growth_Factor = gene_list %in% growth_factor_genes
)
print(growth_factor_result)

#Codes for Immune regulatory genes
immune_go_terms <- c("GO:0002376", "GO:0006955", "GO:0002250", "GO:0045087", "GO:0006952")
immune_query <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "go_id", "name_1006"),
  filters = "go",
  values = immune_go_terms,
  mart = ensembl
)
immune_genes <- unique(immune_query$hgnc_symbol)
gene_df <- read.csv("my_genes.csv", stringsAsFactors = FALSE)
gene_list <- unique(gene_df$genes)
immune_result <- data.frame(
  Gene = gene_list,
  Is_Immune_Gene = gene_list %in% immune_genes
)
print(immune_result)

#Codes for Oncogenes
gene_df <- read.csv("my_genes.csv", stringsAsFactors = FALSE)
gene_list <- unique(gene_df$genes)
cgc <- read.csv("cancer_gene_census.csv", stringsAsFactors = FALSE)
oncogene_list <- cgc[cgc$Role.in.Cancer %in% c("oncogene", "oncogene;TSG", "TSG;oncogene"), "Gene.Symbol"]
oncogene_list <- unique(oncogene_list)
oncogene_result <- data.frame(
  Gene = gene_list,
  Is_Oncogene = gene_list %in% oncogene_list
)
print(oncogene_result)

#Codes for Tumor supressor genes
tsg_list <- cgc[cgc$Role.in.Cancer %in% c("TSG", "TSG;oncogene", "oncogene;TSG"), "Gene.Symbol"]
tsg_list <- unique(tsg_list)
tsg_result <- data.frame(
  Gene = gene_list,
  Is_Tumor_Suppressor = gene_list %in% tsg_list
)
print(tsg_result)
