Methylation data : Bigwig file to Bedgraph:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rtracklayer")
library(rtracklayer)

bw_data <- import("lungcancer-2.bw", format = "BigWig")
head(bw_data)
df <- as.data.frame(bw_data)
head(df)
export(bw_data, "LungCancer2.bedgraph", format = "bedGraph")
