Methylation data : IDAT file to Bedgraph:
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

idat_dir <- "file path/idat_folder"
rgSet <- read.metharray.exp(base = idat_dir)
mSet <- preprocessNoob(rgSet)

qc <- getQC(mSet)
plotQC(qc)

beta <- getBeta(mSet)
annots <- getAnnotation(mSet)

annots_matched <- annots[match(rownames(beta), rownames(annots)), ]

bedgraph_data <- data.frame(
  chrom = annots_matched$chr,
  start = annots_matched$pos - 1,
  end = annots_matched$pos,
  value = beta[, 1] 
)

bedgraph_data <- bedgraph_data[!is.na(bedgraph_data$chrom) & !is.na(bedgraph_data$value), ]
bedgraph_data <- bedgraph_data[grepl("^chr", bedgraph_data$chrom), ]
bedgraph_data <- bedgraph_data[order(bedgraph_data$chrom, bedgraph_data$start), ]

write.table(
  bedgraph_data,
  file = "C:/Users/Lenovo/OneDrive/Documents/GSM7658416.bedgraph",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
