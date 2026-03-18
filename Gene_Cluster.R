Codes for Gene Cluster:
library(Biobase)
library(Mfuzz)
library(dplyr)

tpm_data <- read.csv("your_data.csv", row.names = 1)
head(tpm_data)
data_eset <- new("ExpressionSet", exprs = as.matrix(tpm_data))
data_eset <- standardise(data_eset)
m <- mestimate(data_eset)
print(m)

anyNA(exprs(data_eset))
any(is.nan(exprs(data_eset)))
any(is.infinite(exprs(data_eset)))
exprs(data_eset) <- apply(exprs(data_eset), 2, function(x) {
  x[is.na(x) | is.nan(x) | is.infinite(x)] <- mean(x, na.rm = TRUE)
  return(x)
})
data_eset <- standardise(data_eset)
Dmin_values <- Dmin(data_eset, m = m, crange = seq(2, 20, 1))
plot(seq(2, 20, 1), Dmin_values, type = "b", pch = 19, col = "blue",
     xlab = "Number of clusters", ylab = "Dmin", main = "Optimal Number of Clusters_Breast_Treated")
c <- 20
cl <- mfuzz(data_eset, c = c, m = m)
mfuzz.plot(data_eset, cl = cl, mfrow = c(2, 2), time.labels = colnames(tpm_data))
cluster_members <- cl$cluster
write.csv(cluster_members, "cluster_members.csv")
