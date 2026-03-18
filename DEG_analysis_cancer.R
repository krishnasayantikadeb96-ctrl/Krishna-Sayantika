Codes for DEG:
library(gplots)
library(org.Hs.eg.db)
library(RColorBrewer)
library(NMF)
library(edgeR)

seqdata <- read.csv("cancer_corrected.csv", stringsAsFactors = FALSE)
sampleinfo <- read.delim("sampleinfo.txt", stringsAsFactors = TRUE)
head(seqdata)
dim(seqdata)
sampleinfo
dim(seqdata)
countdata <- seqdata[,-(1:1)]
head(countdata)
rownames(countdata) <- seqdata[,1]
head(countdata)
colnames(countdata)
table(colnames(countdata)==sampleinfo$SampleName)
y <- DGEList(countdata)
y
names(y)
group <- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
group
names(y)
y$samples
group <- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
group
group <- factor(group)
group
y$samples$group <- group
y$samples
myCPM <- cpm(countdata)
head(myCPM)
thresh <- myCPM > 0.5
head(thresh)
table(rowSums(thresh))
keep <- rowSums(thresh) >= 2
summary(keep)
plot(myCPM[,1],countdata[,1])
y <- y[keep, keep.lib.sizes=FALSE]
y$samples$lib.size
logcounts <- cpm(y,log=TRUE)
plotMDS(y,top = 200)
y <- calcNormFactors(y)
y$samples
plotMDS(y)
design <- model.matrix(~ 0 + group)
design
colnames(design) <- levels(group)
design
v <- voom(y,design,plot = TRUE)
fit <- lmFit(v)
names(fit)
cont.matrix <- makeContrasts(cancer.TreatedVsControlCervical=Cervical.Treated-Cervical.Control,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
summa.fit <- decideTests(fit.cont)
summary(summa.fit)
volcanoplot(fit.cont,coef=1,highlight=100,names=fit.cont$genes$SYMBOL, main="cancerStage")
degs<-topTable(fit.cont,coef=1,sort.by="p",number = 'inf')
head(degs)
degs['ENSG00000111801.16',]
write.csv(degs,"C:/Users/Lenovo/OneDrive/Documents/DEGs.csv")
