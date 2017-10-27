library(edgeR)
library(limma)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(rlang)
library(devtools)
library(digest)
install_github('sinhrks/ggfortify')
library(ggfortify)

###loading data#####
pca1<- prcomp(subtype_TCGA_miRNA[,-c(1:54)])
as.factor(subtype_TCGA_miRNA$GeneExp_Subtype)
autoplot(pca1, data = subtype_TCGA_miRNA, 
         colour = 'GeneExp_Subtype', size = 5) + scale_color_brewer(palette='Set1')+ theme_bw()




# Read the data into R
seqdata <- read.delim("data.txt", stringsAsFactors = FALSE)
# Read the sample information into R
sampleinfo <- read.delim("data/sampleinfo_classical.txt")



countdata <- subtype_TCGA_miRNA[,-c(1:54)]
# Look at the output
head(countdata)
# Store gene id as rownames
rownames(countdata) <- subtype_TCGA_miRNA[,1]
head(countdata)
colnames(countdata)
countdata <- as.matrix(countdata)


countdata <- t(countdata)
logcounts <- countdata
logcounts <- as.matrix(logcounts)
class(logcounts) <- "numeric"

# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts")
abline(h=median(logcounts),col="red")
title("Boxplots of miRNAs over all samples")






levels(sampleinfo$GeneExp_Subtype)
# Let's choose purple for basal and orange for luminal
col.classical <- c("red", "red4")[sampleinfo$GeneExp_Subtype]
data.frame(sampleinfo$GeneExp_Subtype,col.classical)


##############################
##Hierarchical clustering with heatmaps
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

# Plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", 
          main="Top 50 most variable genes across samples",ColSideColors=col.class,scale="row",
          margins = c(5,9), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo$Grade), col = unique(col.grade), lty = 1, lwd= 5, cex=.7)



# Save the heatmap
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", 
          main="Top 50 most variable genes across samples",ColSideColors=col.idh,scale="row",
          margins = c(5,9), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo$IDHstatus), col = unique(col.idh), lty = 1, lwd= 5, cex=.7)
dev.off()



y <- DGEList(logcounts)
plotMDS(y)
##############################


# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,2))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(sampleinfo$IDH1mutation)
plotMDS(y,col=col.idh)


###############
##limma package##
#################
#Three matrix: expression matrix, group matrix, contrast matrix
#Three step: lmfit, eBayes, topTable
group_classical <- paste(sampleinfo$GeneExp_Subtype)
group_classical <- factor(group_classical)
##group matrix
design_classical <- model.matrix(~0+group_classical)
colnames(design_classical) <- levels(group_classical)
rownames(design_classical) <- colnames(logcounts)
design_classical

##contrast matrix
contrast.matrix_classical <- makeContrasts(paste0(unique(group_classical), collapse = "-"), levels = design_classical)
contrast.matrix_classical

##step1
fit_classical <- lmFit(logcounts, design_classical)
##step2
fit2_classical<- contrasts.fit(fit_classical, contrast.matrix_classical)
fit2_classical<- eBayes(fit2_classical)
dim(fit2_classical)
##step3
tempOutput_classical = topTable(fit2_classical, coef=1, n=Inf)
DEmiR_classical = na.omit(tempOutput_classical) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(DEmiR_classical)
write.csv(DEmiR_classical, file = "DEmiR_classical.csv")

##################


# For the volcano plot we have to specify how many of the top genes to hightlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit2,coef=1,highlight=100)


par(mfrow=c(1,1))
with(DEmiR_classical, plot(logFC, -log(P.Value), pch=20, main="Volcano plot of Classical vs non-Classical GBM", xlim=c(-2.5,2.5)), cex = .9 )
with(subset(DEmiR_classical, abs(logFC)>1), points(logFC,-log(P.Value), pch=20, col="red"))
library(calibrate)
DEmiR_classical$name <- rownames(DEmiR_classical)
with(subset(DEmiR_classical, abs(logFC)>1),  textxy(logFC, -log(P.Value), labs=name, cex=.9))


DEmiR_classical_select <- c("hsa-miR-219", "hsa-miR-204")
DEmiR_classical_select
DEmiR_classical_miR <- logcounts[DEmiR_classical_select,]
dim(DEmiR_classical_miR)
head(DEmiR_classical_miR)

heatmap.2(DEmiR_classical_miR,col=rev(morecols(50)),trace="none", 
          main="Classical vs non-Classical DE miRs",ColSideColors=col.classical,scale="row", 
          margins = c(5,20), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo$GeneExp_Subtype), col = unique(col.classical), lty = 1, lwd= 5, cex=.7)





