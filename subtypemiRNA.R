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
sampleinfo_classical <- read.delim("data/sampleinfo_classical.txt")



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





#######################################
##mesenchymal
sampleinfo_mesenchymal <- read.delim("data/sampleinfo_mesenchymal.txt")

levels(sampleinfo_mesenchymal$GeneExp_Subtype)
# Let's choose purple for basal and orange for luminal
col.mesenchymal <- c("blue4", "blue")[sampleinfo_mesenchymal$GeneExp_Subtype]
data.frame(sampleinfo_mesenchymal$GeneExp_Subtype,col.mesenchymal)


###############
##limma package##
#################
#Three matrix: expression matrix, group matrix, contrast matrix
#Three step: lmfit, eBayes, topTable
group_mesenchymal <- paste(sampleinfo_mesenchymal$GeneExp_Subtype)
group_mesenchymal <- factor(group_mesenchymal)
##group matrix
design_mesenchymal <- model.matrix(~0+group_mesenchymal)
colnames(design_mesenchymal) <- levels(group_mesenchymal)
rownames(design_mesenchymal) <- colnames(logcounts)
design_mesenchymal

##contrast matrix
contrast.matrix_mesenchymal <- makeContrasts(paste0(unique(group_mesenchymal), collapse = "-"), levels = design_mesenchymal)
contrast.matrix_mesenchymal

##step1
fit_mesenchymal <- lmFit(logcounts, design_mesenchymal)
##step2
fit2_mesenchymal<- contrasts.fit(fit_mesenchymal, contrast.matrix_mesenchymal)
fit2_mesenchymal<- eBayes(fit2_mesenchymal)
dim(fit2_mesenchymal)

tempOutput_mesenchymal = topTable(fit2_mesenchymal, coef=1, n=Inf)
DEmiR_mesenchymal = na.omit(tempOutput_mesenchymal) 

head(DEmiR_mesenchymal)
write.csv(DEmiR_mesenchymal, file = "DEmiR_mesenchymal.csv")


volcanoplot(fit2,coef=1,highlight=100)


with(DEmiR_mesenchymal, plot(logFC, -log(P.Value), pch=20, main="Volcano plot of Mesenchymal vs non-Mesenchymal GBM", xlim=c(-2.5,2.5)), cex = 1 )
with(subset(DEmiR_mesenchymal, abs(logFC)>1), points(logFC,-log(P.Value), pch=20, col="red"))
library(calibrate)
rownames(DEmiR_mesenchymal) <- DEmiR_mesenchymal$X1
DEmiR_mesenchymal$name <- rownames(DEmiR_mesenchymal)
with(subset(DEmiR_mesenchymal, abs(logFC)>1),  textxy(logFC, -log(P.Value), labs=name, cex=1))


DEmiR_mesenchymal_select <- c("hsa-miR-21", "hsa-miR-222")
DEmiR_mesenchymal_select
DEmiR_mesenchymal_miR <- logcounts[DEmiR_mesenchymal_select,]
dim(DEmiR_mesenchymal_miR)
head(DEmiR_mesenchymal_miR)

heatmap.2(DEmiR_mesenchymal_miR,col=rev(morecols(50)),trace="none", 
          main="mesenchymal vs non-mesenchymal DE miRs",ColSideColors=col.mesenchymal,scale="row", 
          margins = c(5,20), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo_mesenchymal$GeneExp_Subtype), col = unique(col.mesenchymal), lty = 1, lwd= 5, cex=.7)
#######################################


#######################################
##Neural
sampleinfo_neural <- read.delim("data/sampleinfo_neural.txt")

levels(sampleinfo_neural$GeneExp_Subtype)
# Let's choose purple for basal and orange for luminal
col.neural <- c("green", "green4")[sampleinfo_neural$GeneExp_Subtype]
data.frame(sampleinfo_neural$GeneExp_Subtype,col.neural)


###############
##limma package##
#################
#Three matrix: expression matrix, group matrix, contrast matrix
#Three step: lmfit, eBayes, topTable
group_neural <- paste(sampleinfo_neural$GeneExp_Subtype)
group_neural <- factor(group_neural)
##group matrix
design_neural <- model.matrix(~0+group_neural)
colnames(design_neural) <- levels(group_neural)
rownames(design_neural) <- colnames(logcounts)
design_neural

##contrast matrix
contrast.matrix_neural <- makeContrasts(paste0(unique(group_neural), collapse = "-"), levels = design_neural)
contrast.matrix_neural

##step1
fit_neural <- lmFit(logcounts, design_neural)
##step2
fit2_neural<- contrasts.fit(fit_neural, contrast.matrix_neural)
fit2_neural<- eBayes(fit2_neural)
dim(fit2_neural)

tempOutput_neural = topTable(fit2_neural, coef=1, n=Inf)
DEmiR_neural = na.omit(tempOutput_neural) 

head(DEmiR_neural)
write.csv(DEmiR_neural, file = "DEmiR_neural.csv")


volcanoplot(fit2,coef=1,highlight=100)


with(DEmiR_neural, plot(logFC, -log(P.Value), pch=20, main="Volcano plot of neural vs non-neural GBM", xlim=c(-2.5,2.5)), cex = 1 )
with(subset(DEmiR_neural, abs(logFC)>1), points(logFC,-log(P.Value), pch=20, col="red"))
library(calibrate)
rownames(DEmiR_neural) <- DEmiR_neural$X1
DEmiR_neural$name <- rownames(DEmiR_neural)
with(subset(DEmiR_neural, abs(logFC)>1),  textxy(logFC, -log(P.Value), labs=name, cex=1))


DEmiR_neural_select <- c("hsa-miR-219", "hsa-miR-124a")
DEmiR_neural_select
DEmiR_neural_miR <- logcounts[DEmiR_neural_select,]
dim(DEmiR_neural_miR)
head(DEmiR_neural_miR)

heatmap.2(DEmiR_neural_miR,col=rev(morecols(50)),trace="none", 
          main="Neural vs non-Neural DE miRs",ColSideColors=col.neural,scale="row", 
          margins = c(5,20), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo_neural$GeneExp_Subtype), col = unique(col.neural), lty = 1, lwd= 5, cex=.7)
#######################################


#######################################
##Proneural
sampleinfo_proneural <- read.delim("data/sampleinfo_proneural.txt")

levels(sampleinfo_proneural$GeneExp_Subtype)
# Let's choose purple for basal and orange for luminal
col.proneural <- c("purple4", "purple")[sampleinfo_proneural$GeneExp_Subtype]
data.frame(sampleinfo_proneural$GeneExp_Subtype,col.proneural)


###############
##limma package##
#################
#Three matrix: expression matrix, group matrix, contrast matrix
#Three step: lmfit, eBayes, topTable
group_proneural <- paste(sampleinfo_proneural$GeneExp_Subtype)
group_proneural <- factor(group_proneural)
##group matrix
design_proneural <- model.matrix(~0+group_proneural)
colnames(design_proneural) <- levels(group_proneural)
rownames(design_proneural) <- colnames(logcounts)
design_proneural

##contrast matrix
contrast.matrix_proneural <- makeContrasts(paste0(unique(group_proneural), collapse = "-"), levels = design_proneural)
contrast.matrix_proneural

##step1
fit_proneural <- lmFit(logcounts, design_proneural)
##step2
fit2_proneural<- contrasts.fit(fit_proneural, contrast.matrix_proneural)
fit2_proneural<- eBayes(fit2_proneural)
dim(fit2_proneural)

tempOutput_proneural = topTable(fit2_proneural, coef=1, n=Inf)
DEmiR_proneural = na.omit(tempOutput_proneural) 

head(DEmiR_proneural)
write.csv(DEmiR_proneural, file = "DEmiR_proneural.csv")


volcanoplot(fit2,coef=1,highlight=100)


with(DEmiR_proneural, plot(logFC, -log(P.Value), pch=20, main="Volcano plot of proneural vs non-proneural GBM", xlim=c(-2.5,2.5)), cex = 1 )
with(subset(DEmiR_proneural, abs(logFC)>1), points(logFC,-log(P.Value), pch=20, col="red"))
library(calibrate)
rownames(DEmiR_proneural) <- DEmiR_proneural$X1
DEmiR_proneural$name <- rownames(DEmiR_proneural)
with(subset(DEmiR_proneural, abs(logFC)>1),  textxy(logFC, -log(P.Value), labs=name, cex=0.8))


DEmiR_proneural_select <- c("hsa-miR-124a", "hsa-miR-219", "hsa-miR-10b","hsa-miR-338",
                         "hsa-miR-34a","hsa-miR-204","hsa-miR-222")
DEmiR_proneural_select
DEmiR_proneural_miR <- logcounts[DEmiR_proneural_select,]
dim(DEmiR_proneural_miR)
head(DEmiR_proneural_miR)

heatmap.2(DEmiR_proneural_miR,col=rev(morecols(50)),trace="none", 
          main="Proneural vs non-Proneural DE miRs",ColSideColors=col.proneural,scale="row", 
          margins = c(5,20), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo_proneural$GeneExp_Subtype), col = unique(col.proneural), lty = 1, lwd= 5, cex=.7)
#######################################



DEmiR_all_select <- c("hsa-miR-124a", "hsa-miR-219", "hsa-miR-10b","hsa-miR-338",
                      "hsa-miR-34a","hsa-miR-204","hsa-miR-222","hsa-miR-21")
DEmiR_all_select
DEmiR_all_miR <- logcounts[DEmiR_all_select,]
dim(DEmiR_all_miR)
head(DEmiR_all_miR)

sampleinfo_all <- read.delim("data/sampleinfo.txt")
col.all <- c("red", "blue", "green", "purple")[sampleinfo_all$GeneExp_Subtype]
data.frame(sampleinfo_all$GeneExp_Subtype,col.all)

heatmap.2(DEmiR_all_miR,col=rev(morecols(50)),trace="none", 
          main="Subtype specific miRs",ColSideColors=col.all,scale="row", 
          margins = c(5,20), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo_all$GeneExp_Subtype), col = unique(col.all), lty = 1, lwd= 5, cex=.7)



####PCA
subtype_TCGA_miRNA_all_select <- subtype_TCGA_miRNA[,DEmiR_all_select]
pca2<- prcomp(subtype_TCGA_miRNA_all_select)
subtypeinfo <- subtype_TCGA_miRNA$GeneExp_Subtype
subtype_TCGA_miRNA_all_select <- cbind(subtype_TCGA_miRNA_all_select, subtypeinfo)
as.factor(subtype_TCGA_miRNA_all_select$subtypeinfo)
autoplot(pca2, data = subtype_TCGA_miRNA_all_select, 
         colour = 'subtypeinfo', size = 5) + scale_color_brewer(palette='Set1')+ theme_bw()


####tSNE
library(caret)  
library(Rtsne)

##loading data
data_tsne <- subtype_TCGA_miRNA_all_select
class(data_tsne) <- "numeric"
##Rtsne function
set.seed(9)
tsne_model1 <- Rtsne(data_tsne, check_duplicates = F, pca = T, perplexity = 50, theta = 0.5, dims = 2)

## getting the two dimension matrix
d_tsne_1 <- as.data.frame(tsne_model1$Y)

## plotting the results without clustering
ggplot(d_tsne_1, aes(x = V1, y = V2))+
  geom_point(size = 2)+
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_brewer(palette = "Set2")

## binding clinical info from original data frame
sampleinfo_all 
d_tsne_2 <- cbind(d_tsne_1, sampleinfo_all)

## Plotting
ggplot(d_tsne_2, aes(x = V1, y = V2, colour = GeneExp_Subtype))+
  geom_point(size = 7)+
  theme_light(base_size=20)  +
  ggtitle("t-SNE")+
  guides(colour=guide_legend(override.aes=list(size=5)))+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) 
  
