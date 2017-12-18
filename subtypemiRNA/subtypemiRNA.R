library(edgeR)
library(limma)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(rlang)
library(devtools)
library(digest)
install_github('sinhrks/ggfortify', force =T)
library(ggfortify)
library(factoextra)

###loading data#####
pca1<- prcomp(subtype_TCGA_miRNA[,-c(1:54)])
as.factor(subtype_TCGA_miRNA$GeneExp_Subtype)
autoplot(pca1, data = subtype_TCGA_miRNA, 
         colour = 'GeneExp_Subtype', size = 5) + scale_color_brewer(palette='Set1')+ theme_bw()
fviz_eig(pca1)
fviz_pca_ind(pca1,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T     # Avoid text overlapping
)


# Read the data into R
seqdata <- read.delim("data.txt", stringsAsFactors = FALSE)
# Read the sample information into R
sampleinfo_classical <- read.delim("data/sampleinfo_classical.txt")



countdata <- subtype_TCGA_miRNA[,-c(1:54)]
# Look at the output
head(countdata)
# Store gene id as rownames
rownames(countdata) <- subtype_TCGA_miRNA[,1]#invalid???#
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
         colour = 'subtypeinfo', size = 7) + scale_color_brewer(palette='Set1')+ theme_bw()
fviz_eig(pca2)


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
  





##################################################################################
##each subtype comparision

###loading data#####
pca1<- prcomp(subtype_TCGA_miRNA[,-c(1:54)])
as.factor(subtype_TCGA_miRNA$GeneExp_Subtype)
autoplot(pca1, data = subtype_TCGA_miRNA, 
         colour = 'GeneExp_Subtype', size = 5) + scale_color_brewer(palette='Set1')+ theme_bw()
fviz_eig(pca1)
fviz_pca_ind(pca1,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T     # Avoid text overlapping
)


# Read the sample information into R
sampleinfo <- read.delim("data/sampleinfo.txt")

countdata <- subtype_TCGA_miRNA[,-c(1:54)]
# Look at the output
head(countdata)
# Store gene id as rownames
rownames(countdata) <- subtype_TCGA_miRNA[,1]#invalid???#
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
# Let's choose red for classical, blue for mesenchymal, green for neural and purple for proneural
col.all <- c("red", "blue","green","purple")[sampleinfo$GeneExp_Subtype]
data.frame(sampleinfo$GeneExp_Subtype,col.all)



#################
##limma package##
#################

####C v M####
sampleinfo.CvM <- sampleinfo[which(sampleinfo$GeneExp_Subtype == "Classical" | sampleinfo$GeneExp_Subtype == "Mesenchymal"),]
library(readr)
subtype.CvM <- read_csv("C:/Users/chenkaim/Desktop/TCGA&GEO/statistic/subtype/subtypemiRNA/data/subtype_CvM.csv")
countdata.CvM <- subtype.CvM[,-c(1:54)]
countdata.CvM <- as.matrix(countdata.CvM)
logcounts.CvM <- as.matrix(t(countdata.CvM))
class(logcounts.CvM) <- "numeric"
head(logcounts.CvM)

#Three matrix: expression matrix, group matrix, contrast matrix
#Three step: lmfit, eBayes, topTable
group.CvM <- paste(sampleinfo.CvM$GeneExp_Subtype)
group.CvM <- factor(group.CvM)
##group matrix
design.CvM <- model.matrix(~0+group.CvM)
colnames(design.CvM) <- levels(group.CvM)
rownames(design.CvM) <- colnames(logcounts)
design.CvM

##contrast matrix
contrast.matrix.CvM <- makeContrasts(paste0(unique(group.CvM), collapse = "-"), levels = design.CvM)
contrast.matrix.CvM

##step1
fit.CvM <- lmFit(logcounts.CvM, design.CvM)
##step2
fit2.CvM<- contrasts.fit(fit.CvM, contrast.matrix.CvM)
fit2.CvM<- eBayes(fit2.CvM)
dim(fit2.CvM)
##step3
tempOutput.CvM = topTable(fit2.CvM, coef=1, n=Inf)
DEmiR.CvM = na.omit(tempOutput.CvM) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(DEmiR.CvM)
write.csv(DEmiR.CvM, file = "DEmiR_CvM.csv")
####C v M####


####C v N####
sampleinfo.CvN <- sampleinfo[which(sampleinfo$GeneExp_Subtype == "Classical" | sampleinfo$GeneExp_Subtype == "Neural"),]
library(readr)
subtype.CvN <- read_csv("C:/Users/chenkaim/Desktop/TCGA&GEO/statistic/subtype/subtypemiRNA/data/subtype_CvN.csv")
countdata.CvN <- subtype.CvN[,-c(1:54)]
countdata.CvN <- as.matrix(countdata.CvN)
logcounts.CvN <- as.matrix(t(countdata.CvN))
class(logcounts.CvN) <- "numeric"
head(logcounts.CvN)

#Three matrix: expression matrix, group matrix, contrast matrix
#Three step: lmfit, eBayes, topTable
group.CvN <- paste(sampleinfo.CvN$GeneExp_Subtype)
group.CvN <- factor(group.CvN)
##group matrix
design.CvN <- model.matrix(~0+group.CvN)
colnames(design.CvN) <- levels(group.CvN)
rownames(design.CvN) <- colnames(logcounts.CvN)
design.CvN

##contrast matrix
contrast.matrix.CvN <- makeContrasts(paste0(unique(group.CvN), collapse = "-"), levels = design.CvN)
contrast.matrix.CvN

##step1
fit.CvN <- lmFit(logcounts.CvN, design.CvN)
##step2
fit2.CvN<- contrasts.fit(fit.CvN, contrast.matrix.CvN)
fit2.CvN<- eBayes(fit2.CvN)
dim(fit2.CvN)
##step3
tempOutput.CvN = topTable(fit2.CvN, coef=1, n=Inf)
DEmiR.CvN = na.omit(tempOutput.CvN) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(DEmiR.CvN)
write.csv(DEmiR.CvN, file = "DEmiR_CvN.csv")
####C v N####


####C v P####
sampleinfo.CvP <- sampleinfo[which(sampleinfo$GeneExp_Subtype == "Classical" | sampleinfo$GeneExp_Subtype == "Proneural"),]
library(readr)
subtype.CvP <- read_csv("C:/Users/chenkaim/Desktop/TCGA&GEO/statistic/subtype/subtypemiRNA/data/subtype_CvP.csv")
countdata.CvP <- subtype.CvP[,-c(1:54)]
countdata.CvP <- as.matrix(countdata.CvP)
logcounts.CvP <- as.matrix(t(countdata.CvP))
class(logcounts.CvP) <- "numeric"
head(logcounts.CvP)

#Three matrix: expression matrix, group matrix, contrast matrix
#Three step: lmfit, eBayes, topTable
group.CvP <- paste(sampleinfo.CvP$GeneExp_Subtype)
group.CvP <- factor(group.CvP)
##group matrix
design.CvP <- model.matrix(~0+group.CvP)
colnames(design.CvP) <- levels(group.CvP)
rownames(design.CvP) <- colnames(logcounts.CvP)
design.CvP

##contrast matrix
contrast.matrix.CvP <- makeContrasts(paste0(unique(group.CvP), collapse = "-"), levels = design.CvP)
contrast.matrix.CvP

##step1
fit.CvP <- lmFit(logcounts.CvP, design.CvP)
##step2
fit2.CvP<- contrasts.fit(fit.CvP, contrast.matrix.CvP)
fit2.CvP<- eBayes(fit2.CvP)
dim(fit2.CvP)
##step3
tempOutput.CvP = topTable(fit2.CvP, coef=1, n=Inf)
DEmiR.CvP = na.omit(tempOutput.CvP) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(DEmiR.CvP)
write.csv(DEmiR.CvP, file = "DEmiR_CvP.csv")
####C v P####


####M v N####
sampleinfo.MvN <- sampleinfo[which(sampleinfo$GeneExp_Subtype == "Mesenchymal" | sampleinfo$GeneExp_Subtype == "Neural"),]
library(readr)
subtype.MvN <- read_csv("C:/Users/chenkaim/Desktop/TCGA&GEO/statistic/subtype/subtypemiRNA/data/subtype_MvN.csv")
countdata.MvN <- subtype.MvN[,-c(1:54)]
countdata.MvN <- as.matrix(countdata.MvN)
logcounts.MvN <- as.matrix(t(countdata.MvN))
class(logcounts.MvN) <- "numeric"
head(logcounts.MvN)

#Three matrix: expression matrix, group matrix, contrast matrix
#Three step: lmfit, eBayes, topTable
group.MvN <- paste(sampleinfo.MvN$GeneExp_Subtype)
group.MvN <- factor(group.MvN)
##group matrix
design.MvN <- model.matrix(~0+group.MvN)
colnames(design.MvN) <- levels(group.MvN)
rownames(design.MvN) <- colnames(logcounts.MvN)
design.MvN

##contrast matrix
contrast.matrix.MvN <- makeContrasts(paste0(unique(group.MvN), collapse = "-"), levels = design.MvN)
contrast.matrix.MvN

##step1
fit.MvN <- lmFit(logcounts.MvN, design.MvN)
##step2
fit2.MvN<- contrasts.fit(fit.MvN, contrast.matrix.MvN)
fit2.MvN<- eBayes(fit2.MvN)
dim(fit2.MvN)
##step3
tempOutput.MvN = topTable(fit2.MvN, coef=1, n=Inf)
DEmiR.MvN = na.omit(tempOutput.MvN) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(DEmiR.MvN)
write.csv(DEmiR.MvN, file = "DEmiR_MvN.csv")
####M v N####




####M v P####
sampleinfo.MvP <- sampleinfo[which(sampleinfo$GeneExp_Subtype == "Mesenchymal" | sampleinfo$GeneExp_Subtype == "Proneural"),]
library(readr)
subtype.MvP <- read_csv("C:/Users/chenkaim/Desktop/TCGA&GEO/statistic/subtype/subtypemiRNA/data/subtype_MvP.csv")
countdata.MvP <- subtype.MvP[,-c(1:54)]
countdata.MvP <- as.matrix(countdata.MvP)
logcounts.MvP <- as.matrix(t(countdata.MvP))
class(logcounts.MvP) <- "numeric"
head(logcounts.MvP)

#Three matrix: expression matrix, group matrix, contrast matrix
#Three step: lmfit, eBayes, topTable
group.MvP <- paste(sampleinfo.MvP$GeneExp_Subtype)
group.MvP <- factor(group.MvP)
##group matrix
design.MvP <- model.matrix(~0+group.MvP)
colnames(design.MvP) <- levels(group.MvP)
rownames(design.MvP) <- colnames(logcounts.MvP)
design.MvP

##contrast matrix
contrast.matrix.MvP <- makeContrasts(paste0(unique(group.MvP), collapse = "-"), levels = design.MvP)
contrast.matrix.MvP

##step1
fit.MvP <- lmFit(logcounts.MvP, design.MvP)
##step2
fit2.MvP<- contrasts.fit(fit.MvP, contrast.matrix.MvP)
fit2.MvP<- eBayes(fit2.MvP)
dim(fit2.MvP)
##step3
tempOutput.MvP = topTable(fit2.MvP, coef=1, n=Inf)
DEmiR.MvP = na.omit(tempOutput.MvP) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(DEmiR.MvP)
write.csv(DEmiR.MvP, file = "DEmiR_MvP.csv")
####M v P####




####N v P####
sampleinfo.NvP <- sampleinfo[which(sampleinfo$GeneExp_Subtype == "Neural" | sampleinfo$GeneExp_Subtype == "Proneural"),]
library(readr)
subtype.NvP <- read_csv("C:/Users/chenkaim/Desktop/TCGA&GEO/statistic/subtype/subtypemiRNA/data/subtype_NvP.csv")
countdata.NvP <- subtype.NvP[,-c(1:54)]
countdata.NvP <- as.matrix(countdata.NvP)
logcounts.NvP <- as.matrix(t(countdata.NvP))
class(logcounts.NvP) <- "numeric"
head(logcounts.NvP)

#Three matrix: expression matrix, group matrix, contrast matrix
#Three step: lmfit, eBayes, topTable
group.NvP <- paste(sampleinfo.NvP$GeneExp_Subtype)
group.NvP <- factor(group.NvP)
##group matrix
design.NvP <- model.matrix(~0+group.NvP)
colnames(design.NvP) <- levels(group.NvP)
rownames(design.NvP) <- colnames(logcounts.NvP)
design.NvP

##contrast matrix
contrast.matrix.NvP <- makeContrasts(paste0(unique(group.NvP), collapse = "-"), levels = design.NvP)
contrast.matrix.NvP

##step1
fit.NvP <- lmFit(logcounts.NvP, design.NvP)
##step2
fit2.NvP<- contrasts.fit(fit.NvP, contrast.matrix.NvP)
fit2.NvP<- eBayes(fit2.NvP)
dim(fit2.NvP)
##step3
tempOutput.NvP = topTable(fit2.NvP, coef=1, n=Inf)
DEmiR.NvP = na.omit(tempOutput.NvP) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(DEmiR.NvP)
write.csv(DEmiR.NvP, file = "DEmiR_NvP.csv")
####M v P####


######plot DE miR CvM######
par(mfrow=c(1,1))
with(DEmiR.CvM, plot(logFC, -log(adj.P.Val), pch=19, main="Volcano plot of Classical vs Mesenchymal GBM", xlim=c(-2,2)), cex = 1.5 )
with(subset(DEmiR.CvM, abs(logFC)>1), points(logFC,-log(adj.P.Val), pch=19, col="red"))
library(calibrate)
DEmiR.CvM$name <- rownames(DEmiR.CvM)
with(subset(DEmiR.CvM, abs(logFC)>1),  textxy(logFC, -log(adj.P.Val), labs=name, cex=1.5))


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
######plot DE miR CvM######

######plot DE miR CvN######
par(mfrow=c(1,1))
with(DEmiR.CvN, plot(logFC, -log(adj.P.Val), pch=19, main="Volcano plot of Classical vs Neural GBM", xlim=c(-2,2)), cex = 0.9 )
with(subset(DEmiR.CvN, abs(logFC)>1), points(logFC,-log(adj.P.Val), pch=19, col="red"))
library(calibrate)
DEmiR.CvN$name <- rownames(DEmiR.CvN)
with(subset(DEmiR.CvN, abs(logFC)>1),  textxy(logFC, -log(adj.P.Val), labs=name, cex=0.9))
abline(v=c(-1,1), col="red", lty=2, lwd=2)


DEmiR_CvN_select <- c("hsa-miR-219", "hsa-miR-124a", "hsa-miR-338")
DEmiR_CvN_select
DEmiR_CvN_miR <- logcounts.CvN[DEmiR_CvN_select,]
dim(DEmiR_CvN_miR)
head(DEmiR_CvN_miR)
col.CvN <- c("red", "blue","green", "purple")[sampleinfo.CvN$GeneExp_Subtype]
data.frame(sampleinfo.CvN$GeneExp_Subtype,col.CvN)

heatmap.2(DEmiR_CvN_miR,col=rev(morecols(50)),trace="none", 
          main="Classical vs Neural DE miRs",ColSideColors=col.CvN,scale="row", 
          margins = c(5,20), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo.CvN$GeneExp_Subtype), col = unique(col.CvN), lty = 1, lwd= 5, cex=.7)
######plot DE miR CvN######



######plot DE miR CvP######
par(mfrow=c(1,1))
with(DEmiR.CvP, plot(logFC, -log(adj.P.Val), pch=19, main="Volcano plot of Classical vs Proneural GBM", xlim=c(-2.5,2.5)), cex = 0.9 )
with(subset(DEmiR.CvP, abs(logFC)>1), points(logFC,-log(adj.P.Val), pch=19, col="red"))
library(calibrate)
DEmiR.CvP$name <- rownames(DEmiR.CvP)
with(subset(DEmiR.CvP, abs(logFC)>1),  textxy(logFC, -log(adj.P.Val), labs=name, cex=0.9))
abline(v=c(-1,1), col="red", lty=2, lwd=2)


DEmiR_CvP_select <- c("hsa-miR-219", "hsa-miR-124a", "hsa-miR-338", 
                      "hsa-miR-10b", "hsa-miR-204", "hsa-miR-222", "hsa-miR-34a")
DEmiR_CvP_select
DEmiR_CvP_miR <- logcounts.CvP[DEmiR_CvP_select,]
dim(DEmiR_CvP_miR)
head(DEmiR_CvP_miR)
col.CvP <- c("red", "blue","green", "purple")[sampleinfo.CvP$GeneExp_Subtype]
data.frame(sampleinfo.CvP$GeneExp_Subtype,col.CvP)

heatmap.2(DEmiR_CvP_miR,col=rev(morecols(50)),trace="none", 
          main="Classical vs Proneural DE miRs",ColSideColors=col.CvP,scale="row", 
          margins = c(5,20), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo.CvP$GeneExp_Subtype), col = unique(col.CvP), lty = 1, lwd= 5, cex=.7)
######plot DE miR CvP######




######plot DE miR MvN######
par(mfrow=c(1,1))
with(DEmiR.MvN, plot(logFC, -log(adj.P.Val), pch=19, main="Volcano plot of Mesenchymal vs Neural GBM", xlim=c(-2.5,2.5)), cex = 0.9 )
with(subset(DEmiR.MvN, abs(logFC)>1), points(logFC,-log(adj.P.Val), pch=19, col="red"))
library(calibrate)
DEmiR.MvN$name <- rownames(DEmiR.MvN)
with(subset(DEmiR.MvN, abs(logFC)>1),  textxy(logFC, -log(adj.P.Val), labs=name, cex=0.9))
abline(v=c(-1,1), col="red", lty=2, lwd=2)


DEmiR_MvN_select <- c("hsa-miR-124a", "hsa-miR-21")
DEmiR_MvN_select
DEmiR_MvN_miR <- logcounts.MvN[DEmiR_MvN_select,]
dim(DEmiR_MvN_miR)
head(DEmiR_MvN_miR)
col.MvN <- c("red", "blue","green", "purple")[sampleinfo.MvN$GeneExp_Subtype]
data.frame(sampleinfo.MvN$GeneExp_Subtype,col.MvN)

heatmap.2(DEmiR_MvN_miR,col=rev(morecols(50)),trace="none", 
          main="Mesenchymal vs Neural DE miRs",ColSideColors=col.MvN,scale="row", 
          margins = c(5,20), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo.MvN$GeneExp_Subtype), col = unique(col.MvN), lty = 1, lwd= 5, cex=.7)
######plot DE miR MvN######




######plot DE miR MvP######
par(mfrow=c(1,1))
with(DEmiR.MvP, plot(logFC, -log(adj.P.Val), pch=19, main="Volcano plot of Mesenchymal vs Proneural GBM", xlim=c(-2.5,2.5)), cex = 0.9 )
with(subset(DEmiR.MvP, abs(logFC)>1), points(logFC,-log(adj.P.Val), pch=19, col="red"))
library(calibrate)
DEmiR.MvP$name <- rownames(DEmiR.MvP)
with(subset(DEmiR.MvP, abs(logFC)>1),  textxy(logFC, -log(adj.P.Val), labs=name, cex=0.9))
abline(v=c(-1,1), col="red", lty=2, lwd=2)


DEmiR_MvP_select <- c("hsa-miR-124a", "hsa-miR-21", "hsa-miR-222", "hsa-miR-221", "hsa-miR-34a",
                      "hsa-miR-223", "hsa-miR-204", "hsa-miR-10b", "hsa-miR-338", "hsa-miR-219")
DEmiR_MvP_select
DEmiR_MvP_miR <- logcounts.MvP[DEmiR_MvP_select,]
dim(DEmiR_MvP_miR)
head(DEmiR_MvP_miR)
col.MvP <- c("red", "blue","green", "purple")[sampleinfo.MvP$GeneExp_Subtype]
data.frame(sampleinfo.MvP$GeneExp_Subtype,col.MvP)

heatmap.2(DEmiR_MvP_miR,col=rev(morecols(50)),trace="none", 
          main="Mesenchymal vs Proneural DE miRs",ColSideColors=col.MvP,scale="row", 
          margins = c(5,20), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo.MvP$GeneExp_Subtype), col = unique(col.MvP), lty = 1, lwd= 5, cex=.7)
######plot DE miR MvP######


######plot DE miR NvP######
par(mfrow=c(1,1))
with(DEmiR.NvP, plot(logFC, -log(adj.P.Val), pch=19, main="Volcano plot of Neural vs Proneural GBM", xlim=c(-2.5,2.5)), cex = 0.9 )
with(subset(DEmiR.NvP, abs(logFC)>1), points(logFC,-log(adj.P.Val), pch=19, col="red"))
library(calibrate)
DEmiR.NvP$name <- rownames(DEmiR.NvP)
with(subset(DEmiR.NvP, abs(logFC)>1),  textxy(logFC, -log(adj.P.Val), labs=name, cex=0.9))
abline(v=c(-1,1), col="red", lty=2, lwd=2)


DEmiR_NvP_select <- c( "hsa-miR-204")
DEmiR_NvP_select
DEmiR_NvP_miR <- logcounts.NvP[DEmiR_NvP_select,]
dim(DEmiR_NvP_miR)
head(DEmiR_NvP_miR)
col.NvP <- c("red", "blue","green", "purple")[sampleinfo.NvP$GeneExp_Subtype]
data.frame(sampleinfo.NvP$GeneExp_Subtype,col.NvP)

heatmap.2(DEmiR_NvP_miR,col=rev(morecols(50)),trace="none", 
          main="Proneural vs Neural DE miRs",ColSideColors=col.NvP,scale="row", 
          margins = c(5,20), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo.NvP$GeneExp_Subtype), col = unique(col.NvP), lty = 1, lwd= 5, cex=.7)
######plot DE miR NvP######



##################################################################################