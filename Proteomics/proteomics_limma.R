# running limma on two batches of AREDS proteomics data; tests with late AMD and any AMD per person. PCAs.

library(impute);library(data.table);library(limma);library(Hmisc);library(RUVSeq);library(Biobase)

prot1 <- fread('~/Proteomics/report.pg_matrix_batch2only.txt') # 1270 proteins.
prot2 <- fread('~/Proteomics/report.pg_matrix_original.csv') # 1497 proteins.

# find intersection of genes.----

length(intersect(prot1$Protein.Names, prot2$Protein.Names)) #  1049 genes are shared, in terms of specific gene symbol assigned.

length(intersect(prot1$Protein.Group, prot2$Protein.Group)) #  1049

# more in depth search for matching protein names----

library(dplyr)
library(tidyverse)
#View(prot1 %>% 
#       separate_rows(Protein.Names,  sep = ";"))

prot1 <- prot1 %>% separate_rows(Genes,  sep = ";") # 1336 genes
prot2 <- prot2 %>% separate_rows(Genes,  sep = ";") # 1567 genes


prot1 <- subset(prot1, Genes %in% intersect(prot1$Genes, prot2$Genes))
prot2 <- subset(prot2, Genes %in% intersect(prot1$Genes, prot2$Genes)) # 1110

# concatenate the batches.----

prot <- merge(prot1[,c(3,4,6:65)], prot2[,c(3,4,6:45)], by=c('Protein.Names','Genes') ) # 1074 genes, 100 samples
names(prot)

# clean up protein names. there's 5 IG* family genes with missing gene symbols here.----
prot <- as.data.table(prot)
prot[Protein.Names=='IGA2_HUMAN',Genes:='IGHA2']
prot[Protein.Names=='IGD_HUMAN',Genes:='IGHD']
prot[Protein.Names=='IGE_HUMAN',Genes:='IGHE']
prot[Protein.Names=='IGG1_HUMAN',Genes:='IGHG1']
prot[Protein.Names=='IGK_HUMAN',Genes:='IGK']

prot[,Protein.Names:=NULL]

# clean up sample names.----

names(prot) <- gsub('_depleteDIA_3uL.raw','',names(prot))
names(prot) <- gsub('NEYI00','NEYI-00', names(prot))


# get metadata.----
metabo.meta <- fread('~/Metabolomics/AREDS/AREDS_allbatches_metadata.csv')
metabo.meta[,c('batch','BOX_NUMBER'):=NULL]

samples <- data.table(PARENT_SAMPLE_NAME = names(prot)[2:101])

metadata <- merge(samples, metabo.meta, by='PARENT_SAMPLE_NAME')

orig <- gsub('_depleteDIA_3uL.raw','',names(prot1)[6:65])
orig <- gsub('NEYI00','NEYI-00',orig)

new <- unique(samples$PARENT_SAMPLE_NAME)[unique(samples$PARENT_SAMPLE_NAME) %nin% orig]

samples[,batch:=ifelse(PARENT_SAMPLE_NAME %in% new, '2', '1')]

## get age and sex----
meta <- fread('~/Proteomics/Metadata.csv')
meta <- unique(meta[,c('PARENT_SAMPLE_NAME','CLIENT_SAMPLE_ID','age','male','smoked','edu')])

metadata <- merge(metadata, meta, by=c('PARENT_SAMPLE_NAME','CLIENT_SAMPLE_ID'))
#fwrite(metadata,'~/Proteomics/Proteomics_metadata.csv')

# pca of samples to verify that there is one very far outlier sample----
pcadata <- as.data.frame(unique(t(prot[,-c(1)])))
names(pcadata) <- prot$Genes

library(PCAtools)
library(factoextra)
library(FactoMineR)

res.pca <- PCA(pcadata, graph = FALSE)
eigenvalues <- res.pca$eig
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")

lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], 
      type="b", pch=19, col = "red")


fviz_screeplot(res.pca, ncp=10)


plot.PCA(res.pca, axes = c(2,3), choix=c("ind"))
plot.PCA(res.pca, axes = c(2,3), choix=c("var"))


# imputation-----
imput <- impute.knn(as.matrix(prot[,-c(1)])) #99 genes with more than 50 % entries missing;mean imputation used for these rows
imput <- imput$data # large numbers like 1*10^7

# run RUVseq to get RUV variable-----
imput.t <- as.data.table(t(imput))

genes <- prot$Genes
names(imput.t) <- genes

assayData <-imput.t
assayData <- as.data.table(t(assayData))
assayData <- assayData %>% mutate_if(is.numeric, round)
names(assayData) <- names(prot)[-1]
assayData <- as.matrix(assayData)
row.names(assayData) <- genes

spikes <- fread('~/Proteomics/Spikes_list.csv') # 418 determined from a previous GLM to find genes where batch is significantly different but not late AMD.
spikes <- unlist(spikes$GeneID, use.names = F)
spikes <- spikes[spikes %in% genes] # 188

corrected <- RUVg(as.matrix(assayData), spikes, k=1)
#save(corrected, file='~/Proteomics/Limma/Limma_DIA100/RUV_corrections.RData')

# before correction:
plotRLE(assayData, outline=FALSE, ylim=c(-1, 1))

# after correction:
plotRLE(corrected$normalizedCounts, outline=FALSE, ylim=c(-1, 1))

plotPCA(corrected$normalizedCounts, cex=1.2, )


# scaling and centering ---------------
imput <- t(scale(t(as.matrix(imput))))

imput <- as.data.frame(imput)
imput$Genes <- genes
imput <- subset(imput, !duplicated(Genes)) # 1070
genes <- imput$Genes
imput <- imput[,1:100]
row.names(imput) <- genes

# limma setup----
all(samples$PARENT_SAMPLE_NAME == names(imput))
samples$order <- 1:nrow(samples)
samples <- merge(metadata,samples, by='PARENT_SAMPLE_NAME')
samples <- samples[order(order)]

metadata <- samples
metadata[,anyAMD_person:=ifelse(SevScaleAvg > 1, '1','0')]
metadata[,lateAMD_person:=ifelse(SevScaleAvg > 9, '1','0')]

pheno <- metadata[,-c('order')]
limmadat <- t(imput)

pheno <- as.data.frame(pheno)
row.names(pheno) <- pheno$PARENT_SAMPLE_NAME
pheno <- subset(pheno, select=c("PARENT_SAMPLE_NAME", 'CLIENT_SAMPLE_ID',"age", "male", 'batch', "lateAMD_person", "anyAMD_person", 'TIME_POINT', 'SevScaleAvg','smoked','edu'))

pheno$RUV <- corrected$W

limmadat <- as.data.frame(limmadat)
row.names(limmadat) <- pheno$PARENT_SAMPLE_NAME
all(row.names(limmadat) == row.names(pheno))

names(limmadat) <- genes

limmadat <- as.data.frame(t(limmadat)) # limma needs samples in columns and genes in rows.

pheno$lateAMD_person <- as.character(pheno$lateAMD_person)
pheno$anyAMD_person <- as.character(pheno$anyAMD_person)
pheno$male <- as.character(pheno$male)
str(pheno)
#fwrite(pheno,'~/Proteomics/Phenotypes.csv')

library(SummarizedExperiment)
eset <- ExpressionSet(assayData = as.matrix(limmadat),
                      phenoData = AnnotatedDataFrame(pheno))
attach(pheno)

# limma, late AMD----
design <- model.matrix(~ lateAMD_person + age + edu + male + smoked + TIME_POINT + RUV)

x <- exprs(eset)
f <- fData(eset)
p <- pData(eset)
fit <- lmFit(eset, design)
fit <- eBayes(fit)
topTable(fit)

summary(fit$s2.post)  # estimated residual variances
sqrt(mean(fit$s2.post))  # 0.7350838

library(RepeatedHighDim)
CIs <- fc_ci(fit, alpha = 0.05)


table(pheno$lateAMD_person) # all timepoints
#0  1 
#91 9 

results <- data.table(colID=row.names(topTable(fit,sort="none",n=Inf, coef=2:8)), topTable(fit,sort="none",n=Inf, coef=2:8))
results <- results[order(adj.P.Val)]

## get limma's logFC for each coefficient----
colnames(design)

logFC_lateAMD_person1 <- as.data.frame(topTable(fit, sort='P', coef=2, n=Inf))
logFC_lateAMD_person1$colID <- row.names(logFC_lateAMD_person1)

results <- merge(results, logFC_lateAMD_person1[,c('logFC','colID')], by='colID')
range(results$logFC, na.rm=T) # -1.0539699  0.9411647

#save.image('~/Proteomics/Proteomics_DIA100_cleanmerge.RData')


# getting confidence intervals-----

CIs <- data.table(colID=row.names(topTable(fit,sort="P",n=Inf, coef=2, confint=T)), topTable(fit,sort="P",n=Inf, coef=2, confint=T))
results <- merge(results, CIs[,c('colID','CI.L','CI.R')], by='colID')
results <- results[order(P.Value)]

## add gene descriptions from refseq----
geneinfo <- fread('~/Proteomics/RefSeq_gene_summaries.csv')
geneinfo <- unique(geneinfo[,c('Symbol','description','Gene_summary')])
geneinfo <- geneinfo[!duplicated(Symbol),]

results <- merge(results, geneinfo[,c('Symbol','description','Gene_summary')], by.x='colID', by.y='Symbol', all.x=T, all.y=F)
results[,absLFC:=abs(lateAMD_person1)]


fwrite(results, '~/Proteomics/Limma/Limma_DIA100/Alltimepts_lateAMD_person_age_male_edu_smoked_RUV_CIs.csv')

# limma, any AMD----
load('~/Proteomics/Proteomics_DIA100_cleanmerge.RData')
fwrite(metadata, '~/Proteomics/Proteomics_DIA100_metadata.csv')
fwrite(limmadat, '~/Proteomics/Proteomics_DIA100_data.csv')

design <- model.matrix(~ anyAMD_person + age + edu + male + smoked + TIME_POINT + RUV)

x <- exprs(eset)
f <- fData(eset)
p <- pData(eset)
fit <- lmFit(eset, design)
fit <- eBayes(fit)
topTable(fit)

summary(fit$s2.post)  # estimated residual variances
sqrt(mean(fit$s2.post)) # 0.734438


table(pheno$anyAMD_person) # all timepoints
#0  1 
#31 69 

results <- data.table(colID=row.names(topTable(fit,sort="none",n=Inf, coef=2:8)), topTable(fit,sort="none",n=Inf, coef=2:8))
results <- results[order(adj.P.Val)]

## get limma's logFC for each coefficient----
colnames(design)

logFC_anyAMD_person1 <- as.data.frame(topTable(fit, sort='P', coef=2, n=Inf))
logFC_anyAMD_person1$colID <- row.names(logFC_anyAMD_person1)

results <- merge(results, logFC_anyAMD_person1[,c('logFC','colID')], by='colID')
range(results$logFC, na.rm=T) # -0.4983209  0.4281593

CIs <- data.table(colID=row.names(topTable(fit,sort="P",n=Inf, coef=2, confint=T)), topTable(fit,sort="P",n=Inf, coef=2, confint=T))
results <- merge(results, CIs[,c('colID','CI.L','CI.R')], by='colID')
results <- results[order(P.Value)]

## add gene descriptions from refseq----
geneinfo <- fread('~/Proteomics/RefSeq_gene_summaries.csv')
geneinfo <- unique(geneinfo[,c('Symbol','description','Gene_summary')])
geneinfo <- geneinfo[!duplicated(Symbol),]

results <- merge(results, geneinfo[,c('Symbol','description','Gene_summary')], by.x='colID', by.y='Symbol', all.x=T, all.y=F)
results[,absLFC:=abs(anyAMD_person1)]

fwrite(results, '~/Proteomics/Limma/Limma_DIA100/Alltimepts_anyAMD_person_age_male_edu_smoked_RUV_CIs.csv')


# pca stuff-----

#annotated PCA after correction


pcadata <- as.data.frame(unique(t(corrected$normalizedCounts)))
names(pcadata) <- prot$Genes

library(PCAtools)
library(factoextra)
library(FactoMineR)

res.pca <- PCA(pcadata, graph = FALSE)
eigenvalues <- res.pca$eig
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")

lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], 
      type="b", pch=19, col = "red")


fviz_screeplot(res.pca, ncp=10)


plot.PCA(res.pca, axes = c(1,2), choix=c("ind"), graph.type = 'ggplot', title='PCA of samples', label = "none",cex=3)
plot.PCA(res.pca, axes = c(2,3), choix=c("var"))


library(ggplot2)
library(factoextra)
library(dplyr)

# Extract PCA coordinates
pca_df <- as.data.frame(res.pca$ind$coord)
pca_df$SampleID <- rownames(pca_df)  # Ensure sample IDs match pheno table
pca_df$RUV <- pheno$RUV
# Merge with phenotype data
pca_df <- merge(pca_df, pheno, by.x = "SampleID", by.y='PARENT_SAMPLE_NAME')

# Plot PCA
ggplot(pca_df, aes(x = Dim.1, y = Dim.2, color = SevScaleAvg)) +
  geom_point(size = 4, alpha = 0.6) +  # Adjust point size and transparency
  scale_color_gradient(low = "blue", high = "red") +  # Continuous color gradient
  theme_bw() + guides(fill=NULL)+
  labs(x = "PC1 (19.39%)", y = "PC2 (5.04%)", color = "Mean severity") +
  tune::coord_obs_pred()+ theme(axis.text = element_text(size=12),
                                legend.position = c(0.15,0.8), 
                                legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5)) 

# 10.88 x 5.25
ggplot(pca_df, aes(x = Dim.1, y = Dim.2, color = batch )) +
  geom_point(size = 4, alpha = 0.6) +  theme_bw() +
  labs(x = "PC1 (19.39%)", y = "PC2 (5.04%)", color = "Batch") +
  tune::coord_obs_pred() + theme(axis.text = element_text(size=12),
                                 legend.position = c(0.1,0.9), 
                                 legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5)) 


pca_df$year <- ifelse(pca_df$TIME_POINT == 10, 5, 0)


ggplot(pca_df, aes(x = Dim.1, y = Dim.2, color = factor(year))) +
  geom_point(size = 4, alpha = 0.6) +  
  theme_bw() +
  labs(x = "PC1 (19.39%)", y = "PC2 (5.04%)", color = "Year") +
  scale_color_manual(values = c("0" = "grey", "5" = "black"), 
                     breaks = c("0", "5")) +  # Ensures only two values appear in the legend
  theme(axis.text = element_text(size=12),
        legend.position = c(0.1, 0.9),  # Move legend inside the plot
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5)) +
  coord_fixed() +  tune::coord_obs_pred()


pca_df$sex <- ifelse(pca_df$male == 1, 'male', 'female')


ggplot(pca_df, aes(x = Dim.1, y = Dim.2, color = factor(sex) )) +
  geom_point(size = 4, alpha = 0.6) +  theme_bw() +
  labs(x = "PC1 (19.39%)", y = "PC2 (5.04%)", color = "Sex") +
  tune::coord_obs_pred() + theme(axis.text = element_text(size=12),
                                 legend.position = c(0.1,0.9), 
                                 legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5)) +
  scale_color_manual(values = c("female" = "grey", "male" = "black"), breaks = c("female", "male")) 


# metadata pca----
library(FactoMineR)
library(factoextra)

# Select only numeric metadata variables
metadata_numeric <- pheno[,c('age','male','batch','lateAMD_person','anyAMD_person',
                             'TIME_POINT','SevScaleAvg','smoked','edu','RUV')]

# Perform PCA
res.pca_meta <- PCA(metadata_numeric, scale.unit = TRUE, graph = FALSE)

# Visualize variance explained by each PC
fviz_eig(res.pca_meta, addlabels = TRUE, ylim = c(0, 100))

fviz_pca_var(res.pca_meta, col.var = "contrib", gradient.cols = c("blue", "red"), repel = TRUE)
