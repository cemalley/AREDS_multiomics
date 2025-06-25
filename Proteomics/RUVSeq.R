# get genes that are low covariance for ruvseq-----
setwd('~/Proteomics/Limma/Limma_DIA100_prenorm/')
b0 <- fread('Baseline_amd2_person_age_male_batch.csv')
b5 <- fread('Year5_amd2_person_age_male_batch.csv')

b0[adj.P.Val < 0.05,] # 454
# find significant, relevant genes.
b0[abs(amd2_person1) >=0.01 &  abs(batchnew) <= 0.075 & P.Value < 0.05,] # 0.075 is the 1st and 3rd quartiles of the batch model estimates
b0[abs(batchnew) <= 0.075 & adj.P.Val > 0.05,] # still 804 genes.
b0[abs(amd2_person1) >=0.01 & adj.P.Val < 0.05,] 

b0[(abs(amd2_person1) < abs(batchnew)) & P.Value > 0.05,]

b5[(abs(amd2_person1) < abs(batchnew) | (abs(amd2_person2) < abs(batchnew))) & P.Value > 0.05,] 

spikes <- unique(intersect(b0[(abs(amd2_person1) < abs(batchnew)) & P.Value > 0.05,colID] ,
                           b5[(abs(amd2_person1) < abs(batchnew) | (abs(amd2_person2) < abs(batchnew))) & P.Value > 0.05,colID]  ))

length(spikes) # 418.
fwrite(data.table(GeneID=spikes), '~/Proteomics/Spikes_list.csv')
#



# start ruvseq----
setwd('~/Proteomics/')
library(data.table);library(preprocessCore);library(impute);library(readxl);library(zoo);library(Biobase);library(limma);library(tidyverse)

metadata <- fread('Metadata.csv')

prot <- fread('report.pg_matrix_DIA_100.tsv')

prot[Genes=="",Genes:=tstrsplit(Protein.Names, '_HUMAN')[1]]
prot[Genes=="IGG1",Genes:='IGHG1_1']
prot[Genes=="IGA2",Genes:='IGHA2_1']
prot[Genes=="IGD",Genes:='IGHD_1']
prot[Genes=="IGE",Genes:='IGHE_1']


norm <- as.matrix(prot[,6:105])

imput <- impute.knn(norm) # 507 rows with more than 50 % 
imput <- imput$data

imput <- as.data.table(imput)
names(imput) <- gsub('_depleteDIA_3uL.raw','',names(imput))
names(imput) <- gsub('C:\\\\Users\\\\tharakanr2\\\\Documents\\\\NEI_Collaboration\\\\DIA\\\\20231207_NEI_DIA\\\\','', names(imput))
names(imput) <- gsub('C:\\\\Users\\\\tharakanr2\\\\Documents\\\\NEI_Collaboration\\\\DIA\\\\20231002_AMD_DIA\\\\','', names(imput))

names(imput) <- gsub('NEYI00','NEYI-00', names(imput))
names(imput)[names(imput) %in% metadata$PARENT_SAMPLE_NAME]
imput$colID <- prot$Genes
imput[,colID:=tstrsplit(colID, ';')[1]]

imput.t <- as.data.table(t(imput[,-c('colID')]))
names(imput.t) <- imput$colID
imput.t$PARENT_SAMPLE_NAME <- names(imput[,-c('colID')])

metadata <- metadata[PARENT_SAMPLE_NAME %in% names(imput),]

imput.t <- imput.t[PARENT_SAMPLE_NAME !='NEYI-00877',] # far outlier in batch 2

names <- names(imput.t[,-c(1683)])
names <- gsub('-', '_', names)

names(imput.t)[1:1682] <- names

###

assayData <-imput.t[,-c(1683)]
assayData <- as.data.table(t(assayData))
#names(assayData) <- imput.t$PARENT_SAMPLE_NAME
assayData <- assayData %>% mutate_if(is.numeric, round)
names(assayData) <- imput.t$PARENT_SAMPLE_NAME
assayData <- as.matrix(assayData)
rownames(assayData) <- names(imput.t[,-c(1683)])

spikes <- fread('~/Proteomics/Spikes_list.csv') # 418



corrected <- RUVg(assayData, spikes, k=1)

# before correction:
plotRLE(assayData, outline=FALSE, ylim=c(-1, 1))

# after correction:
plotRLE(corrected$normalizedCounts, outline=FALSE, ylim=c(-1, 1))

plotPCA(corrected$normalizedCounts, cex=1.2)

# pca----
pcadata <- t(corrected$normalizedCounts)

res.pca <- PCA(pcadata, graph = FALSE)
eigenvalues <- res.pca$eig
barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], 
      type="b", pch=19, col = "red")


fviz_screeplot(res.pca, ncp=10)


plot.PCA(res.pca, axes = c(1,2), choix=c("ind"))


plot(res.pca, choix = "ind")

pcadata <- as.data.frame(res.pca$ind$coord)

pcadata$PARENT_SAMPLE_NAME <- row.names(pcadata)

orig <- #XXXX
new <- unique(pcadata$PARENT_SAMPLE_NAME)[unique(pcadata$PARENT_SAMPLE_NAME) %nin% orig]

pcadata <- as.data.table(pcadata)

pcadata[,batch:=ifelse(PARENT_SAMPLE_NAME %in% new, 'new', 'original')]

ggplot(pcadata, aes(x=Dim.1, y=Dim.2, color=batch)) + 
  geom_point(size=2)+theme_clean() +
  labs(x='Dim.1: 8% variance explained',
       y='Dim.2: 6% variance explained')

save(corrected, file='RUVSeq_corrections_DIA100_preSpectronaut.RData')

imput.t <- fread('~/Proteomics/imput.t.csv')

pcadata <- merge(pcadata, imput.t[,1:7], by='PARENT_SAMPLE_NAME')

pcadata$TIME_POINT <- as.character(pcadata$TIME_POINT)
pcadata$male <- as.character(pcadata$male)
pcadata$anyAMD_person <- as.character(pcadata$anyAMD_person)
pcadata$amd2_person <- as.character(pcadata$amd2_person)

ggplot(pcadata, aes(x=Dim.1, y=Dim.2, color=amd2_person)) + 
  geom_point(size=2)+theme_clean() +
  labs(x='Dim.1: 8% variance explained',
       y='Dim.2: 6% variance explained')



#
