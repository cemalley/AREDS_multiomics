# WGCNA clustering and correlations of modules and metabolites to Mediterranean diet index and its components (AREDS only)


library(data.table);library(WGCNA);library(readxl);library(Hmisc)

setwd('~/Metabolomics/AREDS/Release4/WGCNA/')

data <- fread('~/Metabolomics/AREDS/Release4/Data_with_metadata_prepared_AREDS_release4_person.csv')



visits <- unique(data[,c('PARENT_SAMPLE_NAME','CLIENT_SAMPLE_ID','TIME_POINT')])

table(visits$TIME_POINT)
#  0   2   4   6   8  10  12  13  14  15  16  17  18  19  20  21 
#973 799 850 920 889 873 814   2 670   6 542  15 502   7 227   1 

chemical_names <- fread('~/Metabolomics/AREDS/chemical_names.csv')


xenos <- chemical_names[SUPER_PATHWAY==''| SUPER_PATHWAY=='Xenobiotics' | SUPER_PATHWAY =='Partially Characterized Molecules',CHEM_ID_NEW]
data <- subset(data, select=c(names(data) %nin% xenos))


# start wgcna for all visits----

#########################
library(WGCNA);library(data.table);library(ggthemes);library(ggplot2)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# wgcna--------
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
#sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

datExprM <- subset(data, select=c(grep('CHEM_', names(data))))

sft <- pickSoftThreshold(datExprM, powerVector = powers, verbose=5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Choose the lowest power where R² > 0.85
net = blockwiseModules(datExprM, power = 10,
                       TOMType = "unsigned", minModuleSize = 7,
                       reassignThreshold = 1e-6, mergeCutHeight = 0.5,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "aredsMetaboTOM",
                       verbose = 3)
str(net)
table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
         file = "AREDS_all-networkConstruction.RData")
     
modules <- as.data.frame(moduleLabels)
modules$id <- row.names(modules)
modules <- merge(modules, chemical_names[,c('CHEMICAL_NAME','CHEM_ID_NEW', 'SUPER_PATHWAY', 'SUB_PATHWAY','HMDB','KEGG')], by.x='id', by.y='CHEM_ID_NEW')
modules <- as.data.table(modules)
modules <- modules[order(moduleLabels)]
fwrite(modules, 'AREDS_all_modules.csv')

# add colors-----
colors <- data.table(CHEM_ID_NEW = names(datExprM),
                     moduleColor=labels2colors(net$colors),
                     moduleNumber=net$colors)
colors <- merge(colors, chemical_names[,c('CHEM_ID_NEW','CHEMICAL_NAME','HMDB','KEGG', 'SUPER_PATHWAY','SUB_PATHWAY')],
                by='CHEM_ID_NEW')

fwrite(colors,'WGCNA_module_colors.csv')


# correlation to clinical features------

# ----- Load data -----
# Full WGCNA script with module eigengenes named by color (power = 10)
# -----START: Load and clean data

library(data.table)
library(WGCNA)
library(corrplot)
library(RColorBrewer)
library(Hmisc)           # for rcorr()

# ----- Load and clean data -----
data <- fread("~/Metabolomics/AREDS/Release4/Data_with_metadata_prepared_AREDS_release4_person_including_57219_and_53393.csv")

# Remove xenobiotics/partially characterized compounds
chemical_names <- fread('~/Metabolomics/AREDS/chemical_names.csv')
xenos <- chemical_names[
  SUPER_PATHWAY %in% c("", "Xenobiotics", "Partially Characterized Molecules"),
  CHEM_ID_NEW
]
data <- data[, .SD, .SDcols = !names(data) %in% xenos]

# ----- Build trait matrix -----
diet <- fread('~/AREDS2/AREDS1_Clinical_data/NUTRITION/DATABASES/amedi_areds_allcols.csv')
diet[, REGNO := as.numeric(REGNO)]

datTraitsM <- unique(data[, .(
  CLIENT_SAMPLE_ID,
  PARENT_SAMPLE_NAME,
  TIME_POINT,
  age, male, smoked, edu,
  HighestSevScale, LateAMD_person, BMI
)])
datTraitsM <- merge(
  datTraitsM, diet,
  by.x = "CLIENT_SAMPLE_ID", by.y = "REGNO", all.x = TRUE
)
datTraitsM <- datTraitsM[, .(
  PARENT_SAMPLE_NAME,
  age, male, smoked, BMI, HighestSevScale,
  DT_KCAL, alcohol_gm, alcopts, redmeat, fish, mufasfat,
  whole_fruit, vegetables, whole_grains, nuts, legumes
)]
datTraitsM <- datTraitsM[!duplicated(PARENT_SAMPLE_NAME)]

# ----- Extract only metabolite columns for expression -----
met_cols <- grep("^CHEM_", names(data), value = TRUE)
datExprM <- data[, c("PARENT_SAMPLE_NAME", ..met_cols)]
setkey(datExprM, PARENT_SAMPLE_NAME)
datExprM <- datExprM[!duplicated(PARENT_SAMPLE_NAME)]
datTraitsM <- datTraitsM[!duplicated(PARENT_SAMPLE_NAME)]

# Align samples between expression and traits
datExprM <- datExprM[datTraitsM, on = "PARENT_SAMPLE_NAME"]
expr_matrix <- as.matrix(datExprM[, -1])
rownames(expr_matrix) <- datExprM$PARENT_SAMPLE_NAME

# Final objects
datExprM <- as.data.frame(expr_matrix)
datTraitsM <- as.data.frame(datTraitsM)
rownames(datTraitsM) <- datTraitsM$PARENT_SAMPLE_NAME
datTraitsM$PARENT_SAMPLE_NAME <- NULL
stopifnot(all(rownames(datExprM) == rownames(datTraitsM)))

# ----- WGCNA module detection -----
softPower <- 10

net = blockwiseModules(datExprM, power = softPower,
                       TOMType = "unsigned", minModuleSize = 7,
                       reassignThreshold = 1e-6, mergeCutHeight = 0.5,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "aredsMetaboTOM",
                       verbose = 3)

# Convert numeric module labels to colors
moduleColors <- labels2colors(net$colors)

# Build and rename eigengenes by color
MEs0 <- net$MEs
MEs <- orderMEs(MEs0)
numericIDs <- as.numeric(sub("ME", "", colnames(MEs)))
colorNames <- labels2colors(numericIDs)
colnames(MEs) <- paste0("ME", colorNames)

# ----- Module–trait Spearman correlation (with p-values) -----
common_samples <- intersect(rownames(MEs), rownames(datTraitsM))
MEs <- MEs[common_samples, ]
trait_data <- datTraitsM[common_samples, ]
stopifnot(all(rownames(MEs) == rownames(trait_data)))

# Define module and trait names
module_names <- colnames(MEs)
trait_names  <- colnames(trait_data)

# Remove MEgrey from module names
module_names <- setdiff(module_names, "MEgrey")

# Use rcorr() for Spearman
rc <- rcorr(as.matrix(MEs), as.matrix(trait_data), type = "spearman")

# Extract correlation and p-values between modules and traits
moduleTraitCor    <- rc$r[module_names, trait_names]
moduleTraitPvalue <- rc$P[module_names, trait_names]




# ----- Plot correlation heatmap ----- WORKS-----
myCol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(200)
corLimit <- max(abs(moduleTraitCor), na.rm = TRUE)

corrplot(moduleTraitCor,
         method = "color",hclust.method = 'median',
         col = myCol,
         p.mat = moduleTraitPvalue,
         insig = "label_sig",
         sig.level = c(0.001, 0.01, 0.05),
         pch.col = "grey20",
         pch.cex = 0.9,
         tl.col = "black",
         bg = "lightgrey",
         mar = c(0, 0, 1, 0),
         col.lim = c(-corLimit, corLimit))
# ----- Generate multi-plot PDFs for each non-grey module vs all traits -----
modules <- sub("ME", "", colnames(MEs))
modules <- setdiff(modules, "grey")
traits <- colnames(trait_data)

for (mod in modules) {
  MEcol <- paste0("ME", mod)
  MEvec <- MEs[, MEcol]
  
  n_col <- 4
  n_row <- ceiling(length(traits) / n_col)
  pdf(paste0("Module_", mod, "_trait_scatterplots.pdf"), width = 12, height = 9)
  par(mfrow = c(n_row, n_col), mar = c(4, 4, 2, 1))
  
  for (tr in traits) {
    traitVec <- trait_data[[tr]]
    ok <- complete.cases(MEvec, traitVec)
    
    verboseScatterplot(
      MEvec[ok],
      traitVec[ok],
      xlab = paste("Module Eigengene:", mod),
      ylab = paste("Trait:", tr),
      main = paste0(
        mod, " ~ ", tr,
        "\nr (Spearman) = ", round(moduleTraitCor[MEcol, tr], 2),
        ", p = ", signif(moduleTraitPvalue[MEcol, tr], 2)
      ),
      col = mod,
      cex.main = 0.8
    )
  }
  
  dev.off()
}

# get correlations to individual metabolites------------
myCol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(200)
corLimit <- max(abs(moduleTraitCor), na.rm = TRUE)

corrplot(moduleTraitCor,
         method = "color",hclust.method = 'median',
         col = myCol,
         p.mat = moduleTraitPvalue,
         insig = "label_sig",
         sig.level = c(0.001, 0.01, 0.05),
         pch.col = "grey20",
         pch.cex = 0.9,
         tl.col = "black",
         bg = "lightgrey",
         mar = c(0, 0, 1, 0),
         col.lim = c(-corLimit, corLimit))


# getting spearman correlations for all metabolites to the traits------
# Ensure the same rownames and order
stopifnot(all(rownames(datExprM) == rownames(datTraitsM)))

# Compute Spearman correlation
cor_results <- rcorr(as.matrix(datExprM), as.matrix(datTraitsM), type = "spearman")
save(cor_results, file='Alltimepts_cor_results.RData')

# Extract matrices
r_matrix <- cor_results$r
p_matrix <- cor_results$P

# Extract metabolite-to-trait parts
metabolite_names <- colnames(datExprM)
trait_names <- colnames(datTraitsM)

# Subset correlation and p-values to just metabolite-trait pairs
r_sub <- r_matrix[metabolite_names, trait_names, drop=FALSE]
p_sub <- p_matrix[metabolite_names, trait_names, drop=FALSE]

# Convert to long format data.table
library(data.table)
cor_dt <- as.data.table(as.table(r_sub))
setnames(cor_dt, c("metabolite", "trait", "correlation"))
pval_dt <- as.data.table(as.table(p_sub))
setnames(pval_dt, c("metabolite", "trait", "p_value"))

# Merge correlation and p-value
cor_with_p <- merge(cor_dt, pval_dt, by = c("metabolite", "trait"))

# Add chemical names if desired
cor_with_p <- merge(cor_with_p, chemical_names[, .(CHEM_ID_NEW, CHEMICAL_NAME, SUPER_PATHWAY, SUB_PATHWAY)],
                    by.x = "metabolite", by.y = "CHEM_ID_NEW", all.x = TRUE)

# Save
fwrite(cor_with_p, "alltimepts_metabolite_trait_spearman_correlations.csv")
cor_with_p <- fread('alltimepts_metabolite_trait_spearman_correlations.csv')
# keep only aMedi trait correlations.-----

cor_with_p <- cor_with_p[!is.na(CHEMICAL_NAME) & CHEMICAL_NAME !='' & trait %in% c('DT_KCAL','alcohol_gm','alcopts','redmeat','fish',
                                                              'mufasfat','whole_fruit','vegetables','whole_grains',
                                                              'nuts','legumes'),]

library(ggplot2)
library(data.table)

library(corrplot)
library(RColorBrewer)

# Filter to relevant aMedi traits
amedi_traits <- c('DT_KCAL','alcohol_gm','alcopts','redmeat','fish',
                  'mufasfat','whole_fruit','vegetables','whole_grains',
                  'nuts','legumes')

cor_subset <- cor_with_p[!is.na(CHEMICAL_NAME)& CHEMICAL_NAME !='' & trait %in% amedi_traits]

# Get top 10 metabolites by max |correlation|
top_mets <- cor_subset[, .(max_abs_cor = max(abs(correlation))), by = .(metabolite, CHEMICAL_NAME)]
top10 <- top_mets[order(-max_abs_cor)][1:10]

# Subset correlations and p-values to top 10 metabolites and aMedi traits
sub <- cor_subset[metabolite %in% top10$metabolite]

# Prepare correlation and p-value matrices
cor_mat <- dcast(sub, CHEMICAL_NAME ~ trait, value.var = "correlation")
pval_mat <- dcast(sub, CHEMICAL_NAME ~ trait, value.var = "p_value")

# Set row names and remove first column
cor_mat <- as.data.frame(cor_mat)
rownames(cor_mat) <- cor_mat$CHEMICAL_NAME
cor_mat$CHEMICAL_NAME <- NULL

pval_mat <- as.data.frame(pval_mat)
rownames(pval_mat) <- pval_mat$CHEMICAL_NAME
pval_mat$CHEMICAL_NAME <- NULL

# Convert to matrices
cor_mat <- as.matrix(cor_mat)
pval_mat <- as.matrix(pval_mat)

# Plot
myCol <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(200)
corLimit <- max(abs(cor_mat), na.rm = TRUE)

corrplot(cor_mat,
         method = "color",
         hclust.method = 'median',
         col = myCol,
         p.mat = pval_mat,
         insig = "label_sig",
         sig.level = c(0.001, 0.01, 0.05),
         pch.col = "grey20",
         pch.cex = 0.9,
         tl.col = "black",
         bg = "lightgrey",
         mar = c(0, 0, 1, 0),
         col.lim = c(-corLimit, corLimit))
