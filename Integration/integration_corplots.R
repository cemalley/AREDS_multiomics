# integration of AREDS proteomics and metabolomics and differential correlation to late AMD

library(data.table)
library(ggplot2)
library(corrplot)
library(RColorBrewer)
library(circlize)
library(Hmisc)
library(ComplexHeatmap)
library(cluster)

setwd('~/Integration/')

# --- Load data ---
prot <- fread('AREDS_proteomics.csv')
mets <- fread('AREDS_metabolomics.csv')

# --- Chemical annotations ---
chemical.names <- fread('~/Metabolomics/AREDS/chemical_names.csv')
endobiotics <- chemical.names[
  SUPER_PATHWAY %nin% c('Xenobiotics', 'Partially Characterized Molecules') & TYPE == 'NAMED',
  CHEM_ID_NEW
]
keepcols <- c('PARENT_SAMPLE_NAME', 'LateAMD_person', endobiotics)
mets <- subset(mets, select = keepcols)

# --- Match samples ---
common_ids <- intersect(prot$PARENT_SAMPLE_NAME, mets$PARENT_SAMPLE_NAME)
prot <- prot[PARENT_SAMPLE_NAME %in% common_ids][order(PARENT_SAMPLE_NAME)]
mets <- mets[PARENT_SAMPLE_NAME %in% common_ids][order(PARENT_SAMPLE_NAME)]

# --- Merge and subset ---
sample_info <- prot[, .(PARENT_SAMPLE_NAME, lateAMD_person)]
merged <- merge(merge(prot, mets, by = "PARENT_SAMPLE_NAME"), sample_info, by = "PARENT_SAMPLE_NAME")

prot_cols <- names(merged)[14:1083]
met_cols  <- names(merged)[1085:1978]
met_cols <- met_cols[met_cols %nin% c('lateAMD_person.x', 'lateAMD_person.y', 'CHEM_100001316')]
prot_data <- as.matrix(merged[, ..prot_cols])
met_data  <- as.matrix(merged[, ..met_cols])
status    <- merged$LateAMD_person

# --- Remove zero-variance proteins ---
sd_vals <- apply(prot_data, 2, sd, na.rm = TRUE)
nonzero_var <- !is.na(sd_vals) & sd_vals > 0
prot_data <- prot_data[, nonzero_var]

# --- Remove rows with too many NAs ---
good_rows <- rowSums(is.na(prot_data)) < ncol(prot_data)
prot_data <- prot_data[good_rows, ]
met_data  <- met_data[good_rows, ]
status    <- status[good_rows]

# --- Impute missing values with column means ---
prot_imputed <- apply(prot_data, 2, function(x) {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
})
prot_imputed <- as.matrix(prot_imputed)

# --- Correlations by AMD status ---
prot_case <- prot_imputed[status == 1, ]
met_case  <- met_data[status == 1, ]
prot_ctrl <- prot_imputed[status == 0, ]
met_ctrl  <- met_data[status == 0, ]

cor_case  <- cor(met_case, prot_case, use = "pairwise.complete.obs")
cor_ctrl  <- cor(met_ctrl, prot_ctrl, use = "pairwise.complete.obs")
cor_diff  <- cor_case - cor_ctrl

# --- Fisher Z-transform and difference ---
fisher_z <- function(r, n) 0.5 * log((1 + r) / (1 - r)) * sqrt(n - 3)
z_case <- fisher_z(cor_case, nrow(met_case))
z_ctrl <- fisher_z(cor_ctrl, nrow(met_ctrl))
z_diff <- z_case - z_ctrl
p_diff <- 2 * pnorm(-abs(z_diff))  # two-sided p-value

# --- Filter by both correlation difference and p-value ---
mask <- abs(cor_diff) > 0.3 & p_diff < 0.05
cor_diff_masked <- cor_diff
cor_diff_masked[!mask] <- NA

# --- Remove any rows or columns with NA ---
cor_diff_filtered <- cor_diff_masked
cor_diff_filtered <- cor_diff_filtered[complete.cases(cor_diff_filtered), ]
cor_diff_filtered <- cor_diff_filtered[, colSums(is.na(cor_diff_filtered)) == 0]

# --- Final check ---
stopifnot(!any(is.na(cor_diff_filtered)))


# --- Clustering ---
dist_rows <- dist(cor_diff_filtered, method = "euclidean")
dist_cols <- dist(t(cor_diff_filtered), method = "euclidean")

row_clust <- hclust(dist_rows, method = "ward.D2")
col_clust <- hclust(dist_cols, method = "ward.D2")

# --- Define number of clusters manually or by silhouette separately ---
k_rows <- 7
k_cols <- 5

row_groups <- cutree(row_clust, k = k_rows)
col_groups <- cutree(col_clust, k = k_cols)

# --- Annotations ---
row_anno <- rowAnnotation(cluster = as.factor(row_groups))
col_anno <- HeatmapAnnotation(cluster = as.factor(col_groups))

# --- Plot ---
Heatmap(cor_diff_filtered,
        name = "Δ correlation",
        col = colorRamp2(c(-1, 0, 1), rev(brewer.pal(3, "RdBu"))),
        cluster_rows = row_clust,
        cluster_columns = col_clust,
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        left_annotation = row_anno,
        top_annotation = col_anno,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(title = "Case - Control"))

#---------------------------------------------------
# --- Rank and subset top gene–metabolite correlation differences ---

# Flatten correlation difference matrix into a data.table
cor_diff_long <- as.data.table(as.table(cor_diff))
setnames(cor_diff_long, c("Metabolite", "Gene", "CorDiff"))

# Add p-values and rank
cor_diff_long[, AbsCorDiff := abs(CorDiff)]
cor_diff_long[, Pval := as.vector(p_diff)]

# Filter: only finite and non-missing
cor_diff_long <- cor_diff_long[is.finite(CorDiff) & !is.na(CorDiff)]

# Select top N by absolute correlation difference
top_n <- 200
top_pairs <- cor_diff_long[order(-AbsCorDiff)][1:top_n]

# Subset matrix to top rows and columns
top_mets <- unique(top_pairs$Metabolite)
top_genes <- unique(top_pairs$Gene)

cor_diff_filtered <- cor_diff[top_mets, top_genes, drop = FALSE]


dist_rows <- dist(cor_diff_filtered, method = "euclidean")
dist_cols <- dist(t(cor_diff_filtered), method = "euclidean")

row_clust <- hclust(dist_rows, method = "ward.D2")
col_clust <- hclust(dist_cols, method = "ward.D2")

# --- Define number of clusters manually or by silhouette separately ---
k_rows <- 7
k_cols <- 5

row_groups <- cutree(row_clust, k = k_rows)
col_groups <- cutree(col_clust, k = k_cols)

# --- Annotations ---
row_anno <- rowAnnotation(cluster = as.factor(row_groups))
col_anno <- HeatmapAnnotation(cluster = as.factor(col_groups))

# --- Plot ---
Heatmap(cor_diff_filtered,
        name = "Δ correlation",
        cluster_rows = row_clust,
        cluster_columns = col_clust,
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        left_annotation = row_anno,
        top_annotation = col_anno,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(title = "Case - Control"))

# add row annotations----
library(ComplexHeatmap)

meta_anno <- chemical.names[CHEM_ID_NEW %in% rownames(cor_diff_filtered),] # only 135 mets


# Ensure rownames of meta_anno match those in cor_diff_filtered
meta_anno_df <- as.data.frame(meta_anno)
rownames(meta_anno_df) <- meta_anno_df$CHEM_ID_NEW

# Subset annotation to match rows of the heatmap
meta_anno_df <- meta_anno_df[rownames(cor_diff_filtered), ]

# Ensure cluster IDs match
meta_anno_df$Cluster <- factor(full_row_groups[rownames(meta_anno_df)])

# Define colors
super_colors <- setNames(brewer.pal(length(unique(meta_anno_df$SUPER_PATHWAY)), "Set3"),
                         unique(meta_anno_df$SUPER_PATHWAY))

cluster_colors <- setNames(brewer.pal(length(unique(meta_anno_df$Cluster)), "Paired"),
                           levels(meta_anno_df$Cluster))

# Combined row annotation
row_anno <- rowAnnotation(
  Cluster = meta_anno_df$Cluster,
  SUPER_PATHWAY = meta_anno_df$SUPER_PATHWAY,
  col = list(
    Cluster = cluster_colors,
    SUPER_PATHWAY = super_colors
  ),
  annotation_name_side = "top",
  annotation_legend_param = list(Cluster = list(title = "Cluster"),
                                 SUPER_PATHWAY = list(title = "Super Pathway"))
)

# Map CHEM_ID_NEW to chemical names
chem_name_map <- setNames(meta_anno$CHEMICAL_NAME, meta_anno$CHEM_ID_NEW)

# Apply renaming to cor_diff_filtered rownames
rownames(cor_diff_filtered) <- chem_name_map[rownames(cor_diff_filtered)]

# Also update annotation rownames to match
rownames(meta_anno_df) <- chem_name_map[rownames(meta_anno_df)]

Heatmap(cor_diff_filtered,
        name = "Δ correlation",
        cluster_rows = row_clust,
        cluster_columns = col_clust,
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        left_annotation = row_anno,
        top_annotation = col_anno,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(title = "Case - Control"))

#fixing row cluster annotations----
# Extract the correct row order from clustering
ordered_row_names <- rownames(cor_diff_filtered)[row_clust$order]
ordered_col_names <- colnames(cor_diff_filtered)[col_clust$order]

# Apply cutree on full hclust object
row_clusters <- cutree(row_clust, k = k_rows)
col_clusters <- cutree(col_clust, k = k_cols)

# Reorder to match heatmap row order
row_clusters_ordered <- row_clusters[ordered_row_names]
col_clusters_ordered <- col_clusters[ordered_col_names]

# --- Row annotations: Metabolite pathways ---
# Make sure Cluster is a factor with defined levels
row_clusters_ordered <- factor(row_clusters_ordered, levels = as.character(1:k_rows))
meta_anno_df$SUPER_PATHWAY <- factor(meta_anno_df$SUPER_PATHWAY)

# Define color palettes
cluster_colors <- setNames(brewer.pal(k_rows, "Set1"), levels(row_clusters_ordered))
super_colors <- setNames(brewer.pal(length(levels(meta_anno_df$SUPER_PATHWAY)), "Set3"),
                         levels(meta_anno_df$SUPER_PATHWAY))

# Create row annotation with both cluster and super pathway
row_ha <- rowAnnotation(
  Cluster = row_clusters_ordered,
  SUPER_PATHWAY = meta_anno_df$SUPER_PATHWAY,
  col = list(
    Cluster = cluster_colors,
    SUPER_PATHWAY = super_colors
  ),
  annotation_name_side = "top"
)

# --- Final heatmap ---
Heatmap(cor_diff_filtered,
        name = "Δ correlation",
        col = colorRamp2(c(-1, 0, 1), rev(brewer.pal(3, "RdBu"))),
        cluster_rows = row_clust,
        cluster_columns = col_clust,
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        left_annotation = row_ha,
        top_annotation = col_ha,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(title = "Case - Control"))



# table for the pairs-----
# Get intersecting row and column names that exist in cor_case and cor_ctrl
common_mets   <- intersect(rownames(cor_case), names(rownames(cor_diff_filtered)))
common_prots  <- intersect(colnames(cor_case), colnames(cor_diff_filtered))

# Subset to those common ones
cor_case_sub  <- cor_case[common_mets, common_prots, drop = FALSE]
cor_ctrl_sub  <- cor_ctrl[common_mets, common_prots, drop = FALSE]
cor_diff_sub  <- cor_diff[common_mets, common_prots, drop = FALSE]
p_diff_sub    <- p_diff[common_mets, common_prots, drop = FALSE]

# Melt to long format
cor_diff_long <- as.data.table(as.table(cor_diff_sub))
setnames(cor_diff_long, c("Metabolite", "Protein", "cor_diff"))

p_diff_long <- as.data.table(as.table(p_diff_sub))
setnames(p_diff_long, c("Metabolite", "Protein", "p_diff"))

cor_case_long <- as.data.table(as.table(cor_case_sub))
setnames(cor_case_long, c("Metabolite", "Protein", "cor_case"))

cor_ctrl_long <- as.data.table(as.table(cor_ctrl_sub))
setnames(cor_ctrl_long, c("Metabolite", "Protein", "cor_ctrl"))

# Merge all pieces
summary_dt <- merge(cor_diff_long, p_diff_long, by = c("Metabolite", "Protein"))
summary_dt <- merge(summary_dt, cor_case_long, by = c("Metabolite", "Protein"))
summary_dt <- merge(summary_dt, cor_ctrl_long, by = c("Metabolite", "Protein"))

# Add annotations
summary_dt <- merge(summary_dt,
                    meta_anno[, .(CHEM_ID_NEW, CHEMICAL_NAME, SUPER_PATHWAY, SUB_PATHWAY)],
                    by.x = "Metabolite", by.y = "CHEM_ID_NEW", all.x = TRUE)

# Add cluster info if available
if (exists("full_row_groups")) {
  summary_dt[, Cluster := full_row_groups[as.character(Metabolite)]]
}

# Rank and limit
summary_dt[, abs_diff := abs(cor_diff)]
setorder(summary_dt, -abs_diff)
top_200_summary <- summary_dt[1:200]

# Clean column order
top_200_summary <- top_200_summary[, .(
  Metabolite, CHEMICAL_NAME, Protein,
  cor_case, cor_ctrl, cor_diff, p_diff, abs_diff,
  Cluster, SUPER_PATHWAY, SUB_PATHWAY
)]

# View or write
head(top_200_summary)
fwrite(top_200_summary, "top200_metabolite_protein_pairs.csv")



# get consistently high diff pairs------
# Set a threshold for meaningful correlation difference
threshold <- 0.3

# Get absolute differences
abs_diff <- abs(cor_diff_filtered)

# Count how many strong differences each metabolite has (row-wise)
met_strong_counts <- rowSums(abs_diff > threshold, na.rm = TRUE)

# Count how many strong differences each protein has (column-wise)
prot_strong_counts <- colSums(abs_diff > threshold, na.rm = TRUE)

# Also compute average absolute difference per metabolite/protein
met_mean_diff <- rowMeans(abs_diff, na.rm = TRUE)
prot_mean_diff <- colMeans(abs_diff, na.rm = TRUE)

# Combine into data.tables
met_summary <- data.table(
  CHEM_ID_NEW = rownames(cor_diff_filtered),
  n_strong_pairs = met_strong_counts,
  mean_abs_diff = met_mean_diff
)

prot_summary <- data.table(
  PROTEIN = colnames(cor_diff_filtered),
  n_strong_pairs = prot_strong_counts,
  mean_abs_diff = prot_mean_diff
)

# Sort to find top metabolites and proteins
setorder(met_summary, -n_strong_pairs)
setorder(prot_summary, -n_strong_pairs)

# View top 10
head(met_summary, 10)
head(prot_summary, 10)


fwrite(prot_summary, 'top_prots_summary.csv')
fwrite(met_summary, 'top_met_summary.csv')

