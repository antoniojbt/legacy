#source("http://bioconductor.org/biocLite.R")
#biocLite()
#install.packages("dendextend")
library(ggplot2)
library(flashClust)
library(dendextend)


# PCA
plot_normalised_filtered_by_targets <- ("MDS_normalised_filtered_by_targets.png")
png(plot_normalised_filtered_by_targets, width = 4, height = 4, units = 'in', res = 300)
plotMDS(normalised_expressed, pch=1, labels=membership_file_cleaned$treatment)
plotMDS(normalised_expressed, pch=1, labels=membership_file_cleaned$cell_type)
plotMDS(normalised_expressed, pch=1, labels=membership_file_cleaned$individual)
dev.off()


# Plot PCA of normalised samples:
# Compute the PCs, first transpose the expression values and ignore the first column, then run the PCs:
pca_normalised_filtered <- prcomp(t(normalised_expressed$E), center=TRUE, scale=TRUE)

# Obtain values for all PCs output:
pc <- data.frame(round(pca_normalised_filtered$x, 2))
pc$ID <- rownames(pc)
source('/Users/antoniob/Desktop/ifs_projects_antoniob.dir/VD_1ry_cells.dir/scripts.dir/moveme.R')
pc <- pc[, moveme(names(pc), 'ID first')]
names(pc)[1:10]
class(pc)
write.table(pc, 'principal_components_normalised_filtered.tsv', quote = FALSE,
            sep = '\t', row.names = FALSE)

# Explore dimensions and plot first 10 or so components:
dim(pc)
dim(normalised_expressed)
str(pca_normalised_filtered)
head(pc)
# pc[1:5, 1:5]
summary(pca_normalised_filtered)


# Plot PCA results:
plot_PCA_normalised_filtered <- ('plot_PCA_normalised_filtered.png')
png(plot_PCA_normalised_filtered, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
# Histogram of first x PCs:
plot(pca_normalised_filtered, main='Normalised expression values')
# Scatterplot of PC1 and PC2:
biplot(pca_normalised_filtered, main='Normalised expression values')

par(mfrow=c(1,1))
dev.off()


# Run PCA analysis by groups of interest: TO DO: Cross files and IDs first
head(membership_file_cleaned)
tail(membership_file_cleaned)
pc_data <- data.frame(pca_normalised_filtered$x[, 1:13])
str(pc_data)
head(pc_data)

pca_by_groups <- data.frame(merge(membership_file_cleaned, pc_data, by = 'row.names'))
head(pca_by_groups)
dim(pca_by_groups)
dim(pc_data)
# View(pc_data)
# View(pca_by_groups)
head(arrange(pc_data, PC1), 10)
head(arrange(pca_by_groups, PC1), 10)

plot_PCA_normalised_filtered_by_groups_1 <- ('plot_PCA_normalised_filtered_by_groups_1.png')
png(plot_PCA_normalised_filtered_by_groups_1, width = 13, height = 13, units = 'in', res = 300)
p1 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$individual)) + theme(legend.position="bottom")
p2 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p3 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$timepoint)) + theme(legend.position="bottom")
p4 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$cell_type)) + theme(legend.position="bottom")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

p5 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$replicate)) + theme(legend.position="bottom")
p6 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$sentrix_ID)) + theme(legend.position="bottom")


# Timepoints:
qplot(x=PC2, y=PC3, data=pca_by_groups, colour=factor(membership_file_cleaned$timepoint)) + theme(legend.position="bottom")
qplot(x=PC3, y=PC4, data=pca_by_groups, colour=factor(membership_file_cleaned$timepoint)) + theme(legend.position="bottom")
qplot(x=PC4, y=PC5, data=pca_by_groups, colour=factor(membership_file_cleaned$timepoint)) + theme(legend.position="bottom")
qplot(x=PC5, y=PC6, data=pca_by_groups, colour=factor(membership_file_cleaned$timepoint)) + theme(legend.position="bottom")
qplot(x=PC7, y=PC8, data=pca_by_groups, colour=factor(membership_file_cleaned$timepoint)) + theme(legend.position="bottom")
qplot(x=PC8, y=PC9, data=pca_by_groups, colour=factor(membership_file_cleaned$timepoint)) + theme(legend.position="bottom")


# Treatment:
qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
qplot(x=PC2, y=PC3, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
qplot(x=PC3, y=PC4, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
qplot(x=PC4, y=PC5, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
qplot(x=PC5, y=PC6, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")

# Plot dendrogram
# Prepare hierarchical cluster:
correlation <- cor(normalised_expressed$E, method='pearson')
head(correlation)
# View(correlation)

hc <- flashClust(dist(correlation))
str(hc)
head(hc)

# Convert to a dendogram object:
dendro <- as.dendrogram(hc)

# Add colours by groups
colourCodes <- c(treated="red", basal="green")
#colourCodes <- c(CD14="red", CD4="green", CD8="blue", CD19="yellow")

# Assign labels of dendrogram object with new colors:
labels_colors(dendro) <- colourCodes[membership_file_cleaned$treatment][order.dendrogram(dendro)]
#labels_colors(dendro) <- colourCodes[membership_file_cleaned$cell_type][order.dendrogram(dendro)]
head(colourCodes[membership_file_cleaned$treatment][order.dendrogram(dendro)])
#head(colourCodes[membership_file_cleaned$cell_type][order.dendrogram(dendro)])
head(order.dendrogram(dendro))
head(membership_file_cleaned$treatment)
head(membership_file_cleaned$cell_type)

# Plot simple dendrogram, labels at the same level:
plot_dendrogram_normalised_filtered <- ('plot_dendrogram_normalised_filtered.png')
png(plot_dendrogram_normalised_filtered, width = 13, height = 13, units = 'in', res = 300)
par(mfrow=c(1,2), cex = 1)
plot(dendro, hang = -1)
# Zoom in to a sub tree:
plot(dendro[[1]], horiz = TRUE)
par(mfrow=c(1,1))
dev.off()


# Interpretation: There is a strong separation by cell type but not by treatment. There isn't a batch effect or differences between
# replicates or individuals.