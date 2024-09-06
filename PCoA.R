# Clear the current workspace
rm(list = ls())

# Set the working directory
setwd("D:/20240903_16S/")

# Load necessary libraries
library(ggplot2)
library(openxlsx)
library(vegan)

# Read the data
data <- read.csv("otu.csv", row.names = 1)
group <- read.csv("Sample_info.csv", row.names = 1)

# Calculate Bray-Curtis dissimilarity matrix
dist = vegdist(t(data), method = "bray", diag = TRUE, upper = TRUE)

# Convert dist to matrix
dist <- as.matrix(dist)

# Perform classical multidimensional scaling (PCoA)
pcoa = cmdscale(dist, eig = TRUE)
eig = summary(eigenvals(pcoa))
head(eig)

# Set names for each dimension, starting from PCoA1
axis = paste0("PCoA", 1:ncol(eig))
axis

# Axis explained variance
eig = data.frame(Axis = axis, t(eig)[, -3])
head(eig)

# Get the explained variance of the first two dimensions and round to two decimal places
pco1 = round(eig[1, 3] * 100, 2)
pco2 = round(eig[2, 3] * 100, 2)

# Set x-axis and y-axis titles for the plot
xlab = paste0("PCoA1 (", pco1, "%)")
ylab = paste0("PCoA2 (", pco2, "%)")

# Get the coordinates of each sample in the first two dimensions
pcoa_points = as.data.frame(pcoa$points)
pcoa_points

# Merge sample coordinates with corresponding group information
pcoa_points = data.frame(pcoa_points, group = group[4])

# Conduct PERMANOVA analysis
permanova_results <- adonis(dist ~ Group, data = pcoa_points)

# Format PERMANOVA results for the plot
signif_text <- paste("PERMANOVA: R² =",
                     format(paste(round(permanova_results$aov.tab$R2[1], digits = 4)), nsmall = 4),
                     ", p =",
                     format(paste(permanova_results$aov.tab$'Pr(>F)'[1]), nsmall = 4))

# Plot the data
p = ggplot(pcoa_points, aes(V1, V2)) +
  geom_point(size = 4, aes(color = Group)) +     # Draw points by group with different colors and shapes
  labs(x = xlab, y = ylab, subtitle = signif_text) + # Add x, y axis titles and plot subtitle
  stat_ellipse(aes(fill = Group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + # Add confidence ellipse
  theme(plot.title = element_text(size = 12, hjust = 0.5)) +
  geom_hline(yintercept = 0, linetype = 4) +        # Add vertical auxiliary line at origin
  geom_vline(xintercept = 0, linetype = 4) +  # Add horizontal auxiliary line at origin
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Display the plot
p