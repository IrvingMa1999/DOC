# Load necessary libraries
library(vegan)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(openxlsx)

# Extract PCoA results from pcoa_points
pcoa_res = pcoa_points[,1:2]

# Read the data and environmental information
data <- read.csv("Seawater_non_zero.csv", row.names = 1)
env <- read.csv("Sample_info_seawater.csv", row.names = 1)

# Convert environmental variables to numeric
T.Chl <- as.numeric(env$T.Chl)
Chla <- as.numeric(env$Chla)
ZD <- as.numeric(env$ZD)
DO <- as.numeric(env$DO)
T <- as.numeric(env$T)
Temp <- as.numeric(env$Temp)

# Perform RDA analysis
rda_result <- rda(t(data) ~ T.Chl + Chla + ZD + DO + T)
summary(rda_result)

# Get the importance of each variable, i.e., each variable's individual contribution to the model
anova_rda_by_var <- anova(rda_result, by = "term", step = 100)
print(anova_rda_by_var)

# Sample scores
rda.sample = data.frame(rda_result$CCA$u[,1:2])
rda.sample$group = env$Treatment # Add group information

# Environmental variable scores
rda.env = data.frame(rda_result$CCA$biplot[,1:2])

# Percentage of variance explained by the first two axes
rda1 = round(rda_result$CCA$eig[1] / sum(rda_result$CCA$eig) * 100, 2)
rda2 = round(rda_result$CCA$eig[2] / sum(rda_result$CCA$eig) * 100, 2)

# Colors for groups
colors <- c("blue", "red")

# Plot the RDA results
ggplot(rda.sample, aes(RDA1, RDA2)) +
  geom_point(aes(fill = group, color = group), size = 3) + 
  scale_fill_manual(values = colors) +
  xlab(paste("RDA1 ( ", rda1, "%", " )", sep = "")) + 
  ylab(paste("RDA2 ( ", rda2, "%", " )", sep = "")) +
  geom_segment(data = rda.env, aes(x = 0, y = 0, xend = rda.env[,1], yend = rda.env[,2]),
               arrow = arrow(length = unit(0.35, "cm"), type = "closed", angle = 22.5),
               linetype = 1, colour = "black", size = 0.6) + # Add arrows
  geom_text_repel(data = rda.env, segment.colour = "black",
                  aes(x = rda.env[,1], y = rda.env[,2], label = rownames(rda.env)), size = 6) +
  # Add dotted lines
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme_bw() +
  theme(axis.title = element_text(family = "serif", face = "bold", size = 18, colour = "black")) +
  theme(axis.text = element_text(family = "serif", face = "bold", size = 16, color = "black")) +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "right") 

# Permutation test
envfit <- envfit(rda_result, env, permutations = 999)
r <- as.matrix(envfit$vectors$r)
p <- as.matrix(envfit$vectors$pvals)
env.p <- cbind(r, p)
colnames(env.p) <- c("r2", "p-value")
KK <- as.data.frame(env.p)
KK$p.adj = p.adjust(KK$`p-value`, method = 'BH')
KK

# RDA species scores extraction
species_scores <- scores(rda_result, display = "species")
head(species_scores)

# Sort species scores
top_species <- species_scores[order(abs(species_scores[, 1]), decreasing = TRUE), ]

# Extract top 100 species scores
top_100_species <- head(top_species, 100)
top_100_species_df <- as.data.frame(top_100_species)

# Save top 100 species scores to a CSV file
write.csv(top_100_species_df, "polar_100_function_df.csv")
print(top_100_species_df)

# Extract and print environmental loadings
env_loadings <- scores(rda_result, choices = 1:2, display = "bp", scaling = 2)
print(env_loadings)