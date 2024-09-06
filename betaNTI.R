# Load necessary libraries and install dependencies
library(BiocManager)
BiocManager::install("NST")
library(NST)
library(picante)
library(ape)
library(ggplot2)
library(randomForest)
library(tidyr)
library(RColorBrewer)

# Read in the OTU table and transpose it
otu <- read.csv('otu_seawater_rare.csv', row.names = 1)
otu <- data.frame(t(otu))

# Read in sample group information
group <- read.csv('seawater_sample_rare.csv', row.names = 1)

# Read and prune the phylogenetic tree
tree <- read.tree('unrooted.nwk')
tree$tip.label <- gsub("'", "", tree$tip.label)
tree <- prune.sample(otu, tree)

# Calculate phylogenetic community composition metrics
set.seed(123)
pnst <- pNST(comm = otu, tree = tree, group = group, phylo.shuffle = TRUE, 
             pd.wd = tempdir(), abundance.weighted = TRUE, rand = 999, 
             nworker = 80, SES = TRUE, RC = FALSE)

# Extract betaMNTD and betaNTI values
betaMNTD <- pnst$index.pair
head(betaMNTD)

# Save betaMNTD and betaNTI results to a CSV file
write.csv(betaMNTD, 'betaNTI.csv', quote = FALSE, row.names = FALSE)

# Plotting betaNTI values using ggplot2
df <- read.csv("betaNTI_polar_huigui.csv", header = TRUE)
head(df)

ggplot(df, aes(x = DOC_cha, y = bNTI.wt)) +
  geom_point()


# Scatter plot of betaNTI values with custom color scheme
df <- read.csv("betaNTI.csv", header = TRUE)
color_vector <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#85442E")
head(df)

sp1 <- ggplot(df, aes(x = group, y = betaNTI, color = group)) + 
  geom_point(position = "jitter", alpha = 0.7, size = 2) + 
  geom_hline(yintercept = 2, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(-7.5, 7.5))+
  geom_hline(yintercept = -2, linetype = "dashed", color = "black") + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_color_manual(values = color_vector)

print(sp1)
ggsave("betaNTI_point.pdf", sp1, width = 8, height = 6)