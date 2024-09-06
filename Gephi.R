##########################################
###### Generate a beautified template for Gephi #########
##########################################

# Set the working directory
setwd("D:/16S_NEW/")

# Clear the workspace
rm(list = ls())

# Load necessary packages
library(rgexf)
library(microeco)

# Load the data
#otu = read.csv("otu.csv", header = T, row.names = 1)
#tax = read.csv("Taxonomy_table.csv", row.names = 1)
dataset <- microtable$new(otu_table = polar_High
                          #,tax_table = tax
                          )
print(dataset)
dataset$tidy_dataset()

# Create a trans_network object and calculate the network
t1 <- trans_network$new(dataset = dataset, cor_method = "spearman", filter_thres = 0.001)
t1$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)

# Calculate modularity
t1$cal_module()

# Calculate node table
t1$get_node_table(node_roles = )
zipi <- t1$res_node_table
write.csv(zipi, "tropical_low_function_zipi.csv")

# Plot the ZIPI graph
t1$plot_taxa_roles(use_type = 1)

# Calculate topological parameters of the network graph
t1$cal_network_attr()
canshu <- t1$res_network_attr

# Get the adjacency matrix
t1$get_adjacency_matrix()

# Plot the network graph
t1$plot_network()

# Export the network file for Gephi beautification
t1$save_network(filepath = "polar_Climate_High_function.gexf")
write.csv(canshu, "polar_Climate_High_function_canshu.csv")

# Repeat the process for another dataset
#otu = read.csv("Tropical_Climate_Low_function.csv", header = T, row.names = 1)
#tax = read.csv("Taxonomy_table_even.csv", row.names = 1)
dataset <- microtable$new(otu_table = polar_low
                          #,tax_table = tax
                          )
print(dataset)
dataset$tidy_dataset()

# Create a trans_network object and calculate the network
t1 <- trans_network$new(dataset = dataset, cor_method = "spearman", filter_thres = 0.0005)
t1$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)

# Calculate modularity
t1$cal_module()

# Calculate node table
t1$get_node_table(node_roles = )
t1$res_node_table

# Plot the ZIPI graph
t1$plot_taxa_roles(use_type = 1)

# Calculate topological parameters of the network graph
t1$cal_network_attr()
canshu <- t1$res_network_attr

# Get the adjacency matrix
t1$get_adjacency_matrix()

# Plot the network graph
t1$plot_network()

# Export the network file for Gephi beautification
t1$save_network(filepath = "otu.gexf")