# Load necessary libraries
library(igraph)
library(reshape)
library(vegan)
library(bipartite)
library(car)
library(truncnorm)
library(lattice)
library(gplots)
library(ggplot2)
library(tidyr)
library(ggridges)
library(tidyverse)
library(qcmi)

# Read adjacency matrix from the network inference
otu_biotic <- read.csv("rarefied_otu_table.csv", row.names=1)

# Define a function to calculate network motifs
# This includes cycle facilitation, cycle competition, etc.
source("Code6 network_motif_functions.R")  # Assuming the functions are saved in this file

# Extract adjacency matrices from each subgraph for motif calculation
sub_netadj_biotic <- lapply(sub_graph_ig_biotic, function(g) {
  as.matrix(as_adjacency_matrix(g, type = "both", attr = 'weight'))
})

# Calculate all motifs for each subgraph
result_motifs <- lapply(sub_netadj_biotic, function(mat) {
  netmotif(mat)  # Using the netmotif function defined in 'network_motif_functions.R'
})

# Convert list of results into a data frame
result_netmotif <- do.call(rbind, result_motifs)
colnames(result_netmotif) <- c("cycfac", "cyccom", "facmcom", "commfac", "tranfac", "trancom", "trancomfac")

# Write the motif analysis results to file
write.csv(result_netmotif, "result_netmotif.csv")

# Export graphs of subnetworks for each subgraph
lapply(seq_len(length(sub_graph_ig_biotic)), function(x) {
  write_graph(sub_graph_ig_biotic[[x]], paste("edge", x, ".txt", sep = ""), format = "ncol")
})

# Create random networks and calculate motifs for them
random_motif_deposit <- lapply(seq_len(length(sub_netadj_biotic)), function(i) {
    message(paste("Processing network", i))
    edge_data <- read.table(paste("edge", i, ".txt", sep = ""), header = F)
    random_motifs <- replicate(99, {
        shuffled_edges <- edge_data[sample(nrow(edge_data)), ]
        ig_shuffled <- graph_from_edgelist(as.matrix(shuffled_edges[, 1:2]), directed = FALSE)
        ig_shuffled <- set_edge_attr(ig_shuffled, 'weight', index = E(ig_shuffled), as.numeric(edge_data[, 3]))
        mat <- as.matrix(as_adjacency_matrix(ig_shuffled, type = "both", attr = "weight"))
        netmotif(mat)
    }, simplify = FALSE)
    do.call(rbind, random_motifs)
})

# Calculate Z-scores for motifs against random distributions
result_z_scores <- lapply(1:nrow(result_netmotif), function(i) {
  net <- result_netmotif[i, ]  # 提取第 i 行的模体向量
  rand <- random_motif_deposit[[i]]  # 提取第 i 个随机网络矩阵
  mean_rand <- colMeans(rand)
  sd_rand <- apply(rand, 2, sd)
  z_scores <- (net - mean_rand) / sd_rand
  return(z_scores)
})



# Calculate P-values for observed motifs vs. random expectations
result_p_values <- t(sapply(1:nrow(result_netmotif), function(i) {
  net <- result_netmotif[i, ]
  rand <- random_motif_deposit[[i]]
  
  p_values <- sapply(1:7, function(j) {
    o <- net[j]
    r <- rand[, j]
    t.test(r, mu = o, alternative = "two.sided")$p.value
  })
  
  return(p_values)
}))


# Format results and save
df_z_scores <- as.data.frame(do.call(rbind, result_z_scores))
colnames(df_z_scores) <- colnames(result_netmotif)
row.names(df_z_scores) = colnames(otu_biotic) 
write.csv(df_z_scores, "motif_z_scores.csv")

colnames(result_p_values) <- colnames(result_netmotif)
row.names(result_p_values) = colnames(otu_biotic) 
write.csv(result_p_values, "motif_p_values.csv")



