# Load required libraries
library(phyloseq)
library(SpiecEasi)
library(devtools)
library(igraph)
library(Matrix)

# Read and preprocess OTU table
data <- read.delim("table.txt", row.names=1, header=T, sep='\t')
OTU <- otu_table(data, taxa_are_rows = TRUE)

# Filter taxa that appear more than once across samples
OTU1 <- filter_taxa(OTU, function(x) sum(x > 0) > 1, TRUE)

# Filter taxa with total counts less than 10
OTU2 <- prune_taxa(taxa_sums(OTU1) > 9, OTU1)

# Rarefaction to even depth across samples
G <- rarefy_even_depth(OTU2, sample.size = min(sample_sums(OTU2)),
                       rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
write.csv(G, "rarefied_otu_table.csv")

# Read the rarefied OTU table for network construction
otu <- read.csv("rarefied_otu_table.csv", row.names=1)

# Construct microbial network using SpiecEasi
se.mb.amgut <- spiec.easi(t(otu), method='mb', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50))
ig.mb <- adj2igraph(getRefit(se.mb.amgut))

# Extract and summarize the network edges
sebeta <- symBeta(getOptBeta(se.mb.amgut), mode='maxabs')
sebeta1 <- as.matrix(sebeta)
elist.mb <- summary(sebeta)

# Construct an undirected graph from edge list
ig <- graph_from_edgelist(as.matrix(elist.mb[,1:2]), directed = FALSE)
ig <- set_edge_attr(ig, 'weight', index = E(ig), as.numeric(elist.mb[,3]))

# Filter nodes based on OTU counts and create subgraphs
otu_filter <- otu
ig_biotic <- ig
otutab <- otu_filter
g_biotic <- ig_biotic
node <- vertex_attr(g_biotic)[[1]]
name <- intersect(row.names(otutab), node)
otu <- otutab[name,]

# Display the number of nodes and edges in the biotic graph
print(length(V(g_biotic))) # Number of nodes
print(length(E(g_biotic))) # Number of edges

# Generate subgraphs for samples with node values greater than 10
sub_graph_ig_biotic3 <- list()
for (i in names(otu)) {
    sample_i <- otu[i]
    select_node <- rownames(sample_i)[which(sample_i > 10)]
    sub_graph_ig_biotic3[[i]] <- induced_subgraph(g_biotic, select_node)
    print(i)
}
