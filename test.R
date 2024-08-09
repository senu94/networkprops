# Install Bioconductor if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#install EGAD
# Use BiocManager to install EGAD
BiocManager::install("EGAD")

# Use BiocManager to install diffustats
BiocManager::install("diffuStats")

#EGAD
library(EGAD)
genes.labels <- matrix( sample( c(0,1), 1000, replace=TRUE), nrow=100)
rownames(genes.labels) = paste('gene', 1:100, sep='')
colnames(genes.labels) = paste('function', 1:10, sep='')
net <- cor( matrix( rnorm(10000), ncol=100), method='spearman')
rownames(net) <- paste('gene', 1:100, sep='') 
colnames(net) <- paste('gene', 1:100, sep='')

aurocs <- neighbor_voting(genes.labels, net, output = 'AUROC')

avgprcs <- neighbor_voting(genes.labels, net, output = 'PR')

#igraph
library(igraph)

# Parameters for the random network
num_nodes <- 20  # Number of nodes (e.g., genes or proteins)
probability <- 0.1  # Probability of edge creation

# Generate the random graph using the Erdős-Rényi model
g <- erdos.renyi.game(num_nodes, probability, directed = FALSE)

# Create a personalization vector (focus on the first node)
personalization_vector <- rep(0, vcount(g))
personalization_vector[1] <- 1

# Compute Personalized PageRank
ppr <- page_rank(g, personalized = personalization_vector)

# Compute PageRank with a uniform prior
pr = page_rank(g)

# Print the Personalized PageRank scores
print(ppr$vector)

# Plot the graph with Personalized PageRank scores as labels
plot(g, vertex.label = round(ppr$vector, 3), vertex.size = 20,
     vertex.label.cex = 1.2, main = "Personalized PageRank (PPR)")

# using diffuStats
library(diffuStats)

# use an example graph
data("graph_toy")

input_vec <- graph_toy$input_vec

# using raw diffusion
output_vec <- diffuStats::diffuse(
  graph = graph_toy, 
  method = "raw", 
  scores = input_vec)

# using GeneMania-based weights
output_vec <- diffuStats::diffuse(
  graph = graph_toy, 
  method = "gm", 
  scores = input_vec)

# using Monte Carlo normalized scores
output_vec <- diffuStats::diffuse(
  graph = graph_toy, 
  method = "mc", 
  scores = input_vec)

# using z scores
output_vec <- diffuStats::diffuse(
  graph = graph_toy, 
  method = "z", 
  scores = input_vec)

# using ml
output_vec <- diffuStats::diffuse(
  graph = graph_toy, 
  method = "ml", 
  scores = input_vec)

plot(graph_toy)

igraph::plot.igraph(
  graph_toy, 
  vertex.color = diffuStats::scores2colours(output_vec),
  vertex.shape = diffuStats::scores2shapes(input_vec),
  main = "Diffusion scores in the lattice"
)

# using random forests
library(randomForest)

# Create a sample network
edges <- c(1, 2, 2, 3, 3, 4, 4, 1, 2, 4)
g <- graph(edges, directed = FALSE)

# Assign names to the nodes
V(g)$name <- as.character(1:vcount(g))

# Calculate some network-based features for nodes
degree <- degree(g)
betweenness <- betweenness(g)
closeness <- closeness(g)

# Combine features into a data frame
features_df <- data.frame(degree, betweenness, closeness)
rownames(features_df) <- V(g)$name

# Example node labels (for training data; some nodes may have unknown labels)
node_labels <- c("FunctionA", "FunctionB", "FunctionA", NA)
names(node_labels) <- V(g)$name

# Prepare data for training (excluding nodes with missing labels)
train_indices <- !is.na(node_labels)
train_data <- features_df[train_indices, ]
train_labels <- node_labels[train_indices]

# Nodes for which to predict labels
predict_indices <- is.na(node_labels)
predict_data <- features_df[predict_indices, ]

# Train the Random Forest model
rf_model <- randomForest(train_data, as.factor(train_labels), ntree = 100)

# Predict labels for nodes with unknown labels
predicted_labels <- predict(rf_model, predict_data)

# Combine known and predicted labels
node_labels[predict_indices] <- predicted_labels

# Print the results
print(node_labels)

# using an SVM kernel
# Train the SVM model with RBF kernel
svm_model <- ksvm(as.matrix(train_data), as.factor(train_labels), 
                  kernel = "vanilladot", 
                  C = 1)

# Predict labels for nodes with unknown labels
predicted_labels <- predict(svm_model, as.matrix(predict_data))

# Combine known and predicted labels
node_labels[predict_indices] <- predicted_labels

# Print the results
print(node_labels)

# Core
library(igraph)
library(diffuStats)
# Plotting
library(ggplot2)
library(ggsci)
# Data
library(igraphdata)
data(yeast)
set.seed(1)
