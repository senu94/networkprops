# using diffuStats
library(diffuStats)
library(igraph)

network_unfiltered <- read.csv('arabidopsis_PPI', header = TRUE, sep = ' ')
network <- network_unfiltered[network_unfiltered$combined_score >= 700, ]

# creating the graph
g <- graph_from_data_frame(network, directed = FALSE)

# set scores for edges
E(g)$weight <- network$combined_score

# get the functional labels
library(jsonlite)
seed_proteins_unfiltered <- fromJSON('arabidopsis_seven.json')

# remove GO terms with seed proteins less than 10
seed_proteins <- seed_proteins_unfiltered[sapply(seed_proteins_unfiltered, length)>=10]

# get the go terms
go_terms <- names(seed_proteins)

# get the nodes
proteins_go <- V(g)$name

# # create the empty vector
# go_vector <- setNames(numeric(length(proteins_go)), proteins_go)
# print(go_vector)


for (go_term in go_terms){
  associated_proteins <- c(seed_proteins[[go_term]])
  go_vector <- setNames(ifelse(proteins_go %in% associated_proteins, 1, 0), proteins_go)
  output_vec <- diffuStats::diffuse(
    graph = g,
    method = "raw",
    scores = go_vector)
  print(output_vec)
}
