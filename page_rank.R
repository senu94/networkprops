#igraph
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

for (go_term in go_terms){
  # get all the nodes
  all_nodes <- V(g)$name
  # personalization vector is created
  personalization <- rep(0, length(all_nodes))
  names(personalization) <- all_nodes
  # assign a higher value to the seed proteins
  associated_proteins <- c(seed_proteins[[go_term]])
  personalization[associated_proteins] <- 2/ length(associated_proteins)
  # calculate pagerank
  ppr <- page_rank(g, personalized = personalization)
  # Print the Personalized PageRank scores
  print(ppr$vector)
  # Plot the graph with Personalized PageRank scores as labels
  # plot(g, vertex.label = round(ppr$vector, 3), vertex.size = 20,
  #      vertex.label.cex = 1.2, main = "Personalized PageRank (PPR)")
}
