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

# extract the proteins and shuffle
all_proteins <- unlist(seed_proteins)
shuffled_proteins <- sample(all_proteins)
split_proteins <- split(shuffled_proteins, cut(seq_along(shuffled_proteins), breaks=10, labels= FALSE))

# Create an empty dictionary to store the results
results_dict <- list()


for (go_term in go_terms){
  # Open a new key for the GO term
  results_dict[[go_term]] <- list()
  
  # # get all the nodes
  # all_nodes <- V(g)$name
  # # personalization vector is created
  # personalization <- rep(0, length(all_nodes))
  # names(personalization) <- all_nodes
  
  # Iterate through each portion (1 to 10)
  for (portion in seq_along(split_proteins)) {
    # Get the seed proteins from the remaining portions (excluding the selected portion)
    remaining_proteins <- unlist(split_proteins[-portion])
    seed_proteins_unselected <- intersect(go_dict[[go_term]], remaining_proteins)
    
    # Assign seeds as initial probabilities
    initial_probs <- ifelse(V(g)$name %in% seed_proteins_unselected, 1, 0)

    # calculate page rank
    ppr <- page_rank(g, personalized = initial_probs)$vector
    
    # Store the output probabilities for the selected portion
    results_dict[[go_term]][split_proteins[[portion]]] <- ppr[split_proteins[[portion]]]
    # Print the Personalized Page Rank scores
    # print(ppr$vector)
    # extract the values of the selected proteins
  }
                          
  # personalization[associated_proteins] <- 2/ length(associated_proteins)
  # # calculate pagerank
  # ppr <- page_rank(g, personalized = personalization)
  # # Print the Personalized PageRank scores
  # print(ppr$vector)
  # Plot the graph with Personalized PageRank scores as labels
  # plot(g, vertex.label = round(ppr$vector, 3), vertex.size = 20,
  #      vertex.label.cex = 1.2, main = "Personalized PageRank (PPR)")
}

# Save the dictionary as a JSON file
write_json(results_dict, "output_results.json")
