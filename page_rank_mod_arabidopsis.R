#igraph
library(igraph)

# reading in the network and filtering on edge weight
network_unfiltered_arabidopsis <- read.csv('arabidopsis_PPI', header = TRUE, sep = ' ')
network_arabidopsis <- network_unfiltered_arabidopsis[network_unfiltered_arabidopsis$combined_score >= 700, ]

# creating the graph
g_arabidopsis <- graph_from_data_frame(network_arabidopsis, directed = FALSE)

# set scores for edges
E(g_arabidopsis)$weight <- network_arabidopsis$combined_score

# get the functional labels
library(jsonlite)
seed_proteins_unfiltered_arabidopsis <- fromJSON('arabidopsis_seven.json')

# remove GO terms with seed proteins less than 10
seed_proteins_arabidopsis <- seed_proteins_unfiltered_arabidopsis[sapply(seed_proteins_unfiltered_arabidopsis, length)>=10]

# get the go terms
go_terms_arabidopsis <- names(seed_proteins_arabidopsis)

# extract the proteins and shuffle
all_proteins_arabidopsis <- unlist(seed_proteins_arabidopsis)
shuffled_proteins_arabidopsis <- sample(all_proteins_arabidopsis)
split_proteins_arabidopsis <- split(shuffled_proteins_arabidopsis, cut(seq_along(shuffled_proteins_arabidopsis), breaks=10, labels= FALSE))

# Create an empty dictionary to store the results
results_dict_pr_arabidopsis <- list()

start_time <- proc.time()

for (go_term in go_terms_arabidopsis){
  print(go_term)
  # Open a new key for the GO term
  results_dict_pr_arabidopsis[[go_term]] <- list()
  
  # # get all the nodes
  # all_nodes <- V(g)$name
  # # personalization vector is created
  # personalization <- rep(0, length(all_nodes))
  # names(personalization) <- all_nodes
  
  # Iterate through each portion (1 to 10)
  for (portion in seq_along(split_proteins_arabidopsis)) {
    # Get the seed proteins from the remaining portions (excluding the selected portion)
    remaining_proteins <- unlist(split_proteins_arabidopsis[-portion])
    seed_proteins_unselected <- intersect(seed_proteins_arabidopsis[[go_term]], remaining_proteins)
    
    # Assign seeds as initial probabilities
    initial_probs <- ifelse(V(g_arabidopsis)$name %in% seed_proteins_unselected, 1, 0)

    # calculate page rank
    ppr <- page_rank(g_arabidopsis, personalized = initial_probs)$vector
    
    # Store the output probabilities for the selected portion
    results_dict_pr_arabidopsis[[go_term]][split_proteins_arabidopsis[[portion]]] <- ppr #ppr[split_proteins[[portion]]]
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

end_time <- proc.time()
elapsed_time = end_time - start_time
elapsed_seconds = elapsed_time["elapsed"]
cat(sprintf(paste0('Excecution time for page rank in arabidopsis in seconds\n', elapsed_seconds, '\n')),file = 'elapsed_time.txt', append = TRUE)

# Save the dictionary as a JSON file
write_json(results_dict_pr_arabidopsis, "arabidopsis_page_rank_results.json")
