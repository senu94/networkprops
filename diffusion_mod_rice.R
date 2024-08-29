# using diffuStats
library(diffuStats)
library(igraph)

# get the network
# network_unfiltered <- read.csv('rice_PPI', header = TRUE, sep = ' ')
# network <- network_unfiltered[network_unfiltered$combined_score >= 700, ]

# creating the graph
# g <- graph_from_data_frame(network, directed = FALSE)

# set scores for edges
# E(g)$weight <- network$combined_score

# get the functional labels
# library(jsonlite)
# seed_proteins_unfiltered <- fromJSON('rice_seven.json')

# remove GO terms with seed proteins less than 10
# seed_proteins <- seed_proteins_unfiltered[sapply(seed_proteins_unfiltered, length)>=10]

# get the go terms
# go_terms <- names(seed_proteins)

# get the nodes
proteins_go_rice <- V(g_rice)$name

# extract the proteins and shuffle
# all_proteins <- unlist(seed_proteins)
# shuffled_proteins <- sample(all_proteins)
# split_proteins <- split(shuffled_proteins, cut(seq_along(shuffled_proteins), breaks=10, labels= FALSE))

# # create the empty vector
# go_vector <- setNames(numeric(length(proteins_go)), proteins_go)
# print(go_vector)

# Create an empty dictionary to store the results
results_dict_diffml_rice <- list()

start_time <- proc.time()

for (go_term in go_terms_rice){
  print(go_term)
  # Open a new key for the GO term
  results_dict_diffml_rice[[go_term]] <- list()
  # Iterate through each portion (1 to 10)
  for (portion in seq_along(split_proteins_rice)) {
    # Get the seed proteins from the remaining portions (excluding the selected portion)
    remaining_proteins <- unlist(split_proteins_rice[-portion])
    seed_proteins_unselected <- intersect(seed_proteins_rice[[go_term]], remaining_proteins)
  
    go_vector <- setNames(ifelse(proteins_go_rice %in% seed_proteins_unselected, 1, 0), proteins_go_rice)
    output_vec <- diffuStats::diffuse(
      graph = g_rice,
      method = "ml",
      scores = go_vector)
    
    # Store the output probabilities for the selected portion
    results_dict_diffml_rice[[go_term]][split_proteins[[portion]]] <- output_vec
  }
}

end_time <- proc.time()
elapsed_time = end_time - start_time
elapsed_seconds = elapsed_time["elapsed"]
cat(sprintf(paste0('Excecution time for diffusion ml in rice in seconds\n', elapsed_seconds)),file = 'elapsed_time.txt', append = TRUE)

# Save the dictionary as a JSON file
write_json(results_dict_diffml_rice, "rice_diffml_results.json")
