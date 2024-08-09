# #using EGAD on the network
library(EGAD)
# 
# # build a matrix from your network
# network_unfiltered <- read.csv('arabidopsis_PPI', header = TRUE, sep = ' ')
# network <- network_unfiltered[network_unfiltered$combined_score >= 700, ]
# 
# # Extract unique protein names
# proteins <- unique(c(network$protein1, network$protein2))
# 
# # Create an empty matrix
# adj_matrix <- matrix(0, nrow = length(proteins), ncol = length(proteins),
#                      dimnames = list(proteins, proteins))
# 
# # Fill the matrix with edge weights
# for (i in 1:nrow(network)) {
#   p1 <- network[i, 1]
#   p2 <- network[i, 2]
#   weight <- network[i, 3]
#   
#   adj_matrix[p1, p2] <- weight
#   adj_matrix[p2, p1] <- weight # Fill the inverse as well
# }
# 
# write.csv(adj_matrix,'arabidopsis_network_matrix.csv')
# 
# # get the functional labels
# library(jsonlite)
# seed_proteins_unfiltered <- fromJSON('arabidopsis_seven.json')
# 
# # remove GO terms with seed proteins less than 10
# seed_proteins <- seed_proteins_unfiltered[sapply(seed_proteins_unfiltered, length)>=10]
# 
# # get the go terms
# go_terms <- names(seed_proteins)
# 
# # get the unique proteins
# proteins_go <- unique(unlist(seed_proteins))
# 
# # create the seed matrix 
# go_matrix <- matrix(0, nrow=length(proteins_go), ncol = length(go_terms), dimnames = list(proteins_go, go_terms)) 
# 
# # populate the matrix
# for (go_term in go_terms){
#   associated_proteins <- seed_proteins[[go_term]]
#   go_matrix[associated_proteins, go_term] <- 1
# }
# 
# write.csv(go_matrix,'arabidopsis_functions_matrix.csv')

aurocs <- neighbor_voting(go_matrix, adj_matrix, output = 'AUROC')

avgprcs <- neighbor_voting(go_matrix, adj_matrix, output = 'PR')