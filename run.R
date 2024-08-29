# # build a matrix from your network
# network_unfiltered <- read.csv('rice_PPI', header = TRUE, sep = ' ')
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
# write.csv(adj_matrix,'rice_network_matrix.csv')

# source('neighbour_voting.R')

#print('PR_arabidopsis')
#source('page_rank_mod_arabidopsis.R')

# Save the dictionary as a JSON file
write_json(results_dict_pr_arabidopsis, "arabidopsis_page_rank_results.json")
print('PR_rice')
source('page_rank_mod_rice.R')
print('PR_maize')
source('page_rank_mod_maize.R')

print('DF_arabidopsis')
source('diffusion_mod_arabidopsis.R')
print('DF_rice')
source('diffusion_mod_rice.R')
print('DF_maize')
source('diffusion_mod_maize.R')

