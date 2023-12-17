library(inferCSN)

data("example_matrix")
weight_table <- inferCSN(example_matrix, verbose = TRUE)
