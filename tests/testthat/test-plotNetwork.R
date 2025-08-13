m <- matrix(sample(seq(1,10, length.out=10000), 15000*100, replace = TRUE), ncol = 100)
rownames(m) <- paste0('Gene-', seq(1,15000))
colnames(m) <- paste0('Col-', seq(1,100))
net <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = m))
net$SubClass <- rep(c('A', 'B', 'C', 'D'), each = 25)
subclass_palette <- c('A' = 'red', 'B' = 'blue', 'C' = 'green', 'D' = 'yellow')
net$SubClass_colors <- subclass_palette[net$SubClass]
net$X_coord <- sample(seq(1,2, length.out = 1000), size = ncol(net), replace = TRUE)
net$Y_coord <- sample(seq(1,2, length.out = 1000), size = ncol(net), replace = TRUE)
edges_from <- sample(colnames(net), size = 200, replace = TRUE)
edges_to <- sample(colnames(net), size = 200, replace = TRUE)
edges_from_x <- net$X_coord[match(edges_from, colnames(net))]
edges_from_y <- net$Y_coord[match(edges_from, colnames(net))]
edges_to_x <- net$X_coord[match(edges_to, colnames(net))]
edges_to_y <- net$Y_coord[match(edges_to, colnames(net))]
edges_weight <- sample(seq(0,1, length.out=1000), length(edges_from), replace = TRUE)
edges_df <- data.frame('from' = edges_from, 'to' = edges_to, 'weight' = edges_weight, 'from.x' = edges_from_x, 'from.y' = edges_from_y, 'to.x' = edges_to_x, 'to.y' = edges_to_y)
net@metadata$network$edges <- edges_df

test_that("Type output", {
  expect_s3_class(plotNetwork(net), 'ggplot')
})
