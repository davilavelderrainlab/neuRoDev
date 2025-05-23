S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
annotated_M <- reference_signatures_correlation(S, refS)
new_clusterS <- FC_signatures(matrix(runif(80,0,10), ncol = 4))
rownames(new_clusterS) <- paste0('Gene-', seq(1, dim(new_clusterS)[1]))
colnames(new_clusterS) <- paste0('New-', seq(1, dim(new_clusterS)[2]))
new_M <- reference_signatures_correlation(new_clusterS, refS)
res <- add_to_reference(annotated_M, new_M, annotated_M$`Best.Assignment`)
umap_obj <- res$New
new_clusters = rownames(new_M)

test_that("Output length", {
  expect_equal(length(nearest_neighbor_annotation(annotated_M, new_clusters, umap_obj, color_attr = 'Best.Assignment')), 7)
})

test_that("Returns a list", {
  expect_s4_class(nearest_neighbor_annotation(annotated_M, new_clusters, umap_obj, color_attr = 'Best.Assignment')
                  , 'SimpleList')
})
