test_that("Length of 8 of the resulting list", {
  expect_equal(length(map_clusters_to_representation(M = matrix(seq(1,100), ncol=10), group = c(rep(c('A','B','C'), each = 3), 'D'))), 5)
})

test_that("Type of output", {
  expect_s4_class(map_clusters_to_representation(M = matrix(seq(1,100), ncol=10), group = c(rep(c('A','B','C'), each = 3), 'D'))
                  , 'SimpleList')
})
