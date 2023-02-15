test_that("euclidean functions are consistent", {
  point <- euclidean_random_point(c(4, 2))
  expect_true(euclidean_check_point(point))

  tangent <- euclidean_random_tangent(point)
  expect_true(euclidean_check_tangent(tangent, point))

  # Check that tangent scales appropriately
  hd_point <- euclidean_random_point(c(50, 20))
  expect_lt(sum(euclidean_random_tangent(hd_point, sd = 0.001)^2), 0.1)
  expect_gt(sum(euclidean_random_tangent(hd_point, sd = 10)^2), 10)

  # Check projections
  vec <- randn(4, 2)
  proj_point <- euclidean_project_point(vec)
  expect_true(euclidean_check_point(proj_point))
  proj_tang <- euclidean_project_tangent(vec, point)
  expect_true(euclidean_check_tangent(proj_tang, point))

  # Check map and log
  expect_lt(sqrt(sum(tangent^2)), euclidean_injectivity_radius())
  new_point <- euclidean_map(tangent, point)
  new_tangent <- euclidean_log(point, new_point)
  expect_equal(new_tangent, tangent)
})
