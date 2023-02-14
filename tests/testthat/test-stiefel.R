test_that("Stiefel functions are consistent", {
  point <- stiefel_random_point(4, 2)
  expect_true(stiefel_check_point(point))
  expect_equal(t(point) %*% point, diag(nrow = 2))

  tangent <- stiefel_random_tangent(point)
  expect_true(stiefel_check_tangent(tangent, point))
  expect_equal(t(point) %*% tangent + t(tangent) %*% point, matrix(0, nrow = 2, ncol = 2))

  # Check that tangent scales appropriately
  hd_point <- stiefel_random_point(n = 50, k = 20)
  expect_lt(sum(stiefel_random_tangent(hd_point, sd = 0.001)^2), 0.001)
  expect_gt(sum(stiefel_random_tangent(hd_point, sd = 10)^2), 10)

  # Check projections
  vec <- randn(4, 2)
  proj_point <- stiefel_project_point(vec)
  expect_true(stiefel_check_point(proj_point))
  proj_tang <- stiefel_project_tangent(vec, point)
  expect_true(stiefel_check_tangent(proj_tang, point))

  # Check map and log
  # Make sure that the tangent is smaller than the injectivity radius
  tangent <- tangent / (sqrt(sum(tangent^2)) * 5)
  expect_lt(sqrt(sum(tangent^2)), stiefel_injectivity_radius())
  new_point <- stiefel_map(tangent, point)
  expect_true(stiefel_check_point(new_point))
  new_tangent <- stiefel_log(point, new_point)
  expect_true(stiefel_check_tangent(new_tangent, point))
  expect_equal(new_tangent, tangent)
})
