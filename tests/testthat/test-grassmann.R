test_that("Grassmann functions are consistent", {
  point <- grassmann_random_point(4, 2)
  expect_true(grassmann_check_point(point))
  expect_equal(t(point) %*% point, diag(nrow = 2))

  tangent <- grassmann_random_tangent(point)
  expect_true(grassmann_check_tangent(tangent, point))
  expect_equal(t(point) %*% tangent + t(tangent) %*% point, matrix(0, nrow = 2, ncol = 2))

  # Check that tangent scales appropriately
  hd_point <- grassmann_random_point(n = 50, k = 20)
  expect_lt(sum(grassmann_random_tangent(hd_point, sd = 0.001)^2), 0.001)
  expect_gt(sum(grassmann_random_tangent(hd_point, sd = 10)^2), 10)

  # Check projections
  vec <- randn(4, 2)
  proj_point <- grassmann_project_point(vec)
  expect_true(grassmann_check_point(proj_point))
  proj_tang <- grassmann_project_tangent(vec, point)
  expect_true(grassmann_check_tangent(proj_tang, point))

  # Check map and log
  # Make sure that the tangent is smaller than the injectivity radius
  tangent <- tangent / (sqrt(sum(tangent^2)) * 1.2)
  expect_lt(sqrt(sum(tangent^2)), grassmann_injectivity_radius())
  new_point <- grassmann_map(tangent, point)
  new_tangent <- grassmann_log(point, new_point)
  expect_equal(new_tangent, tangent)
})
