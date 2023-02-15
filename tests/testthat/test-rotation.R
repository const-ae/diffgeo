test_that("rotation functions are consistent", {
  point <- rotation_random_point(4)
  expect_true(rotation_check_point(point))
  expect_equal(t(point) %*% point, diag(nrow = 4))

  tangent <- rotation_random_tangent(point)
  expect_true(rotation_check_tangent(tangent, point))
  expect_equal(t(tangent) + tangent, matrix(0, nrow = 4, ncol = 4))

  # Check that tangent scales appropriately
  hd_point <- rotation_random_point(n = 50)
  expect_lt(sum(rotation_random_tangent(hd_point, sd = 0.001)^2), 0.1)
  expect_gt(sum(rotation_random_tangent(hd_point, sd = 10)^2), 10)

  # Check projections
  vec <- randn(4, 4)
  proj_point <- rotation_project_point(vec)
  expect_true(rotation_check_point(proj_point))
  proj_tang <- rotation_project_tangent(vec, point)
  expect_true(rotation_check_tangent(proj_tang, point))

  # Check map and log
  # Make sure that the tangent is smaller than the injectivity radius
  tangent <- tangent / (sqrt(sum(tangent^2)) * 2)
  expect_lt(sqrt(sum(tangent^2)), rotation_injectivity_radius())
  new_point <- rotation_map(tangent, point)
  expect_true(rotation_check_point(new_point))
  new_tangent <- rotation_log(point, new_point)
  expect_true(rotation_check_tangent(new_tangent, point))
  expect_equal(new_tangent, tangent)
})
