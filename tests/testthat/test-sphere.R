test_that("sphere functions are consistent", {
  point <- sphere_random_point(n = 2)
  expect_true(sphere_check_point(point))

  tangent <- sphere_random_tangent(point)
  expect_equal(sum(tangent * point), 0)

  # Check that tangent scales appropriately
  hd_point <- sphere_random_point(n = 50)
  expect_lt(sum(sphere_random_tangent(hd_point, sd = 0.001)^2), 1e-3)
  expect_gt(sum(sphere_random_tangent(hd_point, sd = 10)^2), 10)

  # Check projections
  vec <- randn(3, 1)
  proj_point <- sphere_project_point(vec)
  expect_true(sphere_check_point(proj_point))
  proj_tang <- sphere_project_tangent(vec, point)
  expect_true(sphere_check_tangent(proj_tang, point))

  # Check map and log
  # Make sure that the tangent is smaller than the injectivity radius
  tangent <- tangent / (sqrt(sum(tangent^2)) * 1.2)
  expect_lt(sqrt(sum(tangent^2)), sphere_injectivity_radius())
  new_point <- sphere_map(tangent, point)
  new_tangent <- sphere_log(point, new_point)
  expect_equal(new_tangent, tangent)
})
