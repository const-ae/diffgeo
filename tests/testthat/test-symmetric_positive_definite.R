test_that("spd functions are consistent", {
  point <- spd_random_point(4)
  expect_true(spd_check_point(point))

  tangent <- spd_random_tangent(point)
  expect_true(spd_check_tangent(tangent, point))
  expect_equal(t(tangent) - tangent, matrix(0, nrow = 4, ncol = 4))

  # Check that tangent scales appropriately
  hd_point <- spd_random_point(n = 50)
  expect_lt(sum(spd_random_tangent(hd_point, sd = 0.001)^2), 0.1)
  expect_gt(sum(spd_random_tangent(hd_point, sd = 10)^2), 10)

  # Check projections
  vec <- randn(4, 4)
  proj_point <- spd_project_point(vec)
  expect_true(spd_check_point(proj_point))
  proj_tang <- spd_project_tangent(vec, point)
  expect_true(spd_check_tangent(proj_tang, point))

  # Check map and log
  # Make sure that the tangent is smaller than the injectivity radius
  expect_lt(sqrt(sum(tangent^2)), spd_injectivity_radius()) # is anyways inf
  new_point <- spd_map(tangent, point)
  expect_true(spd_check_point(new_point))
  new_tangent <- spd_log(point, new_point)
  expect_true(spd_check_tangent(new_tangent, point))
  expect_equal(new_tangent, tangent)
})

