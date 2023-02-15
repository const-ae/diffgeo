
<!-- README.md is generated from README.Rmd. Please edit that file -->

# diffgeo

<!-- badges: start -->
<!-- badges: end -->

The `diffgeo` R package provides functions to work with six common
manifolds:

- Euclidean arrays
- Spheres (in $n$ dimensions)
- Rotation matrices (also called special orthogonal manifolds $SO(n)$)
- Stiefel manifold (matrices with orthogonal columns)
- Grassmann manifold (matrices with orthogonal columns but they are only
  considered different if they span separate spaces)
- Symmetric positive definite (SPD) matrices (symmetric matrices with
  positive eigenvalues)

## Installation

You can install `diffgeo` from
[Github](https://github.com/const-ae/diffgeo) with:

``` r
devtools::install_github("const-ae/diffgeo")
```

## Principles

For each manifold, `diffgeo` provides:

- Random points
- Random tangents
- Exponential map (a way to go along the manifold from one point to
  another)
- Logarithm (the inverse exponential map, that is a way calculate the
  direction between points on the manifold)

## Examples

To demonstrate the functionality of the package, I will show the
functions for the rotation matrices. All functions work analoguously for
the other manifolds.

``` r
library(diffgeo)
set.seed(1)
```

To begin, I create a random rotation matrix in 2D

``` r
rotation <- rotation_random_point(2)
```

We can visualize it’s effect on a random set of points

``` r
points <- rbind(rnorm(n = 20), rnorm(20))
rotated_points <- rotation %*% points
plot(t(points), pch = 16, col = "black", asp = 1, xlim = c(-3, 3), ylim = c(-2, 2))
points(t(rotated_points), pch = 16, col = "red")
segments(points[1,], points[2,], rotated_points[1,], rotated_points[2,])
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

We can also create high-dimensional rotations and check the formal
conditions for the rotation manifold

``` r
rotation2 <- rotation_random_point(5)
# Determinant is 1
Matrix::det(rotation2)
#> [1] 1
# All columns are orthogonal
round(t(rotation2) %*% rotation2, 3)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    0    0    0    0
#> [2,]    0    1    0    0    0
#> [3,]    0    0    1    0    0
#> [4,]    0    0    0    1    0
#> [5,]    0    0    0    0    1
```

An important concept in differential geometry are geodesics, i.e. are
lines on the manifold. To find the line that connects two rotations we
use the `rotation_log` function. These directions are from the tangent
space of the manifold. The tangent space of the rotation manifold is the
set of skew-symmetric matrices.

``` r
rotation3 <- rotation_random_point(5)
direction <- rotation_log(base_point = rotation2, target_point = rotation3)
round(direction, 3)
#>        [,1]   [,2]   [,3]   [,4]   [,5]
#> [1,]  0.000 -0.483 -0.749  1.063  1.019
#> [2,]  0.483  0.000  1.771  0.734 -0.333
#> [3,]  0.749 -1.771  0.000 -1.685 -0.925
#> [4,] -1.063 -0.734  1.685  0.000  1.083
#> [5,] -1.019  0.333  0.925 -1.083  0.000
```

We can also reverse the operation and go from a base-point in the
direction of a tangent vector (“exponential map”). If we do this using
the `direction` we will land exactly at `rotation3` again.

``` r
rotation3
#>            [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,]  0.8461783 -0.2831550 -0.2681835  0.3369755  0.1353906
#> [2,] -0.1266562 -0.7590676 -0.3047431 -0.4071623 -0.3861671
#> [3,] -0.4350263 -0.5098987  0.1398897  0.6238921  0.3767557
#> [4,]  0.1644277 -0.2244352  0.3886541 -0.5630449  0.6741816
#> [5,] -0.2272798  0.1824053 -0.8152227 -0.1200224  0.4858793
rotation_map(v = direction, base_point = rotation2)
#>            [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,]  0.8461783 -0.2831550 -0.2681835  0.3369755  0.1353906
#> [2,] -0.1266562 -0.7590676 -0.3047431 -0.4071623 -0.3861671
#> [3,] -0.4350263 -0.5098987  0.1398897  0.6238921  0.3767557
#> [4,]  0.1644277 -0.2244352  0.3886541 -0.5630449  0.6741816
#> [5,] -0.2272798  0.1824053 -0.8152227 -0.1200224  0.4858793
```

Lastly, the injectivity radius is an an important concept of every
manifold. It is the largest distance for which the exponential map is a
diffeomorphism. The injectivity radius of the rotation manifold is
$\sqrt{2} \pi$ (*this appears to be only a lower bound*). That means
that for tangent vectors with a large norm, the exponential map returns
points for which there exists a more direct connection.

``` r
scaled_dir <- direction / sqrt(sum(direction^2)) * 6
sqrt(sum(scaled_dir^2))
#> [1] 6

new_point <- rotation_map(v = scaled_dir, base_point = rotation2)
rotation_log(rotation2, new_point)
#>            [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,]  0.0000000 -1.7348934  0.2927040 -0.5396579  0.4403079
#> [2,]  1.7348934  0.0000000 -0.3078207  0.8330210 -1.1194375
#> [3,] -0.2927040  0.3078207  0.0000000  2.1184384  1.4445680
#> [4,]  0.5396579 -0.8330210 -2.1184384  0.0000000  0.2491035
#> [5,] -0.4403079  1.1194375 -1.4445680 -0.2491035  0.0000000
scaled_dir
#>            [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,]  0.0000000 -0.6010657 -0.9326503  1.3232021  1.2681686
#> [2,]  0.6010657  0.0000000  2.2049703  0.9136184 -0.4148316
#> [3,]  0.9326503 -2.2049703  0.0000000 -2.0971388 -1.1516225
#> [4,] -1.3232021 -0.9136184  2.0971388  0.0000000  1.3479113
#> [5,] -1.2681686  0.4148316  1.1516225 -1.3479113  0.0000000
```

Example for the 1-D Sphere (aka a circle)

``` r
sp <- sphere_random_point(1)
tang <- sphere_random_tangent(sp)
# Make sure that norm of 'tang' is 1
tang <- tang / sqrt(sum(tang^2))


plot(t(do.call(cbind, lapply(seq(0, 2 * pi, l = 31), \(t) sphere_map(t * tang, sp)))), asp = 1)
points(t(sp), col = "red")
points(t(sphere_map(pi * tang, sp)), col = "green")
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

Example for the Grassmann manifold

``` r
poi <- grassmann_random_point(5, 2)
tan <- grassmann_random_tangent(poi)
tan <- tan / sqrt(sum(tan^2))

xg <- seq(0, 10, l = 100)
diffs <- sapply(xg, \(t){
  diff <- grassmann_log(poi, grassmann_map(t * tan, poi))
  sqrt(sum(diff^2))
})

plot(xg, diffs)
abline(v = grassmann_injectivity_radius() * (1:4), col = "red")
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

# Alternatives

This package is heavily inspired and builds on the
[Manifolds](https://juliamanifolds.github.io/Manifolds.jl/latest/index.html)
and [Manopt.jl](https://manoptjl.org/stable/) packages developed by
Ronny Bergmann et al.. The Manopt.jl package in turn is inspired by the
[Manopt](https://www.manopt.org/) toolbox for Matlab.
