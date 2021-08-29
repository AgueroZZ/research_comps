#### Try the spline constructing function:
## what values do you need? use X to denote (we have 101 locations of interest)
X <- seq(0,1, by = 0.01)
## how many knots do you have? use knots to denote knots (we have 11 knots)
knots <- seq(0,1, by = 0.1)
## Define the b-spline basis of order 1: use splineBasis to denote
splineBasis <- splines::bs(X, knots = knots, degree = 1, intercept = F)
## splineBasis is a 101x12 matrix, each row denote the values of each function in each location
## the last column is all zero

### Ploting:
plot(splineBasis[,1]~X, type = "l", col = "blue")
for (i in 2:11) {
  lines(splineBasis[,i]~X, col = "blue")
}


### now suppose we already know the weights vector:
weights <- rnorm(n = length(knots))
fitted_func <- splineBasis[,-12] %*% matrix(weights, ncol = 1)
plot(fitted_func~X, type = "l", col = "purple")


















#### For the ARIMA Approach: There are two interpolating methods based on different 
#### interpretation of the smoothing problem: starting with cubic spline or starting
#### with linear spline approximation for the SDE

### Interpretation 1: Starting with Linear spline (method same as before)
construct_fitted <- function(weights){
  knots <- x
  splineFunc <- splines::bs(location_interest, degree = 1, knots = knots)
  fitted_func <- splineFunc[,-(length(knots) + 1)] %*% matrix(weights, ncol = 1)
  as.numeric(fitted_func)
}


### Interpretation 2: Starting with Cubic spline assumption
construct_fitted_cubic <- function(weights){
  x <- location_interest
  spline_fitted <- spline(x, y = weights, xout = x, method = "natural",
                          xmin = min(x), xmax = max(x), ties = mean)
  spline_fitted$x
}







