#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Read in data
  DATA_VECTOR(y); //response variable
  DATA_SPARSE_MATRIX(X); // Design matrix
  DATA_SPARSE_MATRIX(P); // Penalty matrix
  
  int d = P.cols(); // Number of B-Spline coefficients
  DATA_SCALAR(logPdet); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(u); // pc prior, u param
  DATA_SCALAR(alpha); // pc prior, alpha param

  // Parameter
  PARAMETER_VECTOR(W); // W = c(U,beta), eta = BX * U + X * beta
  PARAMETER(theta); // theta = -2log(sigma)
  
  
  // Transformations
  vector<Type> eta = X * W;
  Type sigma = exp(-0.5*theta);
  REPORT(sigma);
  
  // Log likelihood
  Type ll = 0;
  ll = sum(dnorm(y, eta, 1, TRUE));
  REPORT(ll);
  
  // Log prior on W
  Type lpW = 0;
  // Cross product
  vector<Type> PW = P*W;
  Type WPW = (W * PW).sum();
  lpW += -0.5 * exp(theta) * WPW; // U=W part

  // Log determinant
  Type logdet1 = d * theta + logPdet;
  lpW += 0.5 * logdet1; // P part

  REPORT(logdet1);
  REPORT(lpW);

  
  // Log prior for theta
  Type lpT = 0;
  Type phi = -log(alpha) / u;
  lpT += log(0.5 * phi) - phi*exp(-0.5*theta) - 0.5*theta;
  REPORT(lpT);
  
  // Final result!
  Type logpost = -1 * (ll + lpW + lpT);
  
  return logpost;
}