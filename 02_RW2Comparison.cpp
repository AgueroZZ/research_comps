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
  DATA_SCALAR(u1); // pc prior, u1 param
  DATA_SCALAR(alpha1); // pc prior, alpha1 param
  DATA_SCALAR(u2); // pc prior, u2 param
  DATA_SCALAR(alpha2); // pc prior, alpha2 param

  // Parameter
  PARAMETER_VECTOR(W); // W = c(U,beta), eta = BX * U + X * beta
  PARAMETER(theta1); // theta = -2log(sigma1)
  PARAMETER(theta2); // theta = -2log(sigma2)

  // Transformations
  vector<Type> eta = X * W;
  Type sigma1 = exp(-0.5*theta1);
  REPORT(sigma1);
  Type sigma2 = exp(-0.5*theta2);
  REPORT(sigma2);

  // Log likelihood
  Type ll = 0;
  ll = sum(dnorm(y, eta, sigma2, TRUE));
  REPORT(ll);
  
  // Log prior on W
  Type lpW = 0;

  // Cross product
  vector<Type> PW = P*W;
  Type WPW = (W * PW).sum();
  lpW += -0.5 * exp(theta1) * WPW; // U=W part

  // Log determinant
  Type logdet1 = d * theta1 + logPdet;
  lpW += 0.5 * logdet1; // P part

  REPORT(logdet1);
  REPORT(lpW);

  
  // Log prior for theta
  Type lpT = 0;
  Type phi1 = -log(alpha1) / u1;
  lpT += log(0.5 * phi1) - phi1*exp(-0.5*theta1) - 0.5*theta1;
  Type phi2 = -log(alpha2) / u2;
  lpT += log(0.5 * phi2) - phi2*exp(-0.5*theta2) - 0.5*theta2;
  REPORT(lpT);
  
  // Final result!
  Type logpost = -1 * (ll + lpW + lpT);
  
  return logpost;
}