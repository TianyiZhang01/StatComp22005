#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
using namespace Rcpp;



//' @title a Gibbs sampler using Rcpp
//' @description a Gibbs sampler using Rcpp, which comes from my homework 10.
//' @param inix the initial value of the markov chain
//' @param iniy the initial value of the markov chain
//' @param mu1 the mean of X of the target distribution
//' @param mu2 the mean of Y of the target distribution
//' @param sigma1 the standard deviation of X of the target distribution
//' @param sigma2 the standard deviation of Y of the target distribution
//' @param rho the correlation between X and Y of the target distribution
//' @param n the quantity of samples would like to generate, with initial 10000
//' @param burn the quantity of samples dropped before the chain go to steady
//' @return a random sample of size \code{n+burn}
//' @examples
//' \dontrun{
//' temp = gibbsC(100, 10, 0, 0, 1, 1, 0.9)[-(1:1000), ]
//' par(mfrow=c(2,1));
//' plot(temp[-(1:1000),1],type='l')
//' plot(temp[-(1:1000),2],type='l')
//' }
//' @export
// [[Rcpp::export]]

NumericMatrix gibbsC(double inix, double iniy, double mu1, double mu2, double sigma1, double sigma2, double rho, int n = 10000, int burn = 1000){
  NumericMatrix a(n+burn, 2);
  double s1 = sqrt(1 - rho * rho) * sigma1, s2 = sqrt(1 - rho * rho) * sigma2 ;

  a(0,0) = inix, a(0,1) = iniy;

  for (int i = 1; i <= n+burn-1; i++)
  {
    double tempy = a(i-1,1);
    double m1 = mu1 + rho * (tempy - mu2) * sigma1 / sigma2;
    a(i,0) = rnorm(1, m1, s1)[0];

    double tempx = a(i,0);
    double m2 = mu2 + rho * (tempx - mu1) * sigma2/sigma1;
    a(i,1) = rnorm(1, m2, s2)[0];
  }

  return (a);
}
