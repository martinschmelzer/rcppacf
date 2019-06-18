#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix recAcf(int lmax, List ProdPairs, NumericVector yCummean, NumericVector yCumsum) {
  int n = yCummean.size();
  double yCumsumH = 0.0;
  NumericMatrix racfs(lmax+1,n);
  for(int h = 0; h < lmax+1; ++h) {
    NumericVector vh = ProdPairs[h];
    if (h == 0) {
      yCumsumH = 0;
    } else {
      yCumsumH = yCumsum[h-1];
    }
    for(int j = h; j < n; ++j) {
      racfs(h,j) = vh[j-h] - yCummean[j]/(j+1)*(yCumsum[j] - yCumsumH + yCumsum[j-h] - (j+1-h)*yCummean[j]);
    }
  }
  return racfs;
}