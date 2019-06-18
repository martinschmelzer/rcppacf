#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rollAcf(int lmax, double m, NumericMatrix ProdPairs, NumericVector yRollmean, List yRollsum) {
  int nm = yRollmean.size();
  NumericMatrix rollacfs(nm,lmax+1);
  for(int h = 0; h < lmax+1; ++h) {
    NumericVector RollSumh = yRollsum[h];
    for(int j = 0; j < nm; ++j) {
      IntegerVector idx = IntegerVector::create(j, j+h);
      NumericVector RollSumSel = RollSumh[idx];
      rollacfs(j,h) = ProdPairs(j,h) - yRollmean[j]*sum(RollSumSel)/m + (m-h)*pow(yRollmean[j],2)/m;
    }
  }
  return rollacfs;
}