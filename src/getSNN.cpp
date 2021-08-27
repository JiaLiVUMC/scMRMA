#include <RcppEigen.h>
//#include <cmath>
//#include <unordered_map>
//#include <fstream>
//#include <string>
//#include <iomanip>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

 
typedef Eigen::Triplet<double> T;

//' SNN matrix
//' @export
//'
//[[Rcpp::export]]
Eigen::SparseMatrix<double> getSNN(Eigen::MatrixXd knnMatrix, double prune) {
  std::vector<T> tripletList;
  int k = knnMatrix.cols();
  tripletList.reserve(knnMatrix.rows() * knnMatrix.cols());
  for(int j=0; j<knnMatrix.cols(); ++j){
    for(int i=0; i<knnMatrix.rows(); ++i) {
      tripletList.push_back(T(i, knnMatrix(i, j) - 1, 1));
    }
  }
  Eigen::SparseMatrix<double> SNN(knnMatrix.rows(), knnMatrix.rows());
  SNN.setFromTriplets(tripletList.begin(), tripletList.end());
  SNN = SNN * (SNN.transpose());
  for (int i=0; i < SNN.outerSize(); ++i){
    for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, i); it; ++it){
      it.valueRef() = it.value()/(k + (k - it.value()));
      if(it.value() < prune){
        it.valueRef() = 0;
      }
    }
  }
  SNN.prune(0.0); 
  return SNN;
}
