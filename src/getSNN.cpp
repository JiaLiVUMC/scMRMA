#include <RcppEigen.h>
//#include <cmath>
//#include <unordered_map>
//#include <fstream>
//#include <string>
//#include <iomanip>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::Triplet<double> T;
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

//[[Rcpp::export]]
NumericMatrix getSNN2(Eigen::MatrixXd knnMatrix, double prune) {
  std::vector<T> tripletList;
  int k = knnMatrix.cols();
  tripletList.reserve(knnMatrix.rows() * knnMatrix.cols());
  for(int j=0; j<knnMatrix.cols(); ++j){
    for(int i=0; i<knnMatrix.rows(); ++i) {
      tripletList.push_back(T(i, knnMatrix(i, j) - 1, 1));
    }
  }
  Eigen::SparseMatrix<double> SNNLink(knnMatrix.rows(), knnMatrix.rows());
  SNNLink.setFromTriplets(tripletList.begin(), tripletList.end());
  SNNLink = SNNLink * (SNNLink.transpose());
  for (int i=0; i < SNNLink.outerSize(); ++i){
    for (Eigen::SparseMatrix<double>::InnerIterator it(SNNLink, i); it; ++it){
      it.valueRef() = it.value()/(k + (k - it.value()));
      if(it.value() < prune){
        it.valueRef() = 0;
      }
    }
  }
  SNNLink.prune(0.0); 
  
  int r=0;
  for (int m=0; m < SNNLink.outerSize(); ++m){
  for (Eigen::SparseMatrix<double>::InnerIterator it(SNNLink,m); it; ++it){
    r++;
    }
  }

  NumericMatrix SNNmatrix(r, 3);
  int p = 0;
  for (int n=0; n < SNNLink.outerSize(); ++n){
    for (Eigen::SparseMatrix<double>::InnerIterator it(SNNLink,n); it; ++it){
      SNNmatrix(p,0) = it.row();
      SNNmatrix(p,1) = it.col();
      SNNmatrix(p,2) = it.value();
      p++;
    }
  }
  return SNNmatrix;
}