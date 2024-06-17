// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <RcppEigen.h>

// Returns A %*% B in the form of a matrix
// [[Rcpp::export]]
SEXP eigenMatProduct(Eigen::MatrixXd A, Eigen::MatrixXd B){
  return Rcpp::wrap(A * B);
}

// Returns A %*% t(B) in the form of a matrix
// [[Rcpp::export]]
SEXP eigenMatTcrossprod(Eigen::MatrixXd A, Eigen::MatrixXd B){
  return Rcpp::wrap(A * B.transpose());
}

// Returns t(A) %*% B in the form of a matrix
// [[Rcpp::export]]
SEXP eigenMatCrossprod(Eigen::MatrixXd A, Eigen::MatrixXd B){
  return Rcpp::wrap(A.transpose() * B);
}

// Returns A %*% b in the form of a column vector
// [[Rcpp::export]]
SEXP eigenMatVecProduct(Eigen::MatrixXd A, Eigen::VectorXd b){
  return Rcpp::wrap(A * b);
}

// Returns the transpose of t(a) %*% B in the form of a column vector
// [[Rcpp::export]]
SEXP eigenVecMatProduct(Eigen::VectorXd a, Eigen::MatrixXd B){
  return Rcpp::wrap(B.transpose() * a);
}

// Returns the transpose of a %*% t(b) in the form of a matrix
// [[Rcpp::export]]
SEXP eigenVecVecProduct(Eigen::VectorXd a, Eigen::VectorXd b){
  return Rcpp::wrap(a * b.transpose());
}

// Returns the row-wise L2 norm of a matrix
// [[Rcpp::export]]
SEXP eigenRowL2(Eigen::MatrixXd A){
  return Rcpp::wrap(A.rowwise().squaredNorm());
}

// Returns the column-wise L2 norm of a matrix
// [[Rcpp::export]]
SEXP eigenColL2(Eigen::MatrixXd A){
  return Rcpp::wrap(A.colwise().squaredNorm());
}

// Returns the column-wise L1 norm of a matrix
// [[Rcpp::export]]
SEXP eigenColL1(Eigen::MatrixXd A){
  return Rcpp::wrap(A.colwise().lpNorm<1>());
}

// Multiplies each row of a matrix by a vector (coefficient-wise)
// [[Rcpp::export]]
SEXP eigenRowMult( Eigen::MatrixXd A, Eigen::ArrayXd b ){
  return Rcpp::wrap( (A.array().rowwise()*b.transpose()).matrix() );
}

/*******************************************************/
// Returns the transpose of a %*% t(b) in the form of a matrix
// [[Rcpp::export]]
SEXP solveB(
    Eigen::MatrixXd Y,
    Eigen::MatrixXd Z,
    Eigen::MatrixXd diagK,
    Eigen::VectorXd s, Eigen::VectorXd o ) {
  
  return Rcpp::wrap(
    (Z.transpose()*Z).ldlt().solve(diagK) *
      ( Z.transpose() * ((Y.rowwise()-s.transpose()).colwise()-o) )
    );  
}



// Solve Z
// [[Rcpp::export]]
SEXP solveZ(
	Eigen::MatrixXd Y,
    Eigen::VectorXd s, Eigen::VectorXd o,
    Eigen::MatrixXd B,
    Eigen::MatrixXd diagK,
    Eigen::MatrixXd diagLambda) {
  
  return Rcpp::wrap(
  	(((Y.rowwise()-s.transpose()).colwise()-o)*B.transpose()) *
  		(B*B.transpose()+diagLambda).ldlt().solve(diagK) );
}


// Solve si
// [[Rcpp::export]]
SEXP solveS(
    Eigen::MatrixXd Y,
    Eigen::MatrixXd ZDB,
    Eigen::VectorXd j,
    Eigen::VectorXd o,
    double J ) {
  
  return Rcpp::wrap(
    ( ((Y-ZDB).colwise()-o).transpose()*j ) / J );
}


// Returns the row sum of the residual of Yi after taking into account everything other than o
// [[Rcpp::export]]
SEXP Y_resO_rowSum(
    Eigen::MatrixXd Y,
    Eigen::MatrixXd ZDB,
    Eigen::VectorXd n,
    Eigen::VectorXd s ) {
  
  return Rcpp::wrap( ((Y-ZDB).rowwise()-s.transpose())*n );
}



// Solves Yi as log of unobserved expression
// [[Rcpp::export]]
SEXP solveY(
    Eigen::SparseMatrix<double> M,
    Eigen::ArrayXXd Y,
    Eigen::ArrayXXd ZDB,
    Eigen::ArrayXd s, Eigen::ArrayXd o,
    Eigen::ArrayXd sigma2 ) {
  
  // Calculate the estimate of Yi based on model parameters
  Eigen::ArrayXXd Y_hat = (ZDB.rowwise()+s.transpose()).colwise()+o;
  // Use Halley's method to obtain an approximation of Yi.
  //   Note that with a reasonable starting point, Halley's method
  //   leads to reasonable estimates in one iteration, which should improve
  //   in subsequent iterations
  Eigen::ArrayXXd fpp = Y.exp().rowwise()*sigma2.transpose();
  Eigen::ArrayXXd fp = fpp + 1;
  Eigen::ArrayXXd f = fpp + Y - Y_hat -
    Eigen::ArrayXXd(Eigen::MatrixXd(M)).rowwise()*sigma2.transpose();
  // check for overshoot, and correct
  Eigen::ArrayXXd solution = Y - 2*f*fp/(2*fp.square()-f*fpp);
  int overshoot = 0, undershoot = 0;
  for( int i = 0; i < solution.rows(); i ++ ) {
    for( int j = 0; j < solution.cols(); j ++ ) {
      // the solution must be between log(M) and Y_hat
      if( solution(i,j) > Y_hat(i,j) &&
          solution(i,j) > log(M.coeff(i,j)+1e-50) ) {
          solution(i,j) = ( log(M.coeff(i,j)+1e-50) + Y_hat(i,j) ) / 2;
          overshoot ++;
      } else if( solution(i,j) < Y_hat(i,j) &&
        solution(i,j) < log(M.coeff(i,j)+1e-50) ) {
        solution(i,j) = ( log(M.coeff(i,j)+1e-50) + Y_hat(i,j) ) / 2;
        undershoot ++;
      }
    }
  }
  if( overshoot+undershoot > 0 ) {
    Rcpp::Rcerr << "! " << overshoot << " overshoots and " <<
      undershoot << " undershoots resolved." <<
      std::endl;
  }
  
  return Rcpp::wrap( solution );

}



// Returns the L2 norm of a vector
// [[Rcpp::export]]
double vecL2_noPrior( Eigen::VectorXd x ) {
  return x.squaredNorm();
}
// Returns the L2 norm of a vector after subtracting a prior (squared Euclidean distance)
// [[Rcpp::export]]
double vecL2_wPrior( Eigen::VectorXd x, Eigen::VectorXd prior ) {
  return (x-prior).squaredNorm();
}
// Returns the L2 norm of a matrix
// [[Rcpp::export]]
double matL2_noPrior( Eigen::MatrixXd X ) {
  return X.squaredNorm();
}

// Returns the L2 norm of a matrix after subtracting a prior 
// [[Rcpp::export]]
double matL2_wPrior( Eigen::MatrixXd X, Eigen::MatrixXd A, Eigen::MatrixXd B ) {
  return (X-A*B).squaredNorm();
}

// // Returns the SSE (sum-squared-error) of Yi
// // [[Rcpp::export]]
// double Yi_SSE(
//     Eigen::MatrixXd Yi,
//     Eigen::MatrixXd ZDBi, Eigen::MatrixXd QiDBi,
//     Eigen::VectorXd si, Eigen::VectorXd o, Eigen::VectorXd oi ) {
//   
//   return ( ((Yi-ZDBi-QiDBi).rowwise()-si.transpose()).colwise()-(o+oi) ).squaredNorm();
// }

// double Y_SSE_M(
//     Eigen::MatrixXd Y,
//     Eigen::MatrixXd ZDB,
//     Eigen::VectorXd s, Eigen::VectorXd o,
//     double sigma2 ) {
//   
//   return ( ((Y-ZDB).rowwise()-s.transpose()).colwise()-o ).squaredNorm() +
//     (1/(Eigen::ArrayXXd(Y).exp()+1/sigma2)).sum();
// }

// Returns the variance of Yi-Yhat
// [[Rcpp::export]]
SEXP solveSigma2_M(
    Eigen::ArrayXXd Y,
    Eigen::ArrayXXd ZDB,
    Eigen::ArrayXd s, Eigen::ArrayXd o,
    double J,
    Eigen::ArrayXd sigma2 ) {
  
  return Rcpp::wrap(
    ( ( ((Y-ZDB).rowwise()-s.transpose()).colwise()-o ).square().colwise().sum() +
      (1/(Eigen::ArrayXXd(Y).exp().rowwise()+sigma2.inverse().transpose())).colwise().sum() ) / J );
  
}


// Returns the variance of Yi-Yhat
// [[Rcpp::export]]
SEXP Yi_var( Eigen::MatrixXd Yi, double sigma2 ) {
  return Rcpp::wrap( 1/(Eigen::ArrayXXd(Yi).exp()+1/sigma2) );
}

// Returns the variance of Yi-Yhat
// [[Rcpp::export]]
SEXP Yi_var_paired(
    Eigen::MatrixXd Yi,
    Eigen::SparseMatrix<double> M1i, Eigen::SparseMatrix<double> M2i,
    double sigma2 ) {
  
  Eigen::ArrayXXd expYi = (-Eigen::ArrayXXd(Yi).abs()).exp();
  Eigen::ArrayXXd M = Eigen::ArrayXXd(M1i+M2i);
  
  return Rcpp::wrap( 1/(M*expYi/(1+expYi).square()+1/sigma2) );
}

// Returns the variance of Yi-Yhat
// [[Rcpp::export]]
SEXP predict_Yhat(
    Eigen::MatrixXd ZDBi, Eigen::MatrixXd QiDBi,
    Eigen::VectorXd si, Eigen::VectorXd o, Eigen::VectorXd oi ) { 
  return Rcpp::wrap( ((ZDBi+QiDBi).rowwise()+si.transpose()).colwise()+(o+oi) );
}

// Returns the root mean squared difference of two matrices
// [[Rcpp::export]]
double matRMSD( Eigen::MatrixXd A, Eigen::MatrixXd B ) {
  return sqrt( (A-B).squaredNorm()/A.size() );
}
// Returns the root mean squared difference of two vectors
// [[Rcpp::export]]
double vecRMSD( Eigen::VectorXd A, Eigen::VectorXd B ) {
  return sqrt( (A-B).squaredNorm()/A.size() );
}


// Performs cbind of two matrices
// Benchmarking shows that R cbind is twice as fast. So, let's continue using that.
// SEXP eigenCBind( Eigen::MatrixXd A, Eigen::MatrixXd B ) {
//   
//   Eigen::MatrixXd AB(A.rows(),A.cols()+B.cols()); AB << A,B;
//   return Rcpp::wrap( AB );
// }

// Performs cbind of two matrices
// Benchmarking shows that R cbind is still slightly faster. So, let's continue using that.
// Rcpp::NumericMatrix RcppCBind(Rcpp::NumericMatrix a, Rcpp::NumericMatrix b) {
//   int acoln = a.ncol();
//   int bcoln = b.ncol();
//   Rcpp::NumericMatrix out = Rcpp::no_init_matrix(a.nrow(), acoln + bcoln);
//   for (int j = 0; j < acoln + bcoln; j++) {
//     if (j < acoln) {
//       out(Rcpp::_, j) = a(Rcpp::_, j);
//     } else {
//       out(Rcpp::_, j) = b(Rcpp::_, j - acoln);
//     }
//   }
//   return out;
// }
