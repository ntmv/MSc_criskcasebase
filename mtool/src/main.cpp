// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
#include <iostream>
#include <vector> 
#include <string>
#include <cmath>       
#include <numeric>        
#include <unordered_map>  
#include <limits>
#include "spams.h"
// we only include RcppArmadillo.h which pulls Rcpp.h in for us


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operation on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//


// [[Rcpp::export]]
arma::vec grad_ls_loss(
    arma::rowvec& x,
    double& y,
    arma::vec& param,
    int& p) {
  arma::vec grad(p);
  grad =  arma::vectorise(x) * (arma::dot(x, param) - y);
  return grad;
}

// [[Rcpp::export]]
arma::vec grad_logistic_loss(
    arma::rowvec& x,
    double& y,
    arma::vec& param,
    int& p) {
  arma::vec grad(p);
  double sig;
  sig = 1.0/( 1.0 + exp(-arma::dot(x, param)) );
  grad =  arma::vectorise(x) * (sig - y);
  return grad;
}


// [[Rcpp::export]]
arma::mat grad_multinom_loss(
    arma::rowvec& x,
    int& y,
    int& K,
    double& offset,
    arma::mat& param,
    int& p) {
  arma::mat grad(p, K);
  arma::vec pi(K);
  
  for (int i = 0; i < K; i++) {
    pi(i) = exp(arma::dot(x, param.col(i)) + offset);
  }
  pi = pi/(arma::sum(pi) + 1.0);
  
  if (y == 0) {
    for (int k = 0; k < K; k++) {
      grad.col(k) = pi(k) * arma::vectorise(x);
    }
  } else {
    for (int k = 0; k < K; k++) {
      if (k == y - 1) {
        grad.col(k) = (pi(k) - 1) * arma::vectorise(x);
      } else {
        grad.col(k) = pi(k) * arma::vectorise(x);
      }
    }
  }
  return grad;
}



void proximalFlat(
    arma::mat& U,
    int& p,              // # of rows of U
    int& K,              // # of columns of U
    std::string& regul,
    Rcpp::IntegerVector& grp_id, // grp is a vector for proximalFlat
    int num_threads,
    double lam1,
    double lam2 = 0.0,
    double lam3 = 0.0,
    bool intercept = false,
    bool resetflow = false,
    bool verbose = false,
    bool pos = false,
    bool clever = true,
    bool eval = true,
    int size_group = 1,
    bool transpose = false) {
  // dimensions
  // int p = U.n_rows;
  // int K = U.n_cols;
  int pK = int (p*K);
  
  // read in U and convert to spams::Matrix<double> alpha0
  double* ptrU = new double[pK];
  for (int r = 0; r < p; r++) {
    for (int c = 0; c < K; c++) {
      ptrU[c * p + r] = U(r, c);
    }
  }
  Matrix<double> alpha0(ptrU,p,K);
  
  
  
  // read in regul and convert to char
  int l = regul.length();
  char* name_regul = new char[l];
  for (int i = 0; i < l+1; i++) {
    name_regul[i] = regul[i];
  }
  
  
  
  // Initialize alpha - proximal operator
  Matrix<double> alpha(p,K);
  alpha.setZeros();
  
  // read in grp_id and convert to spams::Vector<int> groups
  int* ptrG = new int[p];
  for (int i = 0; i < p; i++) {
    ptrG[i] = grp_id(i);
  }
  
  
  
  Vector<int> groups(ptrG, p);
  
  _proximalFlat(&alpha0,&alpha,&groups,num_threads,
                lam1,lam2,lam3,intercept,
                resetflow,name_regul,verbose,pos,
                clever,eval,size_group,transpose);
  
  // put updated alpha back into U
  for (int r = 0; r < p; r++) {
    for (int c = 0; c < K; c++) {
      U(r, c) = alpha[p * c + r];
    }
  }
  
  
  // free the dynamic memory
  delete[] ptrU;
  delete[] name_regul;
  delete[] ptrG;
}





void proximalGraph(
    arma::mat& U,
    int& p,              // # of rows of U
    int& K,              // # of columns of U
    std::string& regul,
    arma::mat& grp,      // grp is a matrix for proximalGraph and Tree
    arma::mat& grpV,
    Rcpp::NumericVector& etaG,
    int num_threads,
    double lam1,
    double lam2 = 0.0,
    double lam3 = 0.0,
    bool intercept = false,
    bool resetflow = false,
    bool verbose = false,
    bool pos = false,
    bool clever = true,
    bool eval = true,
    int size_group = 1,
    bool transpose = false) {
  
  // U dimensions
  // int p = U.n_rows;
  // int K = U.n_cols;
  int pK = int (p*K);
  
  // read in U and convert to spams::Matrix<double> alpha0
  double* ptrU = new double[pK];
  for (int r = 0; r < p; r++) {
    for (int c = 0; c < K; c++) {
      ptrU[c * p + r] = U(r, c);
    }
  }
  Matrix<double> alpha0(ptrU,p,K);
  
  
  
  // grp dimensions
  int gr = grp.n_rows;
  int gc = grp.n_cols;
  int grc = int (gr * gc);
  
  // read in grp and convert to spams::Matrix<bool> grp_dense
  // then to spams::SpMatrix<bool> groups
  bool* ptrG = new bool[grc];
  for (int r = 0; r < gr; r++) {
    for (int c = 0; c < gc; c++) {
      ptrG[c * gr + r] = (grp(r, c) != 0.0);
    }
  }
  Matrix<bool> grp_dense(ptrG, gr, gc);
  SpMatrix<bool> groups;
  grp_dense.toSparse(groups);
  
  
  
  // grpV dimensions
  int gvr = grpV.n_rows;
  int gvc = grpV.n_cols;
  int gvrc = int (gvr * gvc);
  
  // read in grpV and convert to spams::Matrix<bool> grpV_dense
  // then to spams::SpMatrix<bool> groups_var
  bool* ptrGV = new bool[gvrc];
  for (int r = 0; r < gvr; r++) {
    for (int c = 0; c < gvc; c++) {
      ptrGV[c * gvr + r] = (grpV(r, c) != 0.0);
    }
  }
  Matrix<bool> grpV_dense(ptrGV, gvr, gvc);
  SpMatrix<bool> groups_var;
  grpV_dense.toSparse(groups_var);
  
  
  
  // read in etaG and convert to spams::Vector<int> eta_g
  int n_etaG = etaG.length();
  double* ptrEG = new double[n_etaG];
  for (int i = 0; i < n_etaG; i++) {
    ptrEG[i] = etaG(i);
  }
  Vector<double> eta_g(ptrEG, n_etaG);
  
  
  
  // read in regul and convert to char
  int l = regul.length();
  char* name_regul = new char[l];
  for (int i = 0; i < l+1; i++) {
    name_regul[i] = regul[i];
  }
  
  
  
  // Initialize alpha - proximal operator
  Matrix<double> alpha(p,K);
  alpha.setZeros();
  
  
  // call _proximalGraph
  _proximalGraph(&alpha0, &alpha,
                 &eta_g, &groups, &groups_var,
                 num_threads, lam1, lam2, lam3,
                 intercept, resetflow, name_regul,
                 verbose, pos, clever, eval,
                 size_group, transpose);
  
  
  
  // put updated alpha back into U
  for (int r = 0; r < p; r++) {
    for (int c = 0; c < K; c++) {
      U(r, c) = alpha[c * p + r];
    }
  }
  
  // free the dynamic memory
  delete[] ptrU;
  delete[] ptrG;
  delete[] ptrGV;
  delete[] ptrEG;
  delete[] name_regul;
}




void proximalTree(
    arma::mat& U,
    int& p,              // # of rows of U
    int& K,              // # of columns of U
    std::string& regul,
    arma::mat& grp,      // grp is a matrix for proximalGraph and Tree
    Rcpp::NumericVector& etaG,
    Rcpp::IntegerVector& own_var,
    Rcpp::IntegerVector& N_own_var,
    int num_threads,
    double lam1,
    double lam2 = 0.0,
    double lam3 = 0.0,
    bool intercept = false,
    bool resetflow = false,
    bool verbose = false,
    bool pos = false,
    bool clever = true,
    bool eval = true,
    int size_group = 1,
    bool transpose = false) {
  
  
  
  // U dimensions
  // int p = U.n_rows;
  // int K = U.n_cols;
  int pK = int (p*K);
  
  // read in U and convert to spams::Matrix<double> alpha0
  double* ptrU = new double[pK];
  for (int r = 0; r < p; r++) {
    for (int c = 0; c < K; c++) {
      ptrU[c * p + r] = U(r, c);
    }
  }
  Matrix<double> alpha0(ptrU,p,K);
  
  
  
  // grp dimensions
  int gr = grp.n_rows;
  int gc = grp.n_cols;
  int grc = int (gr * gc);
  
  // read in grp and convert to spams::Matrix<bool> grp_dense
  // then to spams::SpMatrix<bool> groups
  bool* ptrG = new bool[grc];
  for (int r = 0; r < gr; r++) {
    for (int c = 0; c < gc; c++) {
      ptrG[c * gr + r] = ( grp(r, c) != 0.0 );
    }
  }
  Matrix<bool> grp_dense(ptrG, gr, gc);
  SpMatrix<bool> groups;
  grp_dense.toSparse(groups);
  
  
  
  // read in etaG and convert to spams::Vector<double> eta_g
  int n_etaG = etaG.length();
  double* ptrEG = new double[n_etaG];
  for (int i = 0; i < n_etaG; i++) {
    ptrEG[i] = etaG(i);
  }
  Vector<double> eta_g(ptrEG, n_etaG);
  
  
  
  // read in own_var and convert to spams::Vector<int> own_variables
  int n_ov = own_var.length();
  int* ptrOV = new int[n_ov];
  for (int i = 0; i < n_ov; i++) {
    ptrOV[i] = own_var(i);
  }
  Vector<int> own_variables(ptrOV, n_ov);
  
  
  
  // read in N_own_var and convert to spams::Vector<int> N_own_variables
  int n_Nov = N_own_var.length();
  int* ptrNOV = new int[n_Nov];
  for (int i = 0; i < n_Nov; i++) {
    ptrNOV[i] = N_own_var(i);
  }
  Vector<int> N_own_variables(ptrNOV, n_Nov);
  
  
  
  // read in regul and convert to char
  int l = regul.length();
  char* name_regul = new char[l];
  for (int i = 0; i < l+1; i++) {
    name_regul[i] = regul[i];
  }
  
  
  
  // Initialize alpha - proximal operator
  Matrix<double> alpha(p,K);
  alpha.setZeros();
  
  
  
  // call _proximalTree
  _proximalTree(&alpha0, &alpha,
                &eta_g, &groups,
                &own_variables, &N_own_variables,
                num_threads, lam1, lam2, lam3,
                intercept, resetflow, name_regul,
                verbose, pos, clever, eval,
                size_group, transpose);
  
  
  
  // put updated alpha back into U
  for (int r = 0; r < p; r++) {
    for (int c = 0; c < K; c++) {
      U(r, c) = alpha[c * p + r];
    }
  }
  
  // free the dynamic memory
  delete[] ptrU;
  delete[] ptrG;
  delete[] ptrEG;
  delete[] ptrOV;
  delete[] ptrNOV;
  delete[] name_regul;
}





// [[Rcpp::export]]
Rcpp::List mtool(
    arma::mat X,            // input
    arma::mat Y,            // outcome
    arma::vec wt,
    int K,                  // number of tasks
    int reg_p,              // number of regularized variables
    arma::vec nk_vec,
    arma::vec task_rowid,
    int loss,               // 1 - least squares; 2 - logistic regression; 3 - Cox
    int penalty,            // 1 - proximalFlat; 2 - proximalGraph; 3 - proximalTree
    std::string regul,
    bool transpose,
    Rcpp::IntegerVector grp_id,
    Rcpp::NumericVector etaG,
    arma::mat grp,
    arma::mat grpV,
    Rcpp::IntegerVector own_var,
    Rcpp::IntegerVector N_own_var,
    double lam1,
    double lam2,
    double lam3,
    double learning_rate,
    double tolerance,
    int niter_inner,
    int maxit,
    int ncores) {
  
  // initialize param
  int p = X.n_cols;
  
  // indexing, stochastic sampling, etc..
  int init_k;   // index of the first obs in task k
  int n_k;      // number of obs in task k
  int id_ik;    // index of the i-th obs in task k
  int index;    // index of the stochastic sample
  arma::rowvec x_sample(p);
  double y_sample;
  double wt_sample;
  
  // gradient and related
  arma::mat grad(p, K);
  arma::vec grad_k(p);         // k-th column of grad
  arma::vec temp1(p);
  arma::vec temp2(p);
  
  // parameter and related
  arma::mat param(p, K);
  param.zeros();
  arma::vec param_k(p);     // k-th column of param
  arma::mat param_old(p, K);
  arma::vec param_old_k(p);
  
  arma::mat param_reg(reg_p, K);   // regularized coefficients
  arma::mat param_t(K, reg_p);      // regularized coefficients
  
  // convergence related
  // arma::mat param_update(p, K);
  double diff;
  int counter_outer = 0;
  
  // compute mu: mean gradient at param_old
  while (true) {
    param_old = param;
    init_k = 0;
    
    for (int k = 0; k < K; k++) {
      n_k = nk_vec(k);
      param_old_k = param_old.col(k);
      grad_k.zeros();
      
      for (int i = 0; i < n_k; i++) {
        id_ik = task_rowid(init_k + i) - 1;
        x_sample = X.row(id_ik);
        y_sample = Y(id_ik, k);
        wt_sample = wt(id_ik);
        if (loss == 1) {
          grad_k = grad_k + wt_sample * grad_ls_loss(x_sample, y_sample, param_old_k, p)/n_k;
        }
        if (loss == 2) {
          grad_k = grad_k + wt_sample * grad_logistic_loss(x_sample, y_sample, param_old_k, p)/n_k;
        }
      }
      grad.col(k) = grad_k;
      init_k += n_k;
    }
    
    // inner loop
    for (int i = 0; i < niter_inner; ++i) {
      init_k = 0;
      
      for (int k = 0; k < K; k++) {
        n_k = nk_vec(k);
        index = arma::randi(arma::distr_param(0, n_k - 1));
        id_ik = task_rowid(init_k + index) - 1;
        
        x_sample = X.row(id_ik);
        y_sample = Y(id_ik, k);
        wt_sample = wt(id_ik);
        
        param_k = param.col(k);
        param_old_k = param_old.col(k);
        
        if (loss == 1) {
          temp1 = grad_ls_loss(x_sample, y_sample, param_k, p);
          temp2 = grad_ls_loss(x_sample, y_sample, param_old_k, p);
        }
        
        if (loss == 2) {
          temp1 = grad_logistic_loss(x_sample, y_sample, param_k, p);
          temp2 = grad_logistic_loss(x_sample, y_sample, param_old_k, p);
        }
        
        param_k = param_k - learning_rate * (wt_sample * (temp1 - temp2) + grad.col(k));
        
        param.col(k) = param_k;
        
        init_k += n_k;
      }
      
      // extract only variables involved in the penalization
      param_reg = param.head_rows(reg_p);
      
      // call proximal function
      if (transpose) {
        param_t = param_reg.t();
        
        if (penalty == 1) {
          proximalFlat(param_t, K, reg_p, regul, grp_id, ncores, lam1 * learning_rate, lam2 * learning_rate, lam3 * learning_rate);
        }
        
        if (penalty == 2) {
          proximalGraph(param_t, K, reg_p, regul, grp, grpV, etaG, ncores, lam1 * learning_rate, lam2 * learning_rate);
        }
        
        if (penalty == 3) {
          proximalTree(param_t, K, reg_p, regul, grp, etaG, own_var, N_own_var, ncores, lam1 * learning_rate, lam2 * learning_rate);
        }
        
        param_reg = param_t.t();
      } else {
        if (penalty == 1) {
          proximalFlat(param_reg, reg_p, K, regul, grp_id, ncores, lam1 * learning_rate, lam2 * learning_rate, lam3 * learning_rate);
        }
        
        if (penalty == 2) {
          proximalGraph(param_reg, reg_p, K, regul, grp, grpV, etaG, ncores, lam1 * learning_rate, lam2 * learning_rate);
        }
        
        if (penalty == 3) {
          proximalTree(param_reg, reg_p, K, regul, grp, etaG, own_var, N_own_var, ncores, lam1 * learning_rate, lam2 * learning_rate);
        }
      }
      
      param.head_rows(reg_p) = param_reg;
      
    }
    
    counter_outer += 1;
    Rcpp::Rcout << "\n Iteration " << counter_outer <<"\n";
    
    diff = arma::norm(param - param_old, "fro");
    diff = diff/(p*K);
    Rcpp::Rcout << "Mean Frobenius norm of coefficient update \n" << diff <<"\n";
    
    if (diff < tolerance || counter_outer>= maxit) {
      break;
    }
  }
  
  arma::sp_mat param_sp(param);
  
  Rcpp::List result = Rcpp::List::create(Rcpp::Named("Estimates")                     = param,
                                         Rcpp::Named("Sparse Estimates")              = param_sp);
  return result;
}



// [[Rcpp::export]]
Rcpp::List MultinomLogistic(
    arma::mat X,
    arma::vec Y,
    arma::vec offset,
    int K,
    int reg_p,              // number of regularized variables
    int penalty,            // 1 - proximalFlat; 2 - proximalGraph; 3 - proximalTree
    std::string regul,
    bool transpose,
    Rcpp::IntegerVector grp_id,
    Rcpp::NumericVector etaG,
    arma::mat grp,
    arma::mat grpV,
    Rcpp::IntegerVector own_var,
    Rcpp::IntegerVector N_own_var,
    double lam1,
    double lam2,
    double lam3,
    double learning_rate,
    double tolerance,
    int niter_inner,
    int maxit,
    int ncores) {
  
  int p = X.n_cols;
  int n = X.n_rows;
  
  // indexing, stochastic sampling, etc..
  int index;    // index of the stochastic sample
  arma::rowvec x_sample(p);
  int y_sample;
  double o_sample; // o for offset
  
  // gradient and related
  arma::mat grad(p, K);
  arma::mat temp1(p, K);
  arma::mat temp2(p, K);
  
  // parameter and related
  arma::mat param(p, K);
  param.zeros();
  arma::mat param_old(p, K);

  arma::mat param_reg(reg_p, K);   // regularized coefficients
  arma::mat param_t(K, reg_p);      // regularized coefficients
  
  // convergence related
  double diff;
  int counter_outer = 0;
  
  // compute mu: mean gradient at param_old
  while (true) {
    param_old = param;
    grad.zeros();
    
    for (int i = 0; i < n; i++) {
      x_sample = X.row(i);
      y_sample = Y(i);
      o_sample = offset(i);
      grad = grad + grad_multinom_loss(x_sample, y_sample, K, o_sample, param_old, p)/n;
    }

    // inner loop
    for (int i = 0; i < niter_inner; ++i) {

      index = arma::randi(arma::distr_param(0, n - 1));

      x_sample = X.row(index);
      y_sample = Y(index);
      o_sample = offset(index);

      temp1 = grad_multinom_loss(x_sample, y_sample, K, o_sample, param, p);
      temp2 = grad_multinom_loss(x_sample, y_sample, K, o_sample, param_old, p);

      param = param - learning_rate * (temp1 - temp2 + grad);

      // extract only variables involved in the penalization
      param_reg = param.head_rows(reg_p);

      // call proximal function
      if (transpose) {
        param_t = param_reg.t();
        
        if (penalty == 1) {
          proximalFlat(param_t, K, reg_p, regul, grp_id, ncores, lam1 * learning_rate, lam2 * learning_rate, lam3 * learning_rate);
        }
        
        if (penalty == 2) {
          proximalGraph(param_t, K, reg_p, regul, grp, grpV, etaG, ncores, lam1 * learning_rate, lam2 * learning_rate);
        }
        
        if (penalty == 3) {
          proximalTree(param_t, K, reg_p, regul, grp, etaG, own_var, N_own_var, ncores, lam1 * learning_rate, lam2 * learning_rate);
        }
        
        param_reg = param_t.t();
      } else {
        if (penalty == 1) {
          proximalFlat(param_reg, reg_p, K, regul, grp_id, ncores, lam1 * learning_rate, lam2 * learning_rate, lam3 * learning_rate);
        }
        
        if (penalty == 2) {
          proximalGraph(param_reg, reg_p, K, regul, grp, grpV, etaG, ncores, lam1 * learning_rate, lam2 * learning_rate);
        }
        
        if (penalty == 3) {
          proximalTree(param_reg, reg_p, K, regul, grp, etaG, own_var, N_own_var, ncores, lam1 * learning_rate, lam2 * learning_rate);
        }
      }
      
      param.head_rows(reg_p) = param_reg;
    }
    
    counter_outer += 1;
    Rcpp::Rcout << "\n Iteration " << counter_outer <<"\n";
    
    diff = arma::norm(param - param_old, "fro");
    diff = diff/(p*K);
    Rcpp::Rcout << "Frobenius norm of coefficient update \n" << diff <<"\n";
    
    if (diff < tolerance || counter_outer>= maxit) {
      break;
    }
  }
  
  arma::sp_mat param_sp(param);
  
  Rcpp::List result = Rcpp::List::create(Rcpp::Named("Estimates")                     = param,
                                         Rcpp::Named("Sparse Estimates")              = param_sp);
  return result;
}















void proximalFlat2(
        arma::mat& U,
        int& p,              
        int& K,           
        std::string& regul,
        Rcpp::IntegerVector& grp_id,
        int num_threads,
        double lam1,
        double lam2 = 0.0,
        double lam3 = 0.0,
        bool intercept = false,
        bool resetflow = false,
        bool verbose = false,
        bool pos = false,
        bool clever = true,
        bool eval = true,
        int size_group = 1,
        bool transpose = false) {
    // dimensions
    // int p = U.n_rows;
    // int K = U.n_cols;
    // int pK = int (p*K);
    
    // read in U and convert to spams::Matrix<double> alpha0
    // double* ptrU = new double[pK];
    // for (int r = 0; r < p; r++) {
    //     for (int c = 0; c < K; c++) {
    //         ptrU[c * p + r] = U(r, c);
    //     }
    // }
    // Matrix<double> alpha0(ptrU,p,K);
    Matrix<double> alpha0(const_cast<double*>(U.memptr()), p, K);
    
    
    
    // read in regul and convert to char
    // int l = regul.length();
    // char* name_regul = new char[l];
    // for (int i = 0; i < l+1; i++) {
    //     name_regul[i] = regul[i];
    // }
    std::vector<char> name_regul_buf(regul.begin(), regul.end()); 
    name_regul_buf.push_back('\0'); 
    char* name_regul_ptr = name_regul_buf.data();
    
    
    
    // Initialize alpha - proximal operator
    Matrix<double> alpha(p,K);
    alpha.setZeros();
    
    // read in grp_id and convert to spams::Vector<int> groups
    // int* ptrG = new int[p];
    // for (int i = 0; i < p; i++) {
    //     ptrG[i] = grp_id(i);
    // }
    // Vector<int> groups(ptrG, p);
    
    Vector<int> groups(const_cast<int*>(grp_id.begin()), p); 
    
    
    
    
    _proximalFlat(&alpha0,&alpha,&groups,num_threads,
                   lam1,lam2,lam3,intercept,
                   resetflow,name_regul_ptr,verbose,pos,
                   clever,eval,size_group,transpose);
    
    // put updated alpha back into U
    for (int r = 0; r < p; r++) {
        for (int c = 0; c < K; c++) {
            U(r, c) = alpha[p * c + r];
        }
    }
    
    
    // free the dynamic memory
    // delete[] ptrU;
    // delete[] name_regul;
    // delete[] ptrG;
}


void proximalFlat2_nospams(
        arma::mat& U,             
        int p,                    
        int K,                    
        const std::string& regul, 
        const Rcpp::IntegerVector& grp_id, 
        double lam1,              
        double lam2 = 0.0,        
        double lam3 = 0.0,        
        bool intercept = false,   
        bool resetflow = false,   
        bool verbose = false,     
        bool pos = false,         
        bool clever = true,       
        bool eval = true,         
        int size_group = 1,       
        bool transpose = false) { 
    
    if (lam1 < 0.0 || lam2 < 0.0 || lam3 < 0.0) {
        Rcpp::stop("Lambda penalty parameters cannot be negative.");
    }
    
    
    if (regul == "L1") {
        // Elementwise L1 Soft-thresholding: sign(u) * max(|u| - lam1, 0)
        if (lam1 == 0.0) return; // No thresholding needed
        U = arma::sign(U) % arma::max(arma::abs(U) - lam1, arma::zeros(p, K));
        
    } else if (regul == "L2") {
        // Elementwise L2 Scaling (Ridge): u / (1 + 2*lam2)
        if (lam2 == 0.0) return; 
        U /= (1.0 + 2.0 * lam2);
        
    } else if (regul == "ElasticNet") {
        // Composition: Apply L2 scaling then L1 thresholding
        if (lam1 == 0.0 && lam2 == 0.0) return; // No penalty
        
        //  L2 Scaling (if lam2 > 0)
        if (lam2 > 0.0) {
            U /= (1.0 + 2.0 * lam2);
        }
        // L1 Soft-thresholding (if lam1 > 0)
        if (lam1 > 0.0) {
            U = arma::sign(U) % arma::max(arma::abs(U) - lam1, arma::zeros(p, K));
        }
        
    } else if (regul == "GroupLasso") {
        // Group Lasso Soft-thresholding (Block Soft-thresholding).
        if (lam1 == 0.0) return; 
        if (grp_id.length() != p) {
            Rcpp::stop("grp_id length must match number of rows (p) for GroupLasso.");
        }
        
        std::unordered_map<int, std::vector<arma::uword>> groups;
        for(int i = 0; i < p; ++i) {
            groups[grp_id[i]].push_back(i); 
        }
        
        // Iterate through each identified group
        for (auto const& [id, row_indices_vec] : groups) {
            if (row_indices_vec.empty()) continue; // Skip empty groups if any
            
            arma::uvec row_indices(row_indices_vec);
            
            double group_norm = arma::norm(U.rows(row_indices), "fro");
            
            // Apply block soft-thresholding
            double scaling_factor = std::max(0.0, 1.0 - lam1 / group_norm);
            U.rows(row_indices) *= scaling_factor;
            
        } 
        
    } else if (regul == "FusedLasso" || regul == "Graph" || regul == "Tree") {
        Rcpp::warning("Proximal operator for '%s' not implemented in this function. Skipping proximal step.", regul.c_str());
        
    } else {
        Rcpp::warning("Unsupported 'regul' type in proximalFlat2_nospams: '%s'. Skipping proximal step.", regul.c_str());
    }
    
    if (pos) {
        U.elem( arma::find(U < 0.0) ).zeros(); 
    }
    
}

void proximalGraph2(
        arma::mat& U,
        int& p,              // # of rows of U
        int& K,              // # of columns of U
        std::string& regul,
        arma::mat& grp,      // grp is a matrix for proximalGraph and Tree
        arma::mat& grpV,
        Rcpp::NumericVector& etaG,
        int num_threads,
        double lam1,
        double lam2 = 0.0,
        double lam3 = 0.0,
        bool intercept = false,
        bool resetflow = false,
        bool verbose = false,
        bool pos = false,
        bool clever = true,
        bool eval = true,
        int size_group = 1,
        bool transpose = false) {
    
    // U dimensions
    // int p = U.n_rows;
    // int K = U.n_cols;
    // int pK = int (p*K);
    
    // read in U and convert to spams::Matrix<double> alpha0
    // double* ptrU = new double[pK];
    // for (int r = 0; r < p; r++) {
    //     for (int c = 0; c < K; c++) {
    //         ptrU[c * p + r] = U(r, c);
    //     }
    // }
    //Matrix<double> alpha0(ptrU,p,K);
    Matrix<double> alpha0(const_cast<double*>(U.memptr()), p, K);
    
    
    // grp dimensions
    int gr = grp.n_rows;
    int gc = grp.n_cols;
    int grc = int (gr * gc);
    
    // read in grp and convert to spams::Matrix<bool> grp_dense
    // then to spams::SpMatrix<bool> groups
    bool* ptrG = new bool[grc];
    for (int r = 0; r < gr; r++) {
        for (int c = 0; c < gc; c++) {
            ptrG[c * gr + r] = (grp(r, c) != 0.0);
        }
    }
    Matrix<bool> grp_dense(ptrG, gr, gc);
    SpMatrix<bool> groups;
    grp_dense.toSparse(groups);
    
    
    
    // grpV dimensions
    int gvr = grpV.n_rows;
    int gvc = grpV.n_cols;
    int gvrc = int (gvr * gvc);
    
    // read in grpV and convert to spams::Matrix<bool> grpV_dense
    // then to spams::SpMatrix<bool> groups_var
    bool* ptrGV = new bool[gvrc];
    for (int r = 0; r < gvr; r++) {
        for (int c = 0; c < gvc; c++) {
            ptrGV[c * gvr + r] = (grpV(r, c) != 0.0);
        }
    }
    Matrix<bool> grpV_dense(ptrGV, gvr, gvc);
    SpMatrix<bool> groups_var;
    grpV_dense.toSparse(groups_var);
    
    
    
    // read in etaG and convert to spams::Vector<int> eta_g
    int n_etaG = etaG.length();
    double* ptrEG = new double[n_etaG];
    for (int i = 0; i < n_etaG; i++) {
        ptrEG[i] = etaG(i);
    }
    Vector<double> eta_g(ptrEG, n_etaG);
    
    
    
    // read in regul and convert to char
    // int l = regul.length();
    // char* name_regul = new char[l];
    // for (int i = 0; i < l+1; i++) {
    //     name_regul[i] = regul[i];
    // }
    std::vector<char> name_regul_buf(regul.begin(), regul.end()); // Copy string content
    name_regul_buf.push_back('\0'); // Add the null terminator
    char* name_regul_ptr = name_regul_buf.data();
    
    
    
    // Initialize alpha - proximal operator
    Matrix<double> alpha(p,K);
    alpha.setZeros();
    
    
    // call _proximalGraph
    _proximalGraph(&alpha0, &alpha,
                    &eta_g, &groups, &groups_var,
                    num_threads, lam1, lam2, lam3,
                    intercept, resetflow, name_regul_ptr,
                    verbose, pos, clever, eval,
                    size_group, transpose);
    
    
    
    // put updated alpha back into U
    // for (int r = 0; r < p; r++) {
    //     for (int c = 0; c < K; c++) {
    //         U(r, c) = alpha[c * p + r];
    //     }
    // }
    
    // free the dynamic memory
    // delete[] ptrU;
    delete[] ptrG;
    delete[] ptrGV;
    delete[] ptrEG;
    // delete[] name_regul;
}




void proximalTree2(
        arma::mat& U,
        int& p,              // # of rows of U
        int& K,              // # of columns of U
        std::string& regul,
        arma::mat& grp,      // grp is a matrix for proximalGraph and Tree
        Rcpp::NumericVector& etaG,
        Rcpp::IntegerVector& own_var,
        Rcpp::IntegerVector& N_own_var,
        int num_threads,
        double lam1,
        double lam2 = 0.0,
        double lam3 = 0.0,
        bool intercept = false,
        bool resetflow = false,
        bool verbose = false,
        bool pos = false,
        bool clever = true,
        bool eval = true,
        int size_group = 1,
        bool transpose = false) {
    
    
    
    // U dimensions
    // int p = U.n_rows;
    // int K = U.n_cols;
    // int pK = int (p*K);
    
    // read in U and convert to spams::Matrix<double> alpha0
    // double* ptrU = new double[pK];
    // for (int r = 0; r < p; r++) {
    //     for (int c = 0; c < K; c++) {
    //         ptrU[c * p + r] = U(r, c);
    //     }
    // }
    //Matrix<double> alpha0(ptrU,p,K);
    Matrix<double> alpha0(const_cast<double*>(U.memptr()), p, K);
    
    
    
    // grp dimensions
    int gr = grp.n_rows;
    int gc = grp.n_cols;
    int grc = int (gr * gc);
    
    // read in grp and convert to spams::Matrix<bool> grp_dense
    // then to spams::SpMatrix<bool> groups
    bool* ptrG = new bool[grc];
    for (int r = 0; r < gr; r++) {
        for (int c = 0; c < gc; c++) {
            ptrG[c * gr + r] = ( grp(r, c) != 0.0 );
        }
    }
    Matrix<bool> grp_dense(ptrG, gr, gc);
    SpMatrix<bool> groups;
    grp_dense.toSparse(groups);
    
    
    
    // read in etaG and convert to spams::Vector<double> eta_g
    int n_etaG = etaG.length();
    double* ptrEG = new double[n_etaG];
    for (int i = 0; i < n_etaG; i++) {
        ptrEG[i] = etaG(i);
    }
    Vector<double> eta_g(ptrEG, n_etaG);
    
    
    
    // read in own_var and convert to spams::Vector<int> own_variables
    int n_ov = own_var.length();
    int* ptrOV = new int[n_ov];
    for (int i = 0; i < n_ov; i++) {
        ptrOV[i] = own_var(i);
    }
    Vector<int> own_variables(ptrOV, n_ov);
    
    
    
    // read in N_own_var and convert to spams::Vector<int> N_own_variables
    int n_Nov = N_own_var.length();
    int* ptrNOV = new int[n_Nov];
    for (int i = 0; i < n_Nov; i++) {
        ptrNOV[i] = N_own_var(i);
    }
    Vector<int> N_own_variables(ptrNOV, n_Nov);
    
    
    
    // read in regul and convert to char
    // int l = regul.length();
    // char* name_regul = new char[l];
    // for (int i = 0; i < l+1; i++) {
    //     name_regul[i] = regul[i];
    // }
    std::vector<char> name_regul_buf(regul.begin(), regul.end()); // Copy string content
    name_regul_buf.push_back('\0'); // Add the null terminator
    char* name_regul_ptr = name_regul_buf.data();
    
    
    
    // Initialize alpha - proximal operator
    Matrix<double> alpha(p,K);
    alpha.setZeros();
    
    
    
    // call _proximalTree
    _proximalTree(&alpha0, &alpha,
                   &eta_g, &groups,
                   &own_variables, &N_own_variables,
                   num_threads, lam1, lam2, lam3,
                   intercept, resetflow, name_regul_ptr,
                   verbose, pos, clever, eval,
                   size_group, transpose);
    
    
    
    // put updated alpha back into U
    for (int r = 0; r < p; r++) {
        for (int c = 0; c < K; c++) {
            U(r, c) = alpha[c * p + r];
        }
    }
    
    // free the dynamic memory
    // delete[] ptrU;
    delete[] ptrG;
    delete[] ptrEG;
    delete[] ptrOV;
    delete[] ptrNOV;
    // delete[] name_regul;
}






// [[Rcpp::export]]
arma::vec grad_ls_loss2(
        arma::rowvec& x,
        double& y,
        arma::vec& param,
        int& p) {
    arma::vec grad(p);
    grad =  arma::vectorise(x) * (arma::dot(x, param) - y);
    return grad;
}

// [[Rcpp::export]]
arma::vec grad_logistic_loss2(
        arma::rowvec& x,
        double& y,
        arma::vec& param,
        int& p) {
    arma::vec grad(p);
    double sig;
    sig = 1.0/( 1.0 + exp(-arma::dot(x, param)) );
    grad =  arma::vectorise(x) * (sig - y);
    return grad;
}


// Modified 
void grad_multinom_loss2(
        const arma::rowvec& x,
        int y, 
        int K,
        double offset, 
        const arma::mat& param,
        int p,
        arma::mat& grad_out) { 
    
    
    arma::vec pi(K); 
    
    // Calculate softmax probabilities pi
    for (int k = 0; k < K; k++) {
        pi(k) = exp(arma::dot(x, param.col(k)) + offset);
    }
    
    double sum_pi = arma::sum(pi);
    // Prevent division by zero 
    if (sum_pi > 1e-10) { 
        pi /= (sum_pi + 1.0);
    } else {
        pi.fill(1.0 / K);
    }
    
    arma::vec x_col = x.t(); 
    
    // Calculate gradient columns
    if (y == 0) { 
        for (int k = 0; k < K; k++) {
            grad_out.col(k) = pi(k) * x_col; 
        }
    } else { 
        int y_idx = y - 1; 
        for (int k = 0; k < K; k++) {
            if (k == y_idx) {
                grad_out.col(k) = (pi(k) - 1.0) * x_col; 
            } else {
                grad_out.col(k) = pi(k) * x_col;     
            }
        }
    }
}


// [[Rcpp::export]]
Rcpp::List MultinomLogistic2(
        arma::mat X,
        arma::vec Y,
        arma::vec offset,
        int K,
        int reg_p,              // number of regularized variables
        int penalty,            // 1 - proximalFlat; 2 - proximalGraph; 3 - proximalTree
        std::string regul,
        bool transpose,
        Rcpp::IntegerVector grp_id,
        Rcpp::NumericVector etaG,
        arma::mat grp,
        arma::mat grpV,
        Rcpp::IntegerVector own_var,
        Rcpp::IntegerVector N_own_var,
        double lam1,
        double lam2,
        double lam3,
        double learning_rate,
        double tolerance,
        int niter_inner,
        int maxit,
        int ncores) {
    
    int p = X.n_cols;
    int n = X.n_rows;
    
    // indexing, stochastic sampling, etc..
    int index;    // index of the stochastic sample
    arma::rowvec x_sample(p);
    int y_sample;
    double o_sample; // o for offset
    
    // gradient and related
    arma::mat grad(p, K);
    arma::mat temp1(p, K);
    arma::mat temp2(p, K);
    
    // parameter and related
    arma::mat param(p, K);
    param.zeros();
    arma::mat param_old(p, K);
    
    arma::mat param_reg(reg_p, K);   // regularized coefficients
    arma::mat param_t(K, reg_p);      // regularized coefficients
    
    // convergence related
    double diff;
    int counter_outer = 0;
    
    // compute mu: mean gradient at param_old
    while (true) {
        param_old = param; 
        arma::mat single_grad(p, K);
        grad.zeros();
        
        for (int i = 0; i < n; i++) {
            const auto& x_sample = X.row(i);
            const int& y_sample = Y(i);
            const double& o_sample = offset(i);

            grad_multinom_loss2(x_sample,
                                       y_sample,
                                       K,
                                       o_sample,
                                       param_old,
                                       p,
                                       single_grad);

            grad += single_grad;
        }

        // Divide the final SUM by n ONCE
        if (n > 0) {
            grad /= n;
        }
        
        
        // inner loop
        for (int i = 0; i < niter_inner; ++i) {
            
            index = arma::randi(arma::distr_param(0, n - 1));
            
            // x_sample = X.row(index);
            const auto& x_sample = X.row(index);
            y_sample = Y(index);
            o_sample = offset(index);
            
            // temp1 = grad_multinom_loss2(x_sample, y_sample, K, o_sample, param, p);
            // temp2 = grad_multinom_loss2(x_sample, y_sample, K, o_sample, param_old, p);
            grad_multinom_loss2(x_sample, y_sample, K, o_sample, param, p, temp1);
            grad_multinom_loss2(x_sample, y_sample, K, o_sample, param_old, p, temp2);
            
            param = param - learning_rate * (temp1 - temp2 + grad);
            
            // extract only variables involved in the penalization
            param_reg = param.head_rows(reg_p);
            
            // call proximal function
            if (transpose) {
                param_t = param_reg.t();
                
                if (penalty == 1) {
                    // proximalFlat2(param_t, K, reg_p, regul, grp_id, ncores, lam1 * learning_rate, lam2 * learning_rate, lam3 * learning_rate);
                    arma::mat param_reg_copy = param_reg; // Explicit copy
                    proximalFlat2_nospams(param_reg_copy, reg_p, K, regul, grp_id,
                                          lam1 * learning_rate,
                                          lam2 * learning_rate,
                                          lam3 * learning_rate);
                    param_reg = param_reg_copy; 
                }
                
                if (penalty == 2) {
                    proximalGraph2(param_t, K, reg_p, regul, grp, grpV, etaG, ncores, lam1 * learning_rate, lam2 * learning_rate);
                }
                
                if (penalty == 3) {
                    proximalTree2(param_t, K, reg_p, regul, grp, etaG, own_var, N_own_var, ncores, lam1 * learning_rate, lam2 * learning_rate);
                }
                
                param_reg = param_t.t();
            } else {
                if (penalty == 1) {
                    proximalFlat2(param_reg, reg_p, K, regul, grp_id, ncores, lam1 * learning_rate, lam2 * learning_rate, lam3 * learning_rate);
                }
                
                if (penalty == 2) {
                    proximalGraph2(param_reg, reg_p, K, regul, grp, grpV, etaG, ncores, lam1 * learning_rate, lam2 * learning_rate);
                }
                
                if (penalty == 3) {
                    proximalTree2(param_reg, reg_p, K, regul, grp, etaG, own_var, N_own_var, ncores, lam1 * learning_rate, lam2 * learning_rate);
                }
            }
            
            param.head_rows(reg_p) = param_reg;
        }
        
        counter_outer += 1;
        Rcpp::Rcout << "\n Iteration " << counter_outer <<"\n";
        
        diff = arma::norm(param - param_old, "fro");
        diff = diff/(p*K);
        Rcpp::Rcout << "Frobenius norm of coefficient update \n" << diff <<"\n";
        
        if (diff < tolerance || counter_outer>= maxit) {
            break;
        }
    }
    
    arma::sp_mat param_sp(param);
    
    Rcpp::List result = Rcpp::List::create(Rcpp::Named("Estimates")                     = param,
                                           Rcpp::Named("Sparse Estimates")              = param_sp);
    return result;
}





// [[Rcpp::export]]
Rcpp::List MultinomLogistic3(
        arma::mat X,
        arma::vec Y,
        arma::vec offset,
        int K,
        int reg_p,                // number of regularized variables
        int penalty,              // 1 - proximalFlat; 2 - proximalGraph; 3 - proximalTree
        std::string regul,
        bool transpose,
        Rcpp::IntegerVector grp_id,
        Rcpp::NumericVector etaG,
        arma::mat grp,
        arma::mat grpV,
        Rcpp::IntegerVector own_var,
        Rcpp::IntegerVector N_own_var,
        double lam1,
        double lam2,
        double lam3,
        double learning_rate,
        double tolerance,
        int niter_inner,
        int maxit,
        int ncores) {
    
    int p = X.n_cols;
    int n = X.n_rows;
    
    int index;
    int y_sample;
    double o_sample;
    
    arma::mat grad(p, K);
    arma::mat temp1(p, K);
    arma::mat temp2(p, K);
    
    arma::mat param(p, K);
    param.zeros(); 
    arma::mat param_old(p, K);
    
    arma::mat param_t;
    if (transpose) {
        param_t.set_size(K, reg_p);
    }
    
    double diff;
    int counter_outer = 0;
    
    std::vector<arma::mat> param_history;
    
    while (true) {
        param_old = param; 
        arma::mat single_grad(p, K);
        grad.zeros();               
        
        for (int i = 0; i < n; i++) {
            const auto& x_sample_view = X.row(i);     
            const int& y_sample_ref = Y(i);           
            const double& o_sample_ref = offset(i);   
            
            grad_multinom_loss2( 
                x_sample_view,
                y_sample_ref,
                K,
                o_sample_ref,
                param_old, 
                p,
                single_grad);
            
            grad += single_grad; 
        }
        if (n > 0) {
            grad /= n; 
        }
        
        
        for (int i = 0; i < niter_inner; ++i) {
            index = arma::randi(arma::distr_param(0, n - 1)); 
            const auto& x_sample_view = X.row(index);  
            y_sample = Y(index);  
            o_sample = offset(index);                  
        
            // Calculate stochastic 
            grad_multinom_loss2(x_sample_view, y_sample, K, o_sample, param, p, temp1);
            grad_multinom_loss2(x_sample_view, y_sample, K, o_sample, param_old, p, temp2);
            
            // parameter update
            param = param - learning_rate * (temp1 - temp2 + grad);
            
            auto param_reg = param.head_rows(reg_p);
            
            // call proximal function
            if (transpose) {
                param_t = param_reg.t();
                
                if (penalty == 1) {
                    proximalFlat2(param_t, K, reg_p, regul, grp_id, ncores, lam1 * learning_rate, lam2 * learning_rate, lam3 * learning_rate);
                }
                
                if (penalty == 2) {
                    proximalGraph2(param_t, K, reg_p, regul, grp, grpV, etaG, ncores, lam1 * learning_rate, lam2 * learning_rate);
                }
                
                if (penalty == 3) {
                    proximalTree2(param_t, K, reg_p, regul, grp, etaG, own_var, N_own_var, ncores, lam1 * learning_rate, lam2 * learning_rate);
                }
                
                param_reg = param_t.t();
            } else {
                arma::mat param_reg_mat = param_reg;
                if (penalty == 1) {
                    proximalFlat2(param_reg_mat, reg_p, K, regul, grp_id, ncores, lam1 * learning_rate, lam2 * learning_rate, lam3 * learning_rate);
                }
                
                if (penalty == 2) {
                    proximalGraph2(param_reg_mat, reg_p, K, regul, grp, grpV, etaG, ncores, lam1 * learning_rate, lam2 * learning_rate);
                }
                
                if (penalty == 3) {
                    proximalTree2(param_reg_mat, reg_p, K,  regul, grp, etaG, own_var, N_own_var, ncores, lam1 * learning_rate, lam2 * learning_rate);
                }
                param_reg = param_reg_mat;
            }
            
        } 
        
        counter_outer += 1; 
        Rcpp::Rcout << "\n Iteration " << counter_outer << "\n"; 
        
        
        diff = arma::norm(param - param_old, "fro");
        
        diff = diff / ( static_cast<double>(p * K) );
        Rcpp::Rcout << "Mean Frobenius norm of coefficient update \n" << diff << "\n";
        

        param_history.push_back(param);
        
        if (diff < tolerance || counter_outer >= maxit) {
            break; 
        }
    } 
    
    // Convert final parameters to sparse format
    arma::sp_mat param_sp(param);
    
    // Create the output list for R
    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("Estimates") = param,          
        Rcpp::Named("Sparse Estimates") = param_sp, 
        Rcpp::Named("CoefficientHistory") = param_history 
    );
    return result; 
}


// [[Rcpp::export]]
Rcpp::List MultinomLogisticExp(
        arma::mat X,
        arma::vec Y,
        arma::vec offset,
        int K,
        int reg_p,                // number of regularized variables
        int penalty,              // 1 - proximalFlat; 2 - proximalGraph; 3 - proximalTree
        std::string regul,
        bool transpose,
        Rcpp::IntegerVector grp_id,
        Rcpp::NumericVector etaG,
        arma::mat grp,
        arma::mat grpV,
        Rcpp::IntegerVector own_var,
        Rcpp::IntegerVector N_own_var,
        double lam1,
        double lam2,
        double lam3,
        double learning_rate,
        double tolerance,
        int niter_inner,
        int maxit,
        int ncores) {
    
    int p = X.n_cols;
    int n = X.n_rows;
    
    int index;
    int y_sample;
    double o_sample;
    
    arma::mat grad(p, K);
    arma::mat temp1(p, K);
    arma::mat temp2(p, K);
    
    arma::mat param(p, K);
    param.zeros(); 
    arma::mat param_old(p, K);
    
    arma::mat param_t;
    if (transpose) {
        param_t.set_size(K, reg_p);
    }
    
    double diff;
    int counter_outer = 0;
   
    std::vector<arma::mat> param_history;

    
    while (true) {
        param_old = param; 
        
        arma::mat single_grad(p, K); 
        grad.zeros();        
        
        for (int i = 0; i < n; i++) {
            const auto& x_sample_view = X.row(i);   
            const int& y_sample_ref = Y(i);         
            const double& o_sample_ref = offset(i); 
            
            grad_multinom_loss2( 
                x_sample_view,
                y_sample_ref,
                K,
                o_sample_ref,
                param_old, 
                p,
                single_grad); 
            
            grad += single_grad; 
        }
        if (n > 0) {
            grad /= n; 
        }
        
        
        // inner Loop 
        for (int i = 0; i < niter_inner; ++i) {
            index = arma::randi(arma::distr_param(0, n - 1)); // Sample index
            
            const auto& x_sample_view = X.row(index);  
            y_sample = Y(index);                       
            o_sample = offset(index);                  
            
            grad_multinom_loss2(x_sample_view, y_sample, K, o_sample, param, p, temp1);
            grad_multinom_loss2(x_sample_view, y_sample, K, o_sample, param_old, p, temp2);
            
            param = param - learning_rate * (temp1 - temp2 + grad);
            
            auto param_reg = param.head_rows(reg_p);
            
            if (transpose) {
                param_t = param_reg.t();
                
                if (penalty == 1) {
                    proximalFlat2(param_t, K, reg_p, regul, grp_id, ncores, lam1 * learning_rate, lam2 * learning_rate, lam3 * learning_rate);
                }
                
                if (penalty == 2) {
                    proximalGraph2(param_t, K, reg_p, regul, grp, grpV, etaG, ncores, lam1 * learning_rate, lam2 * learning_rate);
                }
                
                if (penalty == 3) {
                    proximalTree2(param_t, K, reg_p, regul, grp, etaG, own_var, N_own_var, ncores, lam1 * learning_rate, lam2 * learning_rate);
                }
                
                param_reg = param_t.t();
            } else {
                arma::mat param_reg_mat = param_reg;
                if (penalty == 1) {
                    proximalFlat2(param_reg_mat, reg_p, K, regul, grp_id, ncores, lam1 * learning_rate, lam2 * learning_rate, lam3 * learning_rate);
                }
                
                if (penalty == 2) {
                    proximalGraph2(param_reg_mat, reg_p, K, regul, grp, grpV, etaG, ncores, lam1 * learning_rate, lam2 * learning_rate);
                }
                
                if (penalty == 3) {
                    proximalTree2(param_reg_mat, reg_p, K,  regul, grp, etaG, own_var, N_own_var, ncores, lam1 * learning_rate, lam2 * learning_rate);
                }
            param_reg = param_reg_mat;
            }
            
        }
        
        counter_outer += 1; 
        Rcpp::Rcout << "\n Iteration " << counter_outer << "\n"; 
        
        diff = arma::norm(param - param_old, "fro")/(arma::norm(param_old, "fro")+ 1e-10);
        diff = diff / ( static_cast<double>(p * K) );
        Rcpp::Rcout << "Relative mean Frobenius norm of coefficient update \n" << diff << "\n";
        
      
        param_history.push_back(param); 
        
        if (diff < tolerance || counter_outer >= maxit) {
            break; 
        }
    } 
    
    arma::sp_mat param_sp(param);
    
    // Create the output list for R
    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("Estimates") = param,          
        Rcpp::Named("Sparse Estimates") = param_sp,
        Rcpp::Named("CoefficientHistory") = param_history 
    );
    return result;
}



// [[Rcpp::export]]
Rcpp::List MultinomLogisticAcc( 
        arma::mat X,
        arma::vec Y,
        arma::vec offset,
        int K,
        int reg_p,
        int penalty,
        std::string regul,
        bool transpose,
        Rcpp::IntegerVector grp_id,
        Rcpp::NumericVector etaG,
        arma::mat grp,
        arma::mat grpV,
        Rcpp::IntegerVector own_var,
        Rcpp::IntegerVector N_own_var,
        double lam1,
        double lam2,
        double lam3,
        double learning_rate,
        double momentum_gamma,    // Momentum parameter
        double tolerance,
        int niter_inner,
        int maxit,
        int ncores,
        bool pos = false) {       // Positivity constraint
    
    int p = X.n_cols;
    int n = X.n_rows;
    
    int index;
    int y_sample_int;
    double o_sample;
    
    arma::mat grad(p, K);
    arma::mat temp1(p, K);
    arma::mat temp2(p, K);
    
    arma::mat param(p, K, arma::fill::zeros);
    arma::mat param_old(p, K);
    arma::mat param_unprox(p, K); 
    
    arma::mat velocity(p, K, arma::fill::zeros); // Velocity for momentum
    
    arma::mat param_t; 
    if (transpose) {
        param_t.set_size(K, reg_p); 
    }
    
    double diff;
    int counter_outer = 0;
    std::vector<arma::mat> param_history;
    
    while (true) {
        param_old = param;
        arma::mat single_grad(p, K);
        grad.zeros();
        
        for (int i = 0; i < n; i++) {
            const auto& x_sample_view = X.row(i);
            const int& y_val = Y(i); // Assuming Y elements are suitable for int
            const double& o_val = offset(i);
            grad_multinom_loss2(x_sample_view, y_val, K, o_val, param_old, p, single_grad);
            grad += single_grad;
        }
        if (n > 0) {
            grad /= n;
        }
        
        for (int i = 0; i < niter_inner; ++i) {
            index = arma::randi(arma::distr_param(0, n - 1));
            const auto& x_sample_view = X.row(index);
            y_sample_int = static_cast<int>(Y(index)); 
            o_sample = offset(index);
            
            grad_multinom_loss2(x_sample_view, y_sample_int, K, o_sample, param, p, temp1);
            grad_multinom_loss2(x_sample_view, y_sample_int, K, o_sample, param_old, p, temp2);
            arma::mat svr_grad_current = temp1 - temp2 + grad;
            
            // Update velocity
            velocity = momentum_gamma * velocity - learning_rate * svr_grad_current;
            param_unprox = param + velocity;
            
            auto param_unprox_reg_view = param_unprox.head_rows(reg_p);
            
            if (transpose) {
                param_t = param_unprox_reg_view.t();
                if (penalty == 1) {
                    proximalFlat2(param_t, K, reg_p, regul, grp_id, ncores, lam1 * learning_rate, lam2 * learning_rate, lam3 * learning_rate, pos);
                } else if (penalty == 2) {
                    proximalGraph2(param_t, K, reg_p, regul, grp, grpV, etaG, ncores, lam1 * learning_rate, lam2 * learning_rate, pos);
                } else if (penalty == 3) {
                    proximalTree2(param_t, K, reg_p, regul, grp, etaG, own_var, N_own_var, ncores, lam1 * learning_rate, lam2 * learning_rate, pos);
                }
                param.head_rows(reg_p) = param_t.t(); 
            } else { 
                
                arma::mat param_unprox_reg_copy = param_unprox_reg_view; 
                if (penalty == 1) {

                    proximalFlat2(param_unprox_reg_copy, reg_p, K, regul, grp_id, ncores, lam1 * learning_rate, lam2 * learning_rate, lam3 * learning_rate, pos);
                } else if (penalty == 2) {
                    proximalGraph2(param_unprox_reg_copy, reg_p, K, regul, grp, grpV, etaG, ncores, lam1 * learning_rate, lam2 * learning_rate, pos);
                } else if (penalty == 3) {
                    proximalTree2(param_unprox_reg_copy, reg_p, K, regul, grp, etaG, own_var, N_own_var, ncores, lam1 * learning_rate, lam2 * learning_rate, pos);
                }
                param.head_rows(reg_p) = param_unprox_reg_copy; 
            }
            
            if (reg_p < p) {
                param.tail_rows(p - reg_p) = param_unprox.tail_rows(p - reg_p);
            }
        } 
        
        counter_outer += 1;
        param_history.push_back(param); 
        
        double norm_old = arma::norm(param_old, "fro");
        diff = arma::norm(param - param_old, "fro") / (norm_old + 1e-10);
        
        Rcpp::Rcout << "\n Iteration " << counter_outer << "\n";
        Rcpp::Rcout << "Relative Frobenius norm of coefficient update (approx) \n" << diff << "\n";
        
        if (diff < tolerance || counter_outer >= maxit) {
            break;
        }
    } 
    
    arma::sp_mat param_sp(param);
    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("Estimates") = param,
        Rcpp::Named("Sparse Estimates") = param_sp,
        Rcpp::Named("CoefficientHistory") = param_history
    );
    return result;
}


// [[Rcpp::export]]
Rcpp::List MultinomLogisticSARAH( 
        arma::mat X,
        arma::vec Y,
        arma::vec offset,
        int K,
        int reg_p,
        int penalty,
        std::string regul,
        bool transpose,
        Rcpp::IntegerVector grp_id,
        Rcpp::NumericVector etaG,
        arma::mat grp,
        arma::mat grpV,
        Rcpp::IntegerVector own_var,
        Rcpp::IntegerVector N_own_var,
        double lam1,
        double lam2,
        double lam3,
        double learning_rate,
        double tolerance,
        int niter_inner,
        int maxit,
        int ncores
        ) {
    
    int p = X.n_cols;
    int n = X.n_rows;
    
    int index;
    int y_sample_int;
    double o_sample;
    
    arma::mat current_grad_estimate(p, K);
    // Temporaries for stochastic gradients
    arma::mat stoch_grad_at_w_curr(p, K); 
    arma::mat stoch_grad_at_w_prev(p, K); 
    
    arma::mat param(p, K, arma::fill::zeros);   
    arma::mat param_prev_inner_step(p, K); 
    
    arma::mat param_t; 
    if (transpose) {
        param_t.set_size(K, reg_p);
    }
    
    double diff;
    int counter_outer = 0;
    std::vector<arma::mat> param_history;
    
    while (true) {
        arma::mat param_at_epoch_start = param; // This is w_0^(s) for current epoch s
        
        current_grad_estimate.zeros();
        arma::mat temp_for_full_grad_sum(p, K); 
        for (int i_full = 0; i_full < n; ++i_full) {
            const auto& x_sample_view_full = X.row(i_full);
            const int& y_sample_ref_full = Y(i_full);
            const double& o_sample_ref_full = offset(i_full);
            grad_multinom_loss2(x_sample_view_full, y_sample_ref_full, K, o_sample_ref_full,
                                param_at_epoch_start, p, temp_for_full_grad_sum);
            current_grad_estimate += temp_for_full_grad_sum;
        }
        if (n > 0) {
            current_grad_estimate /= n; 
        }
        
   
        param_prev_inner_step = param_at_epoch_start;
        
        for (int t = 0; t < niter_inner; ++t) {
           
            if (t > 0) {
                index = arma::randi(arma::distr_param(0, n - 1));
                const auto& x_sample_view_inner = X.row(index);
                y_sample_int = static_cast<int>(Y(index));
                o_sample = offset(index);

                grad_multinom_loss2(x_sample_view_inner, y_sample_int, K, o_sample,
                                    param, p, stoch_grad_at_w_curr);  
                grad_multinom_loss2(x_sample_view_inner, y_sample_int, K, o_sample,
                                    param_prev_inner_step, p, stoch_grad_at_w_prev); 
                
                current_grad_estimate = stoch_grad_at_w_curr - stoch_grad_at_w_prev + current_grad_estimate;
            }
            
            // Store current param (w_t^(s)) to be used as w_{t-1}^(s) 
            param_prev_inner_step = param;
            
            param = param - learning_rate * current_grad_estimate; 
            
            auto param_reg_view = param.head_rows(reg_p);
            if (transpose) {
                param_t = param_reg_view.t(); 
                if (penalty == 1) {
                    proximalFlat2(param_t, K, reg_p, regul, grp_id, ncores, lam1 * learning_rate, lam2 * learning_rate, lam3 * learning_rate);
                } else if (penalty == 2) {
                    proximalGraph2(param_t, K, reg_p, regul, grp, grpV, etaG, ncores, lam1 * learning_rate, lam2 * learning_rate);
                } else if (penalty == 3) {
                    proximalTree2(param_t, K, reg_p, regul, grp, etaG, own_var, N_own_var, ncores, lam1 * learning_rate, lam2 * learning_rate);
                }
                param.head_rows(reg_p) = param_t.t(); 
            } else { 
                arma::mat param_reg_copy = param_reg_view; // Explicit copy
                if (penalty == 1) {
                    proximalFlat2(param_reg_copy, reg_p, K, regul, grp_id, ncores, lam1 * learning_rate, lam2 * learning_rate, lam3 * learning_rate );
                } else if (penalty == 2) {
                    proximalGraph2(param_reg_copy, reg_p, K, regul, grp, grpV, etaG, ncores, lam1 * learning_rate, lam2 * learning_rate);
                } else if (penalty == 3) {
                    proximalTree2(param_reg_copy, reg_p, K, regul, grp, etaG, own_var, N_own_var, ncores, lam1 * learning_rate, lam2 * learning_rate);
                }
                param.head_rows(reg_p) = param_reg_copy; 
            }
           
        } 
        
        counter_outer += 1;
        param_history.push_back(param);
        
        double norm_param_at_epoch_start = arma::norm(param_at_epoch_start, "fro");
        diff = arma::norm(param - param_at_epoch_start, "fro") / (norm_param_at_epoch_start + 1e-10);
        
        diff = diff / ( static_cast<double>(p * K) + 1e-10 );
        
        Rcpp::Rcout << "\n Iteration " << counter_outer << "\n";
        Rcpp::Rcout << "Scaled relative Frobenius norm of coefficient update (vs epoch start) \n" << diff << "\n";
        
        if (diff < tolerance || counter_outer >= maxit) {
            break;
        }
    } 
    arma::sp_mat param_sp(param);
    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("Estimates") = param,
        Rcpp::Named("Sparse Estimates") = param_sp,
        Rcpp::Named("CoefficientHistory") = param_history
    );
    return result;
}