#ifndef Rcpp_aux_hpp
#define Rcpp_aux_hpp
#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "lep_aux.hpp"
#include "readPlink.hpp"
using namespace arma;
using namespace Rcpp;

int snps_overlap(vector<SNP>& x_snps, CharacterVector& ss_snps, Col<uword>& xindex, Col<uword>& sindex);


Options* read_opts(SEXP opts){
  Options* lp_opt = NULL;
  if(!Rf_isNull(opts)){
    Rcpp::List opt(opts);
    int n_fold = 5;
    int max_iter = 300;
    int display_gap= 60;
    if(opt.containsElementNamed("n_fold")){
      n_fold = opt["n_fold"];
    }
    if(opt.containsElementNamed("max_iter")){
      max_iter = opt["max_iter"];
    }

    if(opt.containsElementNamed("dis_gap")){
      display_gap = opt["dis_gap"];
    }
    lp_opt = new Options(max_iter, display_gap, n_fold);
  }
  return lp_opt;

}

template<class T>
mat combine_summary_X(Mat<T>& geno, vector<SNP>& X_snps, SEXP SS, double lbPval)
{
  int P = geno.n_cols;
  CharacterVector columnnames;
  CharacterVector snps_from_ss;
  if(!Rf_isNull(SS) && !Rf_isMatrix(SS)){
    Function asMatrix("as.matrix");
    SS = asMatrix(SS);
  }

  if(!Rf_isNull(SS) ){
    if(!Rf_isNull(rownames(SS))){
      snps_from_ss = rownames(SS);
    }
    if(!Rf_isNull(colnames(SS))){
      columnnames = colnames(SS);
    }
  }
  if(!Rf_isNull(SS) && X_snps.size() > 0){
    NumericMatrix SS_Matrix(SS);
    mat lp_summaries(SS_Matrix.begin(), SS_Matrix.nrow(), SS_Matrix.ncol(), false);
    Col<uword> xindex;
    Col<uword> sindex;
    if(snps_from_ss.size() > 0 ){
      snps_overlap(X_snps, snps_from_ss, xindex, sindex);
      lp_summaries = lp_summaries.rows(sindex);
      if(xindex.size() > 0 && xindex.size() < geno.n_cols){
        geno = geno.cols(xindex);
      }else if(xindex.size() == 0){
        // cout << "There are no intersection between genotype data and summary-level data!" << endl;
        // cout <<"LEP is running with only genotype data!" << endl;
      }
    }

    if(lp_summaries.n_cols > 0 && lp_summaries.n_rows > 0){
      uvec q = find(lp_summaries < lbPval);
      lp_summaries.elem(q).fill(lbPval);
    }
    if(lp_summaries.n_rows > lp_summaries.n_cols){
      lp_summaries = lp_summaries.t();
    }
    return lp_summaries;
  }else if(!Rf_isNull(SS) && (X_snps.size() == 0 || snps_from_ss.length() == 0)){
    NumericMatrix SS_Matrix(SS);
    // cout << "The SNPs' names for genotype data or Summary Statistics are "
    //      <<"not properly provided! " << endl;
    if(P == SS_Matrix.nrow() || P == SS_Matrix.ncol()){
     //  cout << "The numbers of SNPs are equal for genotype data and Summary Statistics! We are going to using the p-values!" << endl;
    }else{
       mat lp_summaries;
       lp_summaries.reset();
       return lp_summaries;
    }

    mat lp_summaries(SS_Matrix.begin(), SS_Matrix.nrow(), SS_Matrix.ncol(), false);
    if(lp_summaries.n_cols > 0 && lp_summaries.n_rows > 0){
      uvec q = find(lp_summaries < lbPval);
      lp_summaries.elem(q).fill(lbPval);
    }
    if(lp_summaries.n_rows > lp_summaries.n_cols){
      lp_summaries = lp_summaries.t();
    }
    return lp_summaries;
  }else{
    mat lp_summaries;
    lp_summaries.reset();
    return lp_summaries;
  }

}

void convert2mat(Mat<double>*& obj, SEXP input){
  if(!Rf_isNull(input)){

    if(!Rf_isMatrix(input)){
      Function asMatrix("as.matrix");
      input = asMatrix(input);
    }

    Rcpp::NumericMatrix T(input);
    obj = new Mat<double>(T.begin(), T.nrow(), T.ncol(), false);
  }
}

void get_N_P_type(SEXP Xs, int& N, int& P, int& type){
  if(!Rf_isNull(Xs) && !Rf_isMatrix(Xs)){
    Function asMatrix("as.matrix");
    Xs = asMatrix(Xs);
  }
  if(!Rf_isInteger(Xs)){
  //  cout <<"from double" << endl;
    NumericMatrix nm(Xs);
    N = nm.nrow();
    P = nm.ncol();
    type = type_double;
  }else{
 //   cout <<"from integer" << endl;
    IntegerMatrix nm(Xs);
    type = type_int;
    N = nm.nrow();
    P = nm.ncol();
  }
}

RcppExport SEXP wrap_fit(LEPfit* fit){
  Rcpp::List ret;
  ret["sigma2beta"] = fit -> sigma2beta;
  ret["sigma2e"] = fit -> sigma2e;
  ret["gammas"] = fit -> gammas;
  ret["mu"] = fit -> mu;
  ret["S"] = fit -> S;
  ret["pi"] = fit -> Pi;
  ret["M"] = fit ->P;
  ret["cov"] = fit -> cov;
  ret["L"] = fit -> L;
  ret["iter"] = fit -> iter;
  ret["u"] = fit -> u;
  ret["v"] = fit -> v;
  ret["fdr"] = 1 - fit -> gammas;
  if(fit -> pParam != NULL){
    ret["param_beta"] = (*fit -> pParam);
  }
  return ret;
}


// RcppExport SEXP wrap_fit(LEPfit* fit){
//   Rcpp::List ret;
//   ret["sigma2beta"] = fit -> sigma2beta;
//   ret["sigma2e"] = fit -> sigma2e;
//   ret["gammas"] = fit -> gammas;
//   ret["mu"] = fit -> mu;
//   ret["S"] = fit -> S;
//   ret["pi"] = fit -> Pi;
//   ret["P"] = fit -> P;
//   ret["fdr"] = 1 - fit -> gammas;
//   ret["cov"] = fit -> cov;
//   ret["L"] = fit -> L;
//   ret["iter"] = fit -> iter;
//   ret["u"] = fit -> u;
//   ret["v"] = fit -> v;
//   ret["param_beta"] = fit -> pParam;
//   return ret;
// }

int snps_overlap(vector<SNP>& x_snps, CharacterVector& ss_snps, Col<uword>& xindex, Col<uword>& sindex){
 //  vector<SNP> X_SNPS;
   vector<SNP> SS_SNPS;
   // for(int i = 0; i < x_snps.length(); i++){
   //      SNP snp(std::string(x_snps[i]), i, 1);
   //      X_SNPS.push_back(snp);
   // }
   for(int i = 0; i < ss_snps.length(); i++){
     SNP snp(std::string(ss_snps[i]), i, 1);
     SS_SNPS.push_back(snp);
   }
   return snps_overlap(x_snps, SS_SNPS, xindex, sindex);
}


// int snps_overlap(vector<SNP>& X_SNPS, CharacterVector& ss_snps, Col<uword>& xindex, Col<uword>& sindex){
//  // vector<SNP> X_SNPS;
//   vector<SNP> SS_SNPS;
//   // for(int i = 0; i < x_snps.length(); i++){
//   //   SNP snp(std::string(x_snps[i]), i, 1);
//   //   X_SNPS.push_back(snp);
//   // }
//   for(int i = 0; i < ss_snps.length(); i++){
//     SNP snp(std::string(ss_snps[i]), i, 1);
//     SS_SNPS.push_back(snp);
//   }
//   return snps_overlap(X_SNPS, SS_SNPS, xindex, sindex);
// }


void preprocess_summary(SEXP Xs, SEXP  SS, double lbPval, int& type, int& N, int& P, mat& summary, void* & ptr){
  if(!Rf_isNull(Xs) && !Rf_isMatrix(Xs)){
    Function asMatrix("as.matrix");
    Xs = asMatrix(Xs);
  }
  CharacterVector geno_snps = colnames(Xs);
  vector<SNP> xsnps;
  for(int i = 0; i < geno_snps.length(); i++){
    SNP snp(std::string(geno_snps[i]), i, 1);
    xsnps.push_back(snp);
  }
  get_N_P_type(Xs, N, P, type);
  if(type == type_int)
  {
    int* p = INTEGER(Xs);
    Mat<int> mat_X(p, N, P, false);
    summary = combine_summary_X(mat_X, xsnps, SS, lbPval);
    bool create_new_matrix = false;
    if(mat_X.n_cols < P){
      create_new_matrix = true;
      P = mat_X.n_cols;
    }
    ptr = new Mat<int>(mat_X.memptr(), mat_X.n_rows, P , create_new_matrix);
  }else{
    double* p = REAL(Xs);
    Mat<double> mat_X(p, N, P,false);
    summary = combine_summary_X(mat_X, xsnps, SS, lbPval);
    bool create_new_matrix = false;
    if(mat_X.n_cols < P){
        create_new_matrix = true;
        P = mat_X.n_cols;
    }
    ptr = new Mat<double>(mat_X.memptr(), mat_X.n_rows, P , create_new_matrix);
  }
}


#endif /* lep_hpp */
