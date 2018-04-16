//  Created by DaiMingwei on 16/11/8.
//  Copyright © 2016年 daviddai,eeyang. All rights reserved.
//

#include "lep.hpp"



template <typename T>
double dotX (T* x, double* y, int n) {
  double sum = 0;
  //  #pragma omp parallel for reduction(+:sum)
  for (int i = 0; i < n; i++)
    sum += x[i] * y[i];
  return sum;
}



template<class T>
void addX (double* y, double a, T* x, int n) {
  //   #pragma omp parallel for num_threads(4)
  for (int i = 0; i < n; i++)
    y[i] += a * x[i];
}


template<typename T>
void lep_update(T* x_j, double* gamma, double* mu, double d, double s,  double logPi, double sigma2beta, double sigma2e, int N, double xy, double* ytilde_pt,
                double* u, double* v, double* lpsummay, vec* lpparam, double xmean_j, bool befloat){//

  double r = (*gamma) * (*mu);
  double numerator = xy + d * r -  dotX(x_j, ytilde_pt, N);
  mat y_tilde(ytilde_pt, N, 1, false);

  if(befloat == false){
    double sum_y_tilde = as_scalar(sum(y_tilde));
    numerator += xmean_j * sum_y_tilde;
  }
  (*mu) = s / sigma2e * numerator;

  double SSR_logratio(0);

  SSR_logratio = (*mu) * (*mu) / (2 * s) + 0.5 * log( s / sigma2beta);

  double ai = 0; //additional information provided by the summary statistisc
  if (lpsummay != NULL && lpparam != NULL) {
    uword K = lpparam -> size();
    double* lppara = lpparam -> memptr();
    for (uword i = 0; i < K; i++) {
      double beta_f = lppara[i] * pow(lpsummay[i], lppara[i] - 1);
      ai += log((beta_f * u[i] + 1 - u[i])/(beta_f * (1 - v[i]) + v[i]));//log(lppara[i]) + (lppara[i] - 1) * log(lpsummay[i]);
    }
    SSR_logratio += ai;
  }


  (*gamma) = 1/(1+exp(-(logPi + SSR_logratio)));

  double rnew = (*gamma) * (*mu);

  addX (ytilde_pt, rnew - r, x_j, N);
  if(befloat == false){
    y_tilde -= (rnew - r) * xmean_j;
  }
}



LEPfit* lep(void* lpfX, vec y, int P, mat* lpsummaryinfo, Options* opt, int type, bool efficient, bool verbose)
{
  uword N = y.n_rows;
  void* X_mat;
  Mat<double> Mat_f;
  bool befloat;
  get_matrix_pointer(lpfX, N, P, type, efficient, Mat_f, X_mat, befloat);
  mat SZX; //store column means of X, return by function 'centering'
  double mean_y;// = as_scalar(mean(y)); //mean of y, return by function 'centering'
  centering(X_mat, y, mean_y, SZX, befloat, N, P);

  uword K = lpsummaryinfo != NULL ? (lpsummaryinfo -> n_rows) : 0;

  vec u(K);
  u.fill(0.1);
  vec v(K);
  v.fill(0.9);

  double* q11 = new double[P*K];
  double* q00 = new double[P*K];

  // double* gamma_K= new double[P*K];

  for(int i = 0; i < P*K; i++){
    q11[i] = 0.1;
    q00[i] = 0.9;
  }

  opt = opt != NULL ? opt : new Options();
  int max_iter = opt -> max_iter;
  int display_gap = opt -> display_gap;
  mat xty;
  vec diagXTX = cal_diagXTX(y, X_mat, befloat, SZX, N, xty);

  double pi_p = 0.01; //pi for prior proportion
  double mu0 = 0;
  double alpha0 = 0.5; //initial parameters of beta distribtuion

  Vardist vardist(P, mu0, pi_p);

  double sigma2e = var(y) / 2;
  double sigma2beta = sigma2e;

  vec beta = vardist.gamma % vardist.mu;
  vec ytilde = vec(N);//(*lpfX) * beta;
  ytilde.zeros();

  Col<double>* lpparams = NULL;
  if ( K > 0 ){
    lpparams = new Col<double>(K);
    lpparams -> fill(alpha0); //parameters for beta distribtuions
  }

  double L0 = -INFINITY;
  double L = 0;
  double* lpgamma = vardist.gamma.memptr();
  double* lpmu = vardist.mu.memptr();
  double* lpd = diagXTX.memptr();
  double* lpytilde = ytilde.memptr();
  double* lpxy = xty.memptr();
  double* lpsummary = lpsummaryinfo != NULL ? lpsummaryinfo -> memptr() : NULL;
  uword iter = 0;
  double* mean_x = SZX.memptr();


  for (iter = 0; iter < max_iter; iter ++) {
    clock_t t1 = clock();
    if (verbose) {
       if(iter == 0)  cout <<"Begin Iterations" << endl;
    }
    double logPi = log(pi_p / (1 - pi_p));
    double sigma2e_Sigma2beta = sigma2e / sigma2beta;
    vec xxsigma = diagXTX + sigma2e_Sigma2beta;
    vardist.sigma2beta = sigma2e / xxsigma;
    double* S = vardist.sigma2beta.memptr();
    double gamma_sum = 0;
    if(befloat){
      Mat<double> * mat_f = static_cast<Mat<double> *>(X_mat);
      double* lp_Xf = mat_f -> memptr();
      for (int j = 0; j < P; j++) {
        lep_update(lp_Xf + N*j, lpgamma + j, lpmu + j, lpd[j], S[j], logPi, sigma2beta, sigma2e, (int)N,  lpxy[j], lpytilde, u.memptr(), v.memptr(), lpsummary + K * j, lpparams, mean_x[j], befloat);
        gamma_sum += *(lpgamma + j);
      }
    }else{
      Mat<int> * mat_i = static_cast<Mat<int> *>(X_mat);
      int* lp_Xi = mat_i -> memptr();
      for (int j = 0; j < P; j++) {
        lep_update(lp_Xi + N*j, lpgamma + j, lpmu + j, lpd[j], S[j], logPi, sigma2beta, sigma2e, (int)N,  lpxy[j], lpytilde, u.memptr(), v.memptr(), lpsummary + K * j, lpparams, mean_x[j], befloat);
        gamma_sum += *(lpgamma + j);
      }
    }

    //update alpha_k for beta-distribution
    update_betaparam(P, K, gamma_sum, lpsummary, lpgamma, lpparams, q11, q00);

    // update sigma2beta, sigma2e, u, v
    update_param(N, P, vardist, sum(square(y-ytilde)), diagXTX, sigma2e, sigma2beta, pi_p, q11, q00, u.memptr(), v.memptr(), (int)K, lpgamma, gamma_sum);


    L = lb_linear(ytilde, diagXTX,y, sigma2e, vardist) + lb_gamma(vardist.gamma, logPi) +lb_klbeta(vardist, sigma2beta)
           + lb_pvalue(P, K, lpsummary, lpgamma, lpparams, q11, q00, u.memptr(), v.memptr());

    update_q11_q00(P, K, lpsummary, lpparams, u.memptr(), v.memptr(), q11, q00);


    if(verbose && iter % display_gap == 0){
      cout <<iter <<"th iteration L=" << L << ";sigma2e = " << sigma2e << "; sigma2beta = " << sigma2beta << " ;pi = " <<
        pi_p <<" time = " << (clock() - t1)*1.0/CLOCKS_PER_SEC << endl;
      cout <<"uv for " << K <<" GWAS" << endl;
      cout <<"u =" << u.t() << endl;
      cout <<"v =" << v.t() << endl;
    }

    if(L < L0){
      cout << "Lowerbound decreasing, Error at iteration "<<iter <<"th iteration, diff = " << L-L0 << endl;
      break;
    }else if(L - L0 < 1e-5)//if(fabs((L - L0)/L0) < 1e-8) //
    {
      if(verbose){
        cout<<"Converge at " << iter << "th iteration, L = "<< L << endl;
      }
      break;
    }
    L0 = L;

  }

  // for(int j = 0; j < P; j++){
  //   for(int k = 0; k < K; k++){
  //     gamma_K[P*k + j] = lpgamma[j]*q11[P*k + j] + (1 - lpgamma[j])*(1-q00[P*k + j]);
  //   }
  // }

  mat vv = vardist.gamma % vardist.mu;
  double cov = mean_y - as_scalar(SZX * conv_to<mat>::from(vardist.gamma % vardist.mu));
  LEPfit* fit = new LEPfit(N, P,  K, iter, L,  sigma2e, sigma2beta, pi_p, vardist.gamma, vardist.mu
                                 , vardist.sigma2beta, lpparams,  ytilde, cov,u, v);
  return fit;
}


double lepCV(void* lpfX, vec y, int P, mat* lpsummaryinfo, Options* opt, int type, bool efficient, std::string measure, bool verbose){

    uword N = y.n_rows;
    void* X_mat;
    Mat<double> Mat_f;
    bool befloat;
    get_matrix_pointer(lpfX, N, P, type, efficient, Mat_f, X_mat, befloat);
    opt = opt != NULL ? opt : new Options();
    uword nfold = opt -> n_fold;
    Col<uword> indices = cross_valind(N, nfold);
    vec ylabel = y;
    vec predY(N);
    if (verbose) {
      cout << "Start " << nfold <<" cross validations!" << endl;
    }
    for (uword i = 1; i <= nfold; i++) {
        if (verbose) cout <<".";
        Col<uword> train_idx = find(indices != i);
        Col<uword> test_idx = find(indices == i);
        if( befloat ){
          Mat<double> * mat_f = static_cast<Mat<double> *>(X_mat);
          Mat<double> trainM = mat_f -> rows(train_idx);
          vec ytrain = y(train_idx);
          Mat<double> testM = mat_f -> rows(test_idx);
          vec ytest = y(test_idx);
          double* trainX = trainM.memptr();
          LEPfit* f = lep(trainX, ytrain, P,  lpsummaryinfo, opt, type, efficient, verbose);
          vec predy = f -> predict(&testM);
          predY.elem(test_idx) = predy;
          delete f;
        }else{
          Mat<int> * mat_i = static_cast<Mat<int> *>(X_mat);
          Mat<int> trainM = mat_i -> rows(train_idx);
          vec ytrain = y(train_idx);
          Mat<double> testM = conv_to<mat>::from(mat_i -> rows(test_idx));
          vec ytest = y(test_idx);
          int* trainX = trainM.memptr();
          LEPfit* f = lep(trainX, ytrain, P,  lpsummaryinfo, opt, type, efficient, verbose);
          vec predy = f -> predict(&testM);
          predY.elem(test_idx) = predy;
          delete f;
        }
    }

    double predict = 0;
    if(measure.compare("auc") == 0 ){
       predict = calauc(conv_to<vec>::from(ylabel), conv_to<vec>::from(predY));
    }else if(measure.compare("mse") == 0 ){
       vec diff = y - predY;
       predict = mean(diff % diff);
    }else if(measure.compare("cor") == 0 ){
       predict = as_scalar(cor(y, predY));
    }
    return predict;
}

/**shuffle the index for cross validation*/
arma::Col<uword> cross_valind(arma::uword N, arma::uword nfold){
    arma::Col<uword> indices(N);
    arma_rng::set_seed_random();
    arma::Col<uword> vec_n = arma::shuffle(arma::linspace < arma::Col <uword> >(1, N, N));

    indices.fill(nfold);
    for(uword n = 1; n <= nfold-1; n++){
        arma::Col<uword> in = vec_n.rows((n-1)*N/nfold,n*N/nfold - 1);
        indices.elem(in - 1 ).fill(n);
    }
    return indices;
}
