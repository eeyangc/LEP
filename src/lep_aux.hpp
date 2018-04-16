////  Created by DaiMingwei on 16/11/8.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#ifndef lep_aux_hpp
#define lep_aux_hpp

#include <stdio.h>
#include <RcppArmadillo.h>
using namespace arma;
//#include <omp.h>

#define type_int 0
#define type_double 1

const int kPValueVer = 0;
const int kVbVer = 1;

class PerformanceObject{//object to record the performance
public:
    double FDR;
    double power;
    double auc;

};

class Vardist{
public:
    vec gamma;
    vec mu;
    vec sigma2beta;
    Vardist(uword P, double mu0, double gamma0){
        gamma = vec(P);
        gamma.fill(gamma0);
        mu = vec(P);
        mu.fill(0);
    }
};

class DataModel{ //model to store the original data
public:
    DataModel(mat* X, vec* y, Col<int>* labels, vec* beta, double h){
        this -> X = new Mat<double>(X -> memptr(), X->n_rows, X -> n_cols, false);
        this -> y = new Col<double>(y -> memptr(), y->size(), false);
        this -> labels = new Col<int>(labels -> memptr(), labels->size(), false);
        this -> beta = new Col<double> (beta -> memptr(), beta -> size(), false);
        this -> h = h;

    };
    ~DataModel(){
        delete this -> X;
        delete this -> y;
        delete this -> labels;
        delete this -> beta;
    }
    mat* X;
    vec* y;
    Col<int>* labels;
    vec* beta;
    double h; // heritability

};


class LEPfit{
    //class for the model generated

public:
    LEPfit( uword N, uword P,  uword K, uword iter, double L,  double sigam2e, double sigma2beta, double Pi, vec gammas, vec mu, vec S, vec* pParam, vec Xr, double cov, vec u_, vec v_){
        this -> N = N;
        this -> P = P;
        this -> K = K;
        this -> L = L;
        this -> iter = iter;

        this -> sigma2e = sigam2e;
        this -> sigma2beta = sigma2beta;
        this -> Pi = Pi;

        this -> gammas = gammas;
        this -> mu = mu;
        this -> S = S;
        this -> pParam = pParam;
        this -> Xr = Xr;
        this -> cov = cov;
        this -> u = u_;
        this -> v = v_;
        // this -> gamma_K = mat(gamma_K, P, K);

    }

    LEPfit(vec gammas, vec mu, double cov){
       this -> gammas = gammas;
       this -> mu = mu;
       this -> cov = cov;
    }
    ~LEPfit( ){

    }

    vec predict( mat* X );
    void cal_powerfdr(DataModel* model, double threshold,PerformanceObject* po);
    double cal_auc(DataModel* model);

    double sigma2e;
    double sigma2beta;
    double Pi; //number of variables
    double L; //lower bound
    uword N;
    uword P;
    uword K;
    uword iter;
    vec gammas;
    vec mu;
    vec S;
    vec* pParam;
    vec Xr;
    double cov;
    vec u;
    vec v;
    // mat gamma_K;
};

class Options{
public:
    Options(){
        this -> max_iter = 300;
        this -> display_gap = 60;
        this -> n_fold = 5;
    }

    Options(int max_iter, int display_gap){
        this -> max_iter = max_iter;
        this -> display_gap = display_gap;
    }

    Options(int max_iter, int display_gap, int n_fold){
        this -> max_iter = max_iter;
        this -> display_gap = display_gap;
        this -> n_fold = n_fold;
    }
    int max_iter;
    int display_gap;
    int n_fold;

};


vec fdr2FDR(vec fdr);

double lb_linear(vec ytilde, vec diagXTX, vec y, double sigma2e, Vardist vardist);
double lb_gamma(vec gamma, double log_pi);
double lb_klbeta(Vardist vardist, double sigma2beta);


//double dotX (double* x, double* y, int n);






/*update the parameters of the gamma distributions*/
// void update_betaparam(uword P, uword K, double gamma_sum, double * lpsummary, double* lpgamma, vec* lpparams);
// void update_betaparam(uword P, uword K, double gamma_sum, double * lpsummary, double* lpgamma, vec* lpparams, double* q11, double* q00);

// double lb_pvalue(uword P, uword K, double * lpsummary, double* lpgamma, vec* lpparams);
double lb_pvalue(uword P, uword K, double * lpsummary, double* lpgamma, vec* lpparams,double* q11, double* q00, double* u, double* v);

// void update_param(uword N, uword P, Vardist vardist, double sumyytilde, vec diagXTX, double& sigma2e, double& sigma2beta, double& pi_p);
void update_param(uword N, uword P, Vardist vardist, double sumyytilde, vec diagXTX, double& sigma2e, double& sigma2beta, double& pi_p, double* q11, double* q00, double* u, double* v, int K,double* lpgamma, double gamma_sum);

double clog(double v);
//remove the effects of covariates Z
template<class T>
void remove_cov(T* lpfX, int P, vec& y, mat* SZX, mat* SZy);

double calauc(vec label, vec pred);
//mat MatXfloat(mat& Zt, float* lpfX, int P);

void get_matrix_pointer(void* lpfX, int N, int P, int type, bool efficient,Mat<double> & Mat_f, void* & X_mat, bool & befloat);

Mat<double> cal_means(void* X_mat, bool befloat, int N);

void centering(void* X_mat, vec& y, double& mean_y, mat & SZX, bool befloat, int N, int P);

vec cal_diagXTX(const vec& y, void* X_mat, bool befloat, const mat& SZX, int N, mat& xty);


/*update the parameters of the gamma distributions*/
void update_betaparam(uword P, uword K, double gamma_sum, double * lpsummary, double* lpgamma, vec* lpparams, double* q11, double* q00);


void update_q11_q00(uword P, uword K, double * lpsummary,vec* lpparams, double* u, double* v, double* q11, double* q00);


#endif /* lep_aux_hpp */
