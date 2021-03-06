// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// LEP
RcppExport SEXP LEP(SEXP Xs, arma::vec& y, SEXP SS, SEXP opts, std::string logfile, double lbPval, bool verbose);
RcppExport SEXP _LEP_LEP(SEXP XsSEXP, SEXP ySEXP, SEXP SSSEXP, SEXP optsSEXP, SEXP logfileSEXP, SEXP lbPvalSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< SEXP >::type SS(SSSEXP);
    Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
    Rcpp::traits::input_parameter< std::string >::type logfile(logfileSEXP);
    Rcpp::traits::input_parameter< double >::type lbPval(lbPvalSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(LEP(Xs, y, SS, opts, logfile, lbPval, verbose));
    return rcpp_result_gen;
END_RCPP
}
// LEPCV
RcppExport SEXP LEPCV(SEXP Xs, arma::vec& y, SEXP SS, SEXP opts, std::string logfile, double lbPval, Rcpp::String measure);
RcppExport SEXP _LEP_LEPCV(SEXP XsSEXP, SEXP ySEXP, SEXP SSSEXP, SEXP optsSEXP, SEXP logfileSEXP, SEXP lbPvalSEXP, SEXP measureSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< SEXP >::type SS(SSSEXP);
    Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
    Rcpp::traits::input_parameter< std::string >::type logfile(logfileSEXP);
    Rcpp::traits::input_parameter< double >::type lbPval(lbPvalSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type measure(measureSEXP);
    rcpp_result_gen = Rcpp::wrap(LEPCV(Xs, y, SS, opts, logfile, lbPval, measure));
    return rcpp_result_gen;
END_RCPP
}
// LEP_Plink
RcppExport SEXP LEP_Plink(Rcpp::String genoplinkfile, SEXP SS, SEXP opts, std::string logfile, double lbPval, bool verbose);
RcppExport SEXP _LEP_LEP_Plink(SEXP genoplinkfileSEXP, SEXP SSSEXP, SEXP optsSEXP, SEXP logfileSEXP, SEXP lbPvalSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type genoplinkfile(genoplinkfileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type SS(SSSEXP);
    Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
    Rcpp::traits::input_parameter< std::string >::type logfile(logfileSEXP);
    Rcpp::traits::input_parameter< double >::type lbPval(lbPvalSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(LEP_Plink(genoplinkfile, SS, opts, logfile, lbPval, verbose));
    return rcpp_result_gen;
END_RCPP
}
// LEPCV_Plink
RcppExport SEXP LEPCV_Plink(Rcpp::String genoplinkfile, SEXP SS, SEXP opts, std::string logfile, double lbPval, Rcpp::String measure);
RcppExport SEXP _LEP_LEPCV_Plink(SEXP genoplinkfileSEXP, SEXP SSSEXP, SEXP optsSEXP, SEXP logfileSEXP, SEXP lbPvalSEXP, SEXP measureSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type genoplinkfile(genoplinkfileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type SS(SSSEXP);
    Rcpp::traits::input_parameter< SEXP >::type opts(optsSEXP);
    Rcpp::traits::input_parameter< std::string >::type logfile(logfileSEXP);
    Rcpp::traits::input_parameter< double >::type lbPval(lbPvalSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type measure(measureSEXP);
    rcpp_result_gen = Rcpp::wrap(LEPCV_Plink(genoplinkfile, SS, opts, logfile, lbPval, measure));
    return rcpp_result_gen;
END_RCPP
}
// LEP_Predict
RcppExport SEXP LEP_Predict(SEXP fit_, arma::mat& X);
RcppExport SEXP _LEP_LEP_Predict(SEXP fit_SEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type fit_(fit_SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(LEP_Predict(fit_, X));
    return rcpp_result_gen;
END_RCPP
}
// read_plink
RcppExport SEXP read_plink(Rcpp::String genoplinkfile);
RcppExport SEXP _LEP_read_plink(SEXP genoplinkfileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type genoplinkfile(genoplinkfileSEXP);
    rcpp_result_gen = Rcpp::wrap(read_plink(genoplinkfile));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LEP_LEP", (DL_FUNC) &_LEP_LEP, 7},
    {"_LEP_LEPCV", (DL_FUNC) &_LEP_LEPCV, 7},
    {"_LEP_LEP_Plink", (DL_FUNC) &_LEP_LEP_Plink, 6},
    {"_LEP_LEPCV_Plink", (DL_FUNC) &_LEP_LEPCV_Plink, 6},
    {"_LEP_LEP_Predict", (DL_FUNC) &_LEP_LEP_Predict, 2},
    {"_LEP_read_plink", (DL_FUNC) &_LEP_read_plink, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_LEP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
