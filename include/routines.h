#ifndef _ROUTINES_H
#define _ROUTINES_H

#include<armadillo>
#include<fstream>
#include<gsl/gsl_rng.h>

using namespace std;

void gen_sample(arma::mat &mat, arma::uvec &t, int P,int);
void gen_sample2(arma::mat &mat, arma::uvec &t, int P,int);
void gen_sample3(gsl_rng* r, arma::fmat &m, arma::uvec &t, int P,double eps, double pmu);
void gen_sample4(gsl_rng* r, arma::fmat &m, arma::uvec &t, int P,double,double, double pmu);
void gen_sample4(gsl_rng* r, arma::imat &m, arma::uvec &t, int P,double lambda0,double lambda1, double pmu);
void gen_sample5(gsl_rng* r, arma::imat &m, arma::uvec &t, int P);
void gen_sample6(gsl_rng* r, arma::imat &m, arma::uvec &t, int P,double lambda0,double lambda1, double pmu,bool lastIsFree=false);
void gen_sample6_multimem(gsl_rng* r, arma::imat &m, arma::umat &t, int P,double lambda0,double lambda1, double pmu,bool lasIsFree=false);
void write_sample(ofstream &file, arma::mat);
bool match_prediction(const arma::uvec &, int* pred, int N, int P,int* res,double* acc);
bool match_prediction(const arma::uvec &, int* pred, int N, int P,int* res,double &acc);
bool match_prediction2(const arma::uvec &, int* pred, int N, int P,int* res,double &acc);
void expand_binary(arma::fmat &s,int size);
void filter_binary(const arma::fmat &in,arma::fmat &out, arma::uvec &selec, int threshold,int threshold2=0);
void filter_binary2(const arma::imat &in,arma::imat &out, arma::uvec &selec, int threshold,int threshold2=0);
void filter_binary(const arma::imat &in,arma::imat &out, arma::uvec &selec, int threshold);
void SVDrot(arma::mat &input, arma::uvec &cell_labels, unsigned int nguess);

// GIBBS
arma::vec dirichlet(gsl_rng *r, const arma::vec &alpha);

// Useful functions
bool fileExists(const char* file);
double gsl_ran_beta_pdflog(double x, double alpha, double beta);


#endif
