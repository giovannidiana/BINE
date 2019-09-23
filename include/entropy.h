#ifndef ENTROPY_H
#define ENTROPY_H

#include<armadillo>

void getEns(const arma::mat &m, double &mean, double &max);
double entropy(const arma::vec&,bool norm=true);

#endif // ENTROPY_H

