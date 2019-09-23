#ifndef VARI_H
#define VARI_H

#include<armadillo>

class Vari
{
public:
    Vari(const arma::mat &);
    const arma::mat factors;
    void varimax(double,bool);
    arma::mat rot;
    arma::mat loadings;
};

#endif // VARI_H
