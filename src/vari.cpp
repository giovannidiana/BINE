#include "../include/vari.h"

Vari::Vari(const arma::mat &input) : factors(input)
{

}

void Vari::varimax(double eps=1e-5, bool normalize=true)
{
    int i,j;

    int nc = factors.n_cols;
    int nr = factors.n_rows;

    arma::vec sc(nr);
    arma::mat x(nr,nc);
    arma::mat z;
    arma::mat B;
    arma::rowvec nr1 = arma::ones<arma::rowvec>(nr);
    arma::mat U,V;
    arma::vec s;

    double dpast,d;

    if (nc < 2)
        rot=factors;

    if (normalize) {
        for(i=0;i<nr;++i){
            sc(i)=0;
            for(j=0;j<nc;++j){
                sc(i)+=pow(factors(i,j),2);
            }
            for(j=0;j<nc;++j) x(i,j)=factors(i,j)/sqrt(sc(i));
        }
    }

    rot = arma::mat(nc,nc,arma::fill::eye);

    d = 0;
    for (i=0;i<1000;++i) {
        z = x * rot;
        B = x.t()*(pow(z,3) - z*arma::diagmat(nr1*pow(z,2))/nr);

        arma::svd(U,s,V,B);
        rot = U * V.t();
        dpast = d;
        d = sum(s);
        if (d < dpast * (1 + eps))
            break;
    }

    z = x * rot;
    if (normalize)
        for(i=0;i<nr;++i) z.row(i) = z.row(i)/sc(i);

    loadings = z;

}


