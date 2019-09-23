#include "../include/entropy.h"

void getEns(const arma::mat &m,double &mean,double& max){

    unsigned int N=m.n_cols;
    unsigned int i;

    double enmax=0;
    double enmean=0;
    double E;

    for(i=0;i<N;i++){
        E=entropy(m.col(i));
        if(enmax<E) enmax=E;
        enmean+=E/N;
    }

    mean=enmean;
    max=enmax;

}

double entropy(const arma::vec &vec, bool norm){

    unsigned int n=vec.n_rows;
    unsigned int i;
    double sum;
    double E=0;

    std::vector<double> v(n);

    if(norm){
        sum=pow(arma::norm(vec),2);
        for(i=0; i<n;++i) v[i]=pow(vec(i),2)/sum;
    } else {
        for(i=0; i<n;++i) v[i]=pow(vec(i),2);
    }

    for(i=0;i<n;++i){
        if(v[i]>0) E+=-v[i]*log2(v[i]);
    }

    return(E);

}
