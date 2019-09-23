#ifndef INDEXED_ARRAY_H
#define INDEXED_ARRAY_H
#include<armadillo>

class indexed_array2x2{
public:
    indexed_array2x2(int ind);
    double get(int,int);
    void set(double,double,double,double);
    void set(const arma::mat &m);
    void add(const arma::mat &m);
    void add(const indexed_array2x2 &m);
    void sub(const arma::mat &m);
    void sub(const indexed_array2x2 &m);
    double data[4];
    int index;
};

#endif
