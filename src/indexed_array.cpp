#include"../include/indexed_array.h"

indexed_array2x2::indexed_array2x2(int ind) :
    index(ind)
{
    for(unsigned int i=0;i<4;++i) data[i]=0;
}

double indexed_array2x2::get(int i,int j){
    double value=data[i+2*j];
    return(value);
}

void indexed_array2x2::set(double x1,double x2,double x3,double x4){
    data[0]=x1;
    data[1]=x2;
    data[2]=x3;
    data[3]=x4;
}

void indexed_array2x2::add(const arma::mat &m){
    data[0]+=m(0,0);
    data[1]+=m(1,0);
    data[2]+=m(0,1);
    data[3]+=m(1,1);
}

void indexed_array2x2::add(const indexed_array2x2 &m){
    for(unsigned int i=0;i<4;++i) data[i]+=m.data[i];
}
void indexed_array2x2::sub(const arma::mat &m){
    data[0]-=m(0,0);
    data[1]-=m(1,0);
    data[2]-=m(0,1);
    data[3]-=m(1,1);
}

void indexed_array2x2::sub(const indexed_array2x2 &m){
    for(unsigned int i=0;i<4;++i) data[i]-=m.data[i];
}

void indexed_array2x2::set(const arma::mat &m){
    data[0]=m(0,0);
    data[1]=m(1,0);
    data[2]=m(0,1);
    data[3]=m(1,1);
}
