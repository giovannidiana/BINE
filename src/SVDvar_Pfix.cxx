// SVD/promax Ensemble detection
//
//

#include<iostream>
#include<cstdio>
#define ARMA_NO_DEBUG
#include<armadillo>
#include "../include/routines.h"
#include "../include/indexed_array.h"
#include "../include/vari.h"
#include "../include/entropy.h"
#include<fstream>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_sf.h>
#include<cmath>
#include<omp.h>
#include<sys/stat.h>
#include<ctime>

using namespace std;

int main(int argc, char* argv[]){

    if(argc<11){
        cout<<"Bayesian Neural Assembler"<<endl
            <<"usage:"<<endl
            <<"gibbs <NCELLS> <TIMES> <ASSEMBLIES> <SEED> <LAMBDA0> <LAMBDA1> <ACTIVITY> <READ_BINARY> <outfile> <folder>"<<endl;
        return 0;
    }

    int N=atoi(argv[1]); // cell number
    int M=atoi(argv[2]); // samples (time frames)
    int Ptrue=atoi(argv[3]); // number of assemblies
    int P=Ptrue;
    int RNGSEED=atoi(argv[4]);
    double lambda0=atof(argv[5]);
    double lambda1=atof(argv[6]);
    double activity=atof(argv[7]);
    int READ_BINARY=atoi(argv[8]);
    bool verb=false;

    // Create output folder
    string proc_folder(argv[10]);
    struct stat sb;
    if (stat(proc_folder.c_str(), &sb) != 0){
        const int dir_err = mkdir(proc_folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err){
            printf("Error creating directory!");
            exit(1);
        }
    }

    arma::imat s(N,M);
    arma::uvec t(N);

    gsl_rng *r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set(r,RNGSEED);

    if(READ_BINARY==0){
        gen_sample4(r,s,t,Ptrue,lambda0,lambda1,activity);
        s.save(proc_folder+"/binary_matrix.dat",arma::raw_ascii);
        t.save(proc_folder+"/membership_orig.dat",arma::raw_ascii);
    } else {
        s.load(proc_folder+"/binary_matrix.dat");
        t.load(proc_folder+"/membership_orig.dat");
    }

    // streams
    ofstream membershipOut(proc_folder+"/membership_SVD.dat");
    ofstream entropyOut(proc_folder+"/entropy.dat");

    int i,j,k,l,mu;

    // Final matching labels
    int* labels = new int[Ptrue];

    // minimum accuracy
    double pred_accuracy;

    // Decoding success
    bool good;

    double enMin;

    arma::mat sd=arma::conv_to<arma::mat>::from(s);
    arma::mat U,V;
    arma::vec lambda;
    arma::mat Urot;
    arma::svd(U,lambda,V,sd);
    vector<arma::uword> cell_labels(N);
    vector<int> cell_membership(N);
    int nguess;
    int maxguess=10;

    vector<double> enVec_mean(maxguess);
    vector<double> enVec_max(maxguess);

    arma::mat Ucut=U.cols(0,P-1);
    Vari rot(Ucut);
    rot.varimax(1e-10,true);
    Urot = pow(Ucut * rot.rot,2);
    cout<<rot.rot<<endl;

    for(j=0;j<N;++j) {
        Urot.row(j).max(cell_labels[j]);
        cell_membership[j]=cell_labels[j];
    }

    good=match_prediction2(t,&cell_membership[0],N,P,labels,pred_accuracy);

    //write
    ofstream stdOutput(argv[9],fstream::app);
    for(i=0;i<N;++i) membershipOut<<cell_membership[i]<<endl;
    for(i=0;i<maxguess;++i) entropyOut<<i+2<<' '<<enVec_mean[i]<<' '<<enVec_max[i]<<endl;
    stdOutput<<P<<' '<<lambda0<<' '<<lambda1<<' '<<activity<<' '<<pred_accuracy<<' '<<good<<endl;
    stdOutput.close();

    //free memory
    delete [] labels;

    return 0;
}
