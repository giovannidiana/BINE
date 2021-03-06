// kmeans Ensemble detection
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
#include<sys/stat.h>
#include<ctime>

using namespace std;

int main(int argc, char* argv[]){

    if(argc<9){
        cout<<"Generate surrogate population activity"<<endl
            <<"usage:"<<endl
            <<"./generate_data <NCELLS> <TIMES> <ASSEMBLIES> <SEED> <LAMBDA0> <LAMBDA1> <ACTIVITY> <outfile> <lastIsFree>"<<endl;
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
    bool isLastFree=argc>9;
    
    bool verb=false;

    // Create output folder
    string proc_folder(argv[8]);
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
    gsl_rng_set(r,RNGSEED+(int)(100*lambda0));

    //gen_sample4(r,s,t,Ptrue,lambda0,lambda1,activity);
    //gen_sample5(r,s,t,Ptrue);

    gen_sample6(r,s,t,Ptrue,lambda0,lambda1,activity,isLastFree); // this allows more variable ensembles
    s.save(proc_folder+"/binary_matrix.dat",arma::raw_ascii);
    t.save(proc_folder+"/membership_orig.dat",arma::raw_ascii);
    
    return 0;
}
