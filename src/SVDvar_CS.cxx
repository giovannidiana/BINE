// SVD/promax Ensemble detection
//
//

#include<iostream>
#include<cstdio>
#define ARMA_NO_DEBUG
#include<armadillo>
#include "../include/routines.h"
//#include "../include/crypt.h"
#include "../include/indexed_array.h"
#include "../include/vari.h"
#include<fstream>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_sf.h>
#include<cmath>
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
    int P=0;
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

    int i,j,k,l,mu;

    // Final matching labels
    int* labels = new int[Ptrue];

    // minimum accuracy
    double pred_accuracy;

    // Decoding success
    bool good;

    // random shift
    int random_shift;

    // null-model eigenvectors
    arma::vec lambda_null(N,arma::fill::zeros);
    arma::vec lambda_null_sq(N,arma::fill::zeros);


    arma::mat U,V;
    arma::vec lambda;

    arma::mat sd=arma::conv_to<arma::mat>::from(s);
    arma::mat sd_shuffle=sd;
    int Nshuffles=20;

    for(int counter=0;counter<Nshuffles;counter++){
        std::cout<<"shuffle number "<<counter<<"     \r";
        for(i=0;i<N;i++){
            random_shift=gsl_rng_uniform_int(r,200);
            for(k=random_shift;k<sd.n_cols;k++) sd_shuffle(i,k)=sd(i,k-random_shift);
            for(k=0;k<random_shift;k++) sd_shuffle(i,k)=sd(i,sd.n_cols-random_shift+k);
        }
        arma::svd(U,lambda,V,sd_shuffle);
        lambda_null = lambda_null+lambda;
        lambda_null_sq = lambda_null_sq+pow(lambda,2);
    }

    lambda_null /= Nshuffles;
    lambda_null_sq /= Nshuffles;
    lambda_null_sq = sqrt(lambda_null_sq-pow(lambda_null,2));

    arma::mat Urot;
    arma::svd(U,lambda,V,sd);
    vector<arma::uword> cell_labels(N);
    vector<int> cell_membership(N);
    vector<double> cell_select(N);

    P=0;
    for(i=0;i<lambda.n_elem;++i){
        if(lambda(i)>lambda_null(i)+2*lambda_null_sq(i)) P++;
        else break;
    }
    if(P==0) P++;

    std::cout<<"using "<<P<<" ensembles."<<std::endl;
    arma::mat Ucut=U.cols(0,P-1);
    Vari rot(Ucut);

    rot.varimax(1e-10,true);
    Urot = pow(Ucut * rot.rot,2);

    for(j=0;j<N;++j) {
        Urot.row(j).max(cell_labels[j]);
		if(arma::stddev(Urot.row(j))==0) cout<<"standard deviation = 0"<<endl;
        cell_select[j]=(Urot.row(j).max()-arma::mean(Urot.row(j)))/arma::stddev(Urot.row(j));
        cell_membership[j]=cell_labels[j];
    }

    good=match_prediction2(t,&cell_membership[0],N,P,labels,pred_accuracy);

    //write
    ofstream stdOutput(argv[9],fstream::app);
    for(i=0;i<N;++i) membershipOut<<cell_membership[i]<<' '<<cell_select[i]<<endl;
    stdOutput<<P<<' '<<lambda0<<' '<<lambda1<<' '<<activity<<' '<<pred_accuracy<<' '<<good<<endl;
    stdOutput.close();

    //free memory
    delete [] labels;

    return 0;
}
