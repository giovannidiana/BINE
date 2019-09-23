// copied from gibbs4.cxx
// Implemented scheme of algorithm5 in Neal2000
//

#include<iostream>
#include<cstdio>
#define ARMA_NO_DEBUG
#include<armadillo>
#include "../include/routines.h"
#include "../include/indexed_array.h"
#include<fstream>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_sf.h>
#include<cmath>
#include<omp.h>
#include<sys/stat.h>
#include<ctime>

using namespace std;

int main(int argc, char* argv[]){

    if(argc<13){
        cout<<"Bayesian Neural Assembler"<<endl
            <<"usage:"<<endl
            <<"gibbs <NITER> <BURN_IN> <NCELLS> <TIMES> <ASSEMBLIES> <SEED> <LAMBDA0> <LAMBDA1> <ACTIVITY> <READ_BINARY> <outfile> <folder>"<<endl;
        return 0;
    }

    int NITER=atoi(argv[1]);
    int BURN_IN=atoi(argv[2]);
    int N=atoi(argv[3]); // cell number
    int M=atoi(argv[4]); // samples (time frames)
    int Ptrue=atoi(argv[5]); // number of assemblies
    int P=20;
    int PTRACK=10;
    int RNGSEED=atoi(argv[6]);
    double lambda0=atof(argv[7]);
    double lambda1=atof(argv[8]);
    double activity=atof(argv[9]);
    int READ_BINARY=atoi(argv[10]);
    int end_code=-1;
    bool PFIX=false;
    bool verb=false;

    // Create output folder
    string proc_folder(argv[12]);
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
    arma::umat tt(N,3);

    gsl_rng *r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set(r,RNGSEED);

    if(READ_BINARY==0){
        gen_sample4(r,s,t,Ptrue,lambda0,lambda1,activity);
        s.save(proc_folder+"/binary_matrix.dat",arma::raw_ascii);
        t.save(proc_folder+"/membership_orig.dat",arma::raw_ascii);
    } else if(READ_BINARY==1){
        s.load(proc_folder+"/binary_matrix.dat");
        t.load(proc_folder+"/membership_orig.dat");
    } else if(READ_BINARY==2){
        s.load(proc_folder+"/binary_matrix.dat");
        tt.load(proc_folder+"/membership_orig.dat");
        t=tt.col(0);
    }

    // streams
    ofstream FreeEnergyOutput(proc_folder+"/F.dat");
    ofstream pmuOutput(proc_folder+"/pmu.dat");
    ofstream lambda0Output(proc_folder+"/lambda0.dat");
    ofstream lambda1Output(proc_folder+"/lambda1.dat");
    ofstream nOutput(proc_folder+"/n.dat");
    ofstream membershipTraj(proc_folder+"/membership_traj.dat");
    ofstream PTraj(proc_folder+"/P.dat");
    ofstream txsweep(proc_folder+"/txsweep.dat");
    ofstream omegaTraj;
    ofstream lockFile;

    omegaTraj.open(proc_folder+"/omega_traj.dat");

    int i,j,k,l,mu;

    // predicted membership
    int* predicted_membership = new int[N];

    // Order of membership update (to be shuffled)
    int* order = new int[N];
    for(l=0;l<N;l++) order[l]=l;

    // Final matching labels
    int* labels = new int[Ptrue];

    // minimum accuracy
    double pred_accuracy;

    // Decoding success
    bool good;

    // Neuronal membership (updated each sweep)
    int* ids = new int[N];

    // This is needed to decide the average membership
    arma::mat ids_prob(N,PTRACK,arma::fill::zeros);

    unsigned char* omega_write = new unsigned char[M*P];
    vector<int> omega_crypt;
    double meanFreeEnergy=0,F=0;
    double meanTXsweep=0;
    int counter=0;
    double acc;

    // INIT

    // generate parameters
    double alpha_lambda0=1;
    double beta_lambda0=2;
    double alpha_lambda1=2;
    double beta_lambda1=1;
    double alpha_n=10;

    double alpha_p=1;
    double beta_p=2;

    vector< vector<int> > omega(P,vector<int>(M,0));
    vector<int> omega_test(M);

    int mu_old,tmp_mu;

    for(i=0;i<N;i++) ids[i] = gsl_rng_uniform_int(r,P);    

    double prob_binom;
    double p0L,p1L;
    double p_tmp;

    gsl_ran_discrete_t* gen=NULL;

    int sample=0;
    double z;

    arma::mat ids_mat(P,N,arma::fill::zeros);
    vector<int> group_sizes(P,0);
    vector<int> Act(P,0);

    for(i=0;i<N;++i){
        ids_mat(ids[i],i)=1;
        group_sizes[ids[i]]++;
    }

    arma::mat tmp_mat(2,2,arma::fill::zeros);
    arma::vec tmp_vec(2);    
    //arma::cube T(2,2,P,arma::fill::zeros);
    std::vector<indexed_array2x2> T(P,indexed_array2x2(0));
    double varF=0,eta=1;

    for(k=0;k<M;k++){
        for(mu=0;mu<P;mu++){
            Act[mu]+=omega[mu][k];
        }
    }

    mu=0; do{
        if(group_sizes[mu]==0){
            P--;
            if(verb) cout<<"delete "<<mu<<endl;
            for(l=0;l<N;++l){
                if(ids[l]>mu) ids[l]--;
            }

            T.erase(T.begin()+mu);
            omega.erase(omega.begin()+mu);
            Act.erase(Act.begin()+mu);
            group_sizes.erase(group_sizes.begin()+mu);
            mu=0;
        } else mu++;
    } while(mu<group_sizes.size());

    while(sample<NITER){

        sample++;
        cout<<"sample "<<sample<<' '<<P<<' '<<eta<<"         \r"<<flush;

       
	    if(verb){
			for(mu=0;mu<P;mu++) cout<<group_sizes[mu]<<' ';cout<<endl;
		}

        for(mu=0;mu<P;mu++) T[mu].set(beta_lambda0,beta_lambda1,alpha_lambda0,alpha_lambda1);

        for(i=0;i<N;++i){
            for(k=0;k<M;k++){
                T[ids[i]].data[omega[ids[i]][k]+2*s(i,k)]+=1;
            }
        }        

        for(mu=0;mu<P;mu++){
            for(k=0;k<M;k++){
                tmp_mat.fill(0);

                int active=0,inactive=0;
                for(i=0;i<N;i++){
                    if(ids[i]==mu){
                        if(s(i,k)==0) inactive++;
                        else active++;
                    }
                }

                tmp_mat(omega[mu][k],0)=inactive;
                tmp_mat(omega[mu][k],1)=active;

                tmp_vec(0)=tmp_mat(omega[mu][k],0);
                tmp_vec(1)=tmp_mat(omega[mu][k],1);

                T[mu].sub(tmp_mat);

                Act[mu]+=-omega[mu][k];

                p0L=gsl_sf_lnbeta(Act[mu]+alpha_p,M-Act[mu]+beta_p+1)+
                    gsl_sf_lnbeta(T[mu].get(0,1)+tmp_vec(1),T[mu].get(0,0)+tmp_vec(0))+
                    gsl_sf_lnbeta(T[mu].get(1,1),T[mu].get(1,0));
                p1L=gsl_sf_lnbeta(Act[mu]+alpha_p+1,M-Act[mu]+beta_p)+
                    gsl_sf_lnbeta(T[mu].get(0,1),T[mu].get(0,0))+
                    gsl_sf_lnbeta(T[mu].get(1,1)+tmp_vec(1),T[mu].get(1,0)+tmp_vec(0));
                prob_binom=1/(1+exp(p0L-p1L));
                omega[mu][k]=(gsl_rng_uniform(r)<prob_binom) ? 1 : 0;
                Act[mu]+=omega[mu][k];
                T[mu].data[omega[mu][k]]+=tmp_vec(0);
                T[mu].data[omega[mu][k]+2]+=tmp_vec(1);
            }
        }

        counter=0; 


        gsl_ran_shuffle (r, order, N, sizeof (int));

		for(j=0;j<N;++j){            

            i=order[j];

            vector<indexed_array2x2> tmp_cube(P,indexed_array2x2(0));

            for(k=0;k<M;k++){
                for(mu=0;mu<P;mu++){
                    tmp_cube[mu].data[omega[mu][k]+2*s(i,k)]+=1;
                }
            }

            mu_old=ids[i];
            group_sizes[mu_old]+=-1;
            T[mu_old].sub(tmp_cube[mu_old]);

            double* tmp_prob = new double[P+1];
            for(mu=0;mu<P;++mu) tmp_prob[mu]=group_sizes[mu];
            tmp_prob[P]=alpha_n;

            gen = gsl_ran_discrete_preproc(P+1, tmp_prob);

            tmp_mu=gsl_ran_discrete(r,gen);
            gsl_ran_discrete_free(gen);
            gen=NULL;

            if(tmp_mu==mu_old){
                T[ids[i]].add(tmp_cube[ids[i]]);
            } else if(tmp_mu<P){
                acc=//terms due to the new group
                    gsl_sf_lnbeta(T[tmp_mu].get(0,1)+tmp_cube[tmp_mu].get(0,1),T[tmp_mu].get(0,0)+tmp_cube[tmp_mu].get(0,0))-
                     gsl_sf_lnbeta(T[tmp_mu].get(0,1),T[tmp_mu].get(0,0))+
                    gsl_sf_lnbeta(T[tmp_mu].get(1,1)+tmp_cube[tmp_mu].get(1,1),T[tmp_mu].get(1,0)+tmp_cube[tmp_mu].get(1,0))-
                     gsl_sf_lnbeta(T[tmp_mu].get(1,1),T[tmp_mu].get(1,0))+
                    // terms due to the changes in the old group
                    gsl_sf_lnbeta(T[mu_old].get(0,1),T[mu_old].get(0,0))-
                     gsl_sf_lnbeta(T[mu_old].get(0,1)+tmp_cube[mu_old].get(0,1),T[mu_old].get(0,0)+tmp_cube[mu_old].get(0,0))+
                    gsl_sf_lnbeta(T[mu_old].get(1,1),T[mu_old].get(1,0))-
                     gsl_sf_lnbeta(T[mu_old].get(1,1)+tmp_cube[mu_old].get(1,1),T[mu_old].get(1,0)+tmp_cube[mu_old].get(1,0));

                if(gsl_rng_uniform(r)<exp(acc)){
					counter++;
                    ids[i]=tmp_mu;
                    if(verb) cout<<i<<"("<<mu_old<<")->"<<"("<<tmp_mu<<")"<<endl;
                }

                T[ids[i]].add(tmp_cube[ids[i]]);

            } else if(!PFIX && tmp_mu==P) {
                p_tmp=gsl_ran_beta(r,alpha_p,beta_p);
                for(k=0;k<M;k++) omega_test[k]=(gsl_rng_uniform(r)<p_tmp) ? 1 : 0;
                //for(k=0;k<M;k++) omega_test[k]=s(i,k);

                tmp_mat.fill(0);
                for(k=0;k<M;k++) tmp_mat(omega_test[k],s(i,k))+=1;
                acc=
                        //terms due to the new group
                        gsl_sf_lnbeta(tmp_mat(0,1)+alpha_lambda0,tmp_mat(0,0)+beta_lambda0)-
                        gsl_sf_lnbeta(alpha_lambda0,beta_lambda0)+
                        gsl_sf_lnbeta(tmp_mat(1,1)+alpha_lambda1,tmp_mat(1,0)+beta_lambda1)-
                        gsl_sf_lnbeta(alpha_lambda1,beta_lambda1)+
                        // terms due to the changes in the old group
                        gsl_sf_lnbeta(T[mu_old].get(0,1),T[mu_old].get(0,0))-
                        gsl_sf_lnbeta(T[mu_old].get(0,1)+tmp_cube[mu_old].get(0,1),T[mu_old].get(0,0)+tmp_cube[mu_old].get(0,0))+
                        gsl_sf_lnbeta(T[mu_old].get(1,1),T[mu_old].get(1,0))-
                        gsl_sf_lnbeta(T[mu_old].get(1,1)+tmp_cube[mu_old].get(1,1),T[mu_old].get(1,0)+tmp_cube[mu_old].get(1,0));

                if(gsl_rng_uniform(r)<exp(acc)){
					counter++;
                    P+=1;
                    omega.push_back(omega_test);

                    T.push_back(indexed_array2x2(0));
                    T[P-1].set(beta_lambda0,beta_lambda1,alpha_lambda0,alpha_lambda1);

                    T[P-1].add(tmp_mat);

                    group_sizes.push_back(0);
                    Act.push_back(0);

                    for(k=0;k<M;k++) Act[P-1]+=omega_test[k];
                    ids[i]=P-1;
                    if(verb) cout<<i<<"("<<mu_old<<")->"<<"("<<P-1<<")*"<<endl;
                } else {
                    T[ids[i]].add(tmp_cube[ids[i]]);
                }

            }

            group_sizes[ids[i]]+=1;
            if(!PFIX && group_sizes[mu_old]==0){

                P--;
                if(verb) cout<<"delete "<<mu_old<<endl;
                for(l=0;l<N;++l){
                    if(ids[l]>mu_old) ids[l]--;
                }

                T.erase(T.begin()+mu_old);
                omega.erase(omega.begin()+mu_old);
                Act.erase(Act.begin()+mu_old);
                group_sizes.erase(group_sizes.begin()+mu_old);

            }

            delete [] tmp_prob;

        }


//        // update n
//        group_sizes=arma::sum(part.ids_mat,1);
//        part.n=dirichlet(r,alpha_n+group_sizes);

//        double* pgen = new double[P];
//        double* lambda0gen = new double[P];
//        double* lambda1gen = new double[P];

//        //pmu
//        for(mu=0;mu<P;++mu){
//            pgen[mu] = gsl_ran_beta(r,alpha_p+Act[mu],beta_p+M-Act[mu]);
//        }

//        //lambda
//        for(mu=0;mu<P;++mu){
//            lambda0gen[mu] = gsl_ran_beta(r,alpha_lambda0+T[mu].get(1,0),beta_lambda0+T[mu].get(0,0));
//            lambda1gen[mu] = gsl_ran_beta(r,alpha_lambda1+T[mu].get(1,1),beta_lambda1+T[mu].get(0,1));
//        }

        if(sample>BURN_IN){

            F=0;

            for(mu=0;mu<P;mu++){
                //F+= gsl_sf_lngamma(group_sizes(mu))-gsl_sf_lngamma(alpha_n/1000.);
                F+= gsl_sf_lnbeta(Act[mu]+alpha_p,M-Act[mu]+beta_p)-gsl_sf_lnbeta(alpha_p,beta_p);
                F+= gsl_sf_lnbeta(T[mu].get(0,1),T[mu].get(0,0))-gsl_sf_lnbeta(alpha_lambda0,beta_lambda0);
                F+= gsl_sf_lnbeta(T[mu].get(1,1),T[mu].get(1,0))-gsl_sf_lnbeta(alpha_lambda1,beta_lambda1);
            }

            for(i=0;i<N;++i){
                if(ids[i]<PTRACK){
                    ids_prob(i,ids[i])++;
                }
            }

            meanFreeEnergy+=F;
            meanTXsweep+=(double)counter;

            if(sample-BURN_IN>20){
                varF+=pow(F-meanFreeEnergy/(sample-BURN_IN),2);
                eta=sqrt(varF/(sample-BURN_IN-20))*(sample-BURN_IN)/fabs(meanFreeEnergy);
                //if(eta < 5e-3) break;
            }

            FreeEnergyOutput<<F<<endl;
            txsweep<<(double)counter/N<<endl;

            for(mu=0;mu<min(PTRACK,P);mu++) {
				for(k=0;k<M;k++) omegaTraj<<omega[mu][k]<<' '; omegaTraj<<endl;
			}
			while(mu<PTRACK){
				for(k=0;k<M;k++) omegaTraj<<"NA"<<' '; omegaTraj<<endl;
				mu++;
			}
            
//            for(mu=0;mu<P;mu++) pmuOutput<<pgen[mu]<<' ';
//            while(mu<PTRACK){
//                pmuOutput<<"NA"<<' ';
//                mu++;
//            }
//            pmuOutput<<endl;

//            for(mu=0;mu<P;mu++) nOutput<<part.n(mu)<<' ';
//            nOutput<<endl;

            for(i=0;i<N;++i) membershipTraj<<ids[i]<<' ';
            membershipTraj<<endl;

//            for(mu=0;mu<P;mu++) lambda0Output<<lambda0gen[mu]<< ' ';
//            while(mu<PTRACK) {
//                lambda0Output<<"NA"<<' ';
//                mu++;
//            }
//            lambda0Output<<endl;

//            for(mu=0;mu<P;mu++) lambda1Output<<lambda1gen[mu]<< ' ';
//            while(mu<PTRACK) {
//                lambda1Output<<"NA"<<' ';
//                mu++;
//            }
//            lambda1Output<<endl;

            PTraj<<P<<endl;

//            if(lambda1gen!=NULL) delete [] lambda1gen;
//            if(lambda0gen!=NULL) delete [] lambda0gen;
//            if(pgen!=NULL) delete [] pgen;
        }

    }

    ofstream OutputMembership(proc_folder+"/membership_gibbs.dat");

    for(i=0;i<N;++i){
        predicted_membership[i]=ids_prob.row(i).index_max();
        OutputMembership<<predicted_membership[i]<<std::endl;
    }

    good=match_prediction2(t,ids,N,P,labels,pred_accuracy);

//    // Safe writing to output using a lock file.

    // standby until lock file desappear
    while(fileExists("lock")){};

    // lock
    lockFile.open("lock");
    lockFile<<1<<endl;
    lockFile.close();

    //write
    ofstream stdOutput(argv[11],fstream::app);
    stdOutput<<RNGSEED<<' '<<P<<' '<<lambda0<<' '<<lambda1<<' '<<activity<<' '<<meanFreeEnergy/(sample-BURN_IN)<<' '<<meanTXsweep/(sample-BURN_IN)<<' '<<pred_accuracy<<' '<<good<<endl;
    stdOutput.close();

    //unlock
    remove("lock");
    FreeEnergyOutput.close();
    pmuOutput.close();
    lambda0Output.close();
    lambda1Output.close();
    nOutput.close();
    membershipTraj.close();
    omegaTraj.close();

    //free memory
    delete [] predicted_membership;
    delete [] order;
    delete [] labels;
    delete [] ids;
    delete [] omega_write;


    return 0;
}
