#include<iostream>
#include<cstdio>
#define ARMA_NO_DEBUG
#include<armadillo>
#include "../include/routines.h"
#include "../include/indexed_array.h"
#include "../include/ReadLine.h"
#include<fstream>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_sf.h>
#include<cmath>
#include<sys/stat.h>
#include<ctime>
#include <getopt.h>

using namespace std;

int main(int argc, char* argv[]){

	int c;

    arma::imat s_raw; // binary activity s_{ik}
	string proc_folder;
    int NITER;
    int BURN_IN;
    int TRIM=0;
    int N;
    int M;
    int P; // initial number of assemblies
    int i,j,k;
    int THRESH = 0;
	int THRESH2= 0;

    int RNGSEED=0;
    bool verb=false;
    bool contd=false;
    int MAXP=100;
	int options_required[]={1,1,0,1,0,1,0,0,0,0,0,0};
	int options_provided[]={0,0,0,0,0,0,0,0,0,0,0,0};

    static struct option long_options[] = {
		      {"folder",  required_argument, NULL,  'f' },
			  {"niter",  required_argument, NULL,  'i' },
			  {"trim",  required_argument, NULL,  't' },
			  {"assemblies", required_argument, NULL,  'a' },
			  {"seed",  required_argument, NULL,  's' },
			  {"file",  required_argument, NULL,  'b' },
			  {"min_neur", required_argument,  NULL,  '1' },
			  {"min_act", required_argument,  NULL,  '2' },
			  {"continue", no_argument, NULL,  'c' },
			  {"burn_in", required_argument, NULL, 'u'},
			  {"verbose", no_argument,NULL,'v'},
			  {NULL,         0,                NULL,  0 }
	};

	while(true){
		// index along long_options 
		int option_index = 0;
		
		c = getopt_long(argc, argv, "", long_options, &option_index);

		options_provided[option_index]=1;

		switch(c) {
			case 'f':
			// Create output folder
			proc_folder.assign(optarg);
			struct stat sb;
			if (stat(proc_folder.c_str(), &sb) != 0){
				const int dir_err = mkdir(proc_folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
				if (-1 == dir_err){
					printf("Error creating directory!");
					exit(1);
				}
			}
			break;

			case 'b':
			s_raw.load(optarg);
			break;

			case 'i':
			NITER=atoi(optarg);
			break;

			case 's':
			RNGSEED=atoi(optarg);
			break;

			case 'a':
			P=atoi(optarg);
			break;

			case 'u':
			BURN_IN=atoi(optarg);
			break;

			case 't':
			TRIM==atoi(optarg);
			break;

			case 'c':
			contd=(strcmp(optarg,"1")==0);
			break;

			case '1':
			THRESH=atoi(optarg);
			break;

			case '2':
			THRESH2=atoi(optarg);
			break;
			
			case 'v':
			verb=true;
			break;

			default:
			cerr<<"Bayesian Inference of Neuronal Assemblies"<<endl
			    <<"G. Diana, T. Sainsbury, M. Mayer"<<endl
				<<"bioRxiv 452557; doi: https://doi.org/10.1101/452557"<<endl 
				<<"usage:"<<endl
				<<argv[0]<<endl;
			return 1;
		}
    }

	// check all required options are there
	for(i=0;i<12;i++){
		if(options_required[i]==1 && options_provided[i]==0){
			cerr<<"! missing "<<long_options[i].name<<endl;
			return 1;
		}
	}

    arma::imat s;
    arma::uvec cell_selected;

    gsl_rng *r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set(r,RNGSEED);

//  Restrict to synchronous activity matrix s(N,M)
    filter_binary2(s_raw,s,cell_selected,THRESH,THRESH2);

    N=s.n_rows;
    M=s.n_cols;

    cout<<N<<' '<<M<<endl;

    // streams
    ofstream FreeEnergyOutput;
    ofstream pmuOutput;
    ofstream lambda0Output;
    ofstream lambda1Output;
    ofstream nOutput;
    ofstream membershipTraj;
    ofstream PTraj;
    ofstream txsweep;
    ofstream omegaTraj;
    ofstream lockFile;

    int ids[N];
    if(contd){
        rline(proc_folder+"/membership_traj.dat",ids,N);
		P=*max_element(ids,ids+N)+1;
    } else {
        for(i=0;i<N;i++) ids[i] = gsl_rng_uniform_int(r,P);
    }

    if(!contd){
        FreeEnergyOutput.open(proc_folder+"/F.dat");
        pmuOutput.open(proc_folder+"/pmu.dat");
        lambda0Output.open(proc_folder+"/lambda0.dat");
        lambda1Output.open(proc_folder+"/lambda1.dat");
        nOutput.open(proc_folder+"/n.dat");
        membershipTraj.open(proc_folder+"/membership_traj.dat");
        PTraj.open(proc_folder+"/P.dat");
        txsweep.open(proc_folder+"/txsweep.dat");
    	omegaTraj.open(proc_folder+"/omega_traj.dat");
    } else {
        FreeEnergyOutput.open(proc_folder+"/F.dat",ios_base::app);
        pmuOutput.open(proc_folder+"/pmu.dat",ios_base::app);
        lambda0Output.open(proc_folder+"/lambda0.dat",ios_base::app);
        lambda1Output.open(proc_folder+"/lambda1.dat",ios_base::app);
        nOutput.open(proc_folder+"/n.dat",ios_base::app);
        membershipTraj.open(proc_folder+"/membership_traj.dat",ios_base::app);
        PTraj.open(proc_folder+"/P.dat",ios_base::app);
        txsweep.open(proc_folder+"/txsweep.dat",ios_base::app);
		omegaTraj.open(proc_folder+"/omega_traj.dat",ios_base::app);
    }

    ofstream dimensionOutput(proc_folder+"/dims.dat");
    ofstream cellselectFile(proc_folder+"/selection.dat");

    // write selection
    for(unsigned int i=0;i<cell_selected.n_elem;i++) cellselectFile<<cell_selected(i)<<endl;
    cellselectFile.close();

    // write dimensions
    dimensionOutput<<N<<' '<<M<<' '<<P<<' '<<NITER<<endl;
    dimensionOutput.close();


    // output variables
    unsigned char omega_write[M*P];
    vector<int> omega_crypt;    
    double F;
    int counter=0;
    double acc;

    // INIT

    // generate parameters
    double alpha_lambda0=1;
    double beta_lambda0=5;
    double alpha_lambda1=10;
    double beta_lambda1=1;
	// Fish analysis done with alpha_n=0.1
    double alpha_n=1;

    double alpha_p=1;
    double beta_p=5;


    vector< vector<int> > omega(P,vector<int>(M,0));
    vector<int> omega_test(M);

    int mu,nu;
    int mu_old,tmp_mu;

    double prob_binom;
    double p0L,p1L;
    double p_tmp;

    gsl_ran_discrete_t* gen=NULL;

    int sample=0;

    arma::mat ids_mat(P,N,arma::fill::zeros);
    vector<int> group_sizes(P,0);
    vector<int> Act(P,0);

    for(i=0;i<N;++i){
        ids_mat(ids[i],i)=1;
        group_sizes[ids[i]]++;
    }

    arma::mat tmp_mat(2,2,arma::fill::zeros);
    arma::vec tmp_vec(2);    
    std::vector<indexed_array2x2> T(P,indexed_array2x2(0));

    for(k=0;k<M;k++){
        for(mu=0;mu<P;mu++){
            Act[mu]+=omega[mu][k];
        }
    }

    mu=0; do{
        if(group_sizes[mu]==0){
            P--;
            if(verb) cout<<"delete "<<mu<<endl;
            for(unsigned int l=0;l<N;++l){
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
        cout<<"sample "<<sample<<' '<<P<<"         \r"<<flush;

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
        int order[N];
        for(unsigned int l=0;l<N;l++) order[l]=l;
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

            } else if(tmp_mu==P) {
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
            if(group_sizes[mu_old]==0){

                P--;
                if(verb) cout<<"delete "<<mu_old<<endl;
                for(unsigned int l=0;l<N;++l){
                    if(ids[l]>mu_old) ids[l]--;
                }

                T.erase(T.begin()+mu_old);
                omega.erase(omega.begin()+mu_old);
                Act.erase(Act.begin()+mu_old);
                group_sizes.erase(group_sizes.begin()+mu_old);

            }

            delete [] tmp_prob;

        }

        if(sample>BURN_IN && sample%TRIM==0){

            F=0;
            for(mu=0;mu<P;mu++){
                //F+= gsl_sf_lngamma(group_sizes(mu))-gsl_sf_lngamma(alpha_n/1000.);
                F+= gsl_sf_lnbeta(Act[mu]+alpha_p,M-Act[mu]+beta_p)-gsl_sf_lnbeta(alpha_p,beta_p);
                F+= gsl_sf_lnbeta(T[mu].get(0,1),T[mu].get(0,0))-gsl_sf_lnbeta(alpha_lambda0,beta_lambda0);
                F+= gsl_sf_lnbeta(T[mu].get(1,1),T[mu].get(1,0))-gsl_sf_lnbeta(alpha_lambda1,beta_lambda1);
            }

            FreeEnergyOutput<<F<<endl;
            txsweep<<(double)counter/N<<endl;

            for(mu=0;mu<min(MAXP,P);mu++) {
				for(k=0;k<M;k++) omegaTraj<<omega[mu][k]<<' '; omegaTraj<<endl;
			}
			while(mu<MAXP){
				for(k=0;k<M;k++) omegaTraj<<"NA"<<' '; omegaTraj<<endl;
				mu++;
			}
            
            double* pgen       = new double[MAXP]; std::fill(pgen,pgen+MAXP,0);
            double* lambda0gen = new double[MAXP]; std::fill(lambda0gen,lambda0gen+MAXP,0);
            double* lambda1gen = new double[MAXP]; std::fill(lambda1gen,lambda1gen+MAXP,0);

            // group sizes
            for(mu=0;mu<min(P,MAXP);mu++) nOutput<<group_sizes[mu]<<' ';
            while(mu<MAXP) {
                nOutput<<"NA"<<' ';
                mu++;
            }
            nOutput<<endl;

            //pmu
            for(mu=0;mu<min(MAXP,P);++mu) pgen[mu] = gsl_ran_beta(r,alpha_p+Act[mu],beta_p+M-Act[mu]);
            for(mu=0;mu<min(P,MAXP);mu++) pmuOutput<<pgen[mu]<<' ';
            while(mu<MAXP){
                pmuOutput<<"NA"<<' ';
                mu++;
            }
            pmuOutput<<endl;

            //lambda
            for(mu=0;mu<min(MAXP,P);++mu){
				// before 30/10/2018. This version seems incorrect due due to swapping the arguments of the function .get(omega,s).  
                //lambda0gen[mu] = gsl_ran_beta(r,alpha_lambda0+T[mu].get(1,0),beta_lambda0+T[mu].get(0,0));
                //lambda1gen[mu] = gsl_ran_beta(r,alpha_lambda1+T[mu].get(1,1),beta_lambda1+T[mu].get(0,1));
				// 9/11/2018. removed the alpha's and beta's because already incorporated in the definitions of T[mu]
                lambda0gen[mu] = gsl_ran_beta(r,T[mu].get(0,1),T[mu].get(0,0));
                lambda1gen[mu] = gsl_ran_beta(r,T[mu].get(1,1),T[mu].get(1,0));
            }

// test something
            

            for(mu=0;mu<min(P,MAXP);mu++) lambda0Output<<lambda0gen[mu]<< ' ';
            for(mu=0;mu<min(P,MAXP);mu++) lambda1Output<<lambda1gen[mu]<< ' ';

            while(mu<MAXP) {
                lambda0Output<<"NA"<<' ';
                lambda1Output<<"NA"<<' ';
                mu++;
            }
            lambda0Output<<endl;
            lambda1Output<<endl;

            // membership
            for(i=0;i<N;++i) membershipTraj<<ids[i]<<' ';
            membershipTraj<<endl;

            // Number of assemblies
            PTraj<<P<<endl;

            delete [] lambda1gen;
            delete [] lambda0gen;
            delete [] pgen;
        }

    }

    FreeEnergyOutput.close();
    pmuOutput.close();
    lambda0Output.close();
    lambda1Output.close();
    nOutput.close();
    membershipTraj.close();
    omegaTraj.close();

    return 0;
}
