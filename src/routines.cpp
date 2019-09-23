// Many bugs fixed, everything works.

#include<iostream>
#include <random>
#include <fstream>
#include "../include/routines.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <sys/stat.h>
#include "../include/vari.h"

using namespace std;

void gen_sample(arma::mat &m, arma::uvec &t, int P,int seed=0){
	m.fill(0);
	int i,k;
	int ncells=m.n_rows;
	int ntimes=m.n_cols;
	double lambda[2] = {1e-5,1.-1e-5};
	arma::vec n(P); n.fill(1./P);
    arma::vec p(P); p.fill(.1);
	int omega[P]; 
	std::default_random_engine generator;
	std::uniform_int_distribution<int> rint(0,P-1);
	std::uniform_real_distribution<double> runif(0.0,1.0);
	arma::arma_rng::set_seed(seed);

	for(i=0;i<ncells;++i){
		t(i) = rint(generator);
	}

	for(k=0;k<ntimes;++k){
		for(i=0;i<P;++i) omega[i] = (runif(generator) < p[i]) ? 1 : 0;
		for(i=0;i<ncells;++i){
			if(omega[t(i)]==1){
				if(runif(generator)< lambda[1]) m(i,k) = 1;
			} else { 
				if(runif(generator)<lambda[0]) m(i,k) =1;
			}
		}
	}
}

void gen_sample2(arma::mat &m, arma::uvec &t, int P, int seed=0){
	m.fill(0);
	int i,k;
	int ncells=m.n_rows;
	int ntimes=m.n_cols;
    double eps=1e-1;
	double lambda[2] = {eps,1-eps};
	arma::vec n(P); n.fill(1./P);
    arma::vec p(P); p.fill(.01);
	int omega[P]; 
	std::default_random_engine generator;
	std::uniform_int_distribution<int> rint(0,P-1);
	std::uniform_real_distribution<double> runif(0.0,1.0);
	arma::arma_rng::set_seed(seed);
    
	int label=0;
	for(i=0;i<ncells;++i){
		t(i) = label;
		if(i>0 && i%(ncells/P)==0 && label<P-1) label++;
	}

	for(k=0;k<ntimes;++k){
		for(i=0;i<P;++i) omega[i] = (runif(generator) < p[i]) ? 1 : 0;
		for(i=0;i<ncells;++i){
			if(omega[t(i)]==1){
				if(runif(generator)< lambda[1]) m(i,k) = 1;
			} else { 
				if(runif(generator)<lambda[0]) m(i,k) =1;
			}
		}
	}
}

void gen_sample3(gsl_rng* r, arma::fmat &m, arma::uvec &t, int P,double eps, double pmu){
    m.fill(0);
    int i,k;
    int ncells=m.n_rows;
    int ntimes=m.n_cols;
    double lambda[2] = {eps,1-eps};
    arma::vec n(P); n.fill(1./P);
    arma::vec p(P); p.fill(pmu);
    int omega[P];

    int label=0;
    for(i=0;i<ncells;++i){
        t(i) = label;
        if(i>0 && i%(ncells/P)==0 && label<P-1) label++;
    }

    for(k=0;k<ntimes;++k){
        for(i=0;i<P;++i) omega[i] = (gsl_rng_uniform(r) < p[i]) ? 1 : 0;
        for(i=0;i<ncells;++i){
            if(omega[t(i)]==1){
                if(gsl_rng_uniform(r)< lambda[1]) m(i,k) = 1;
            } else {
                if(gsl_rng_uniform(r)<lambda[0]) m(i,k) = 1;
            }
        }
    }
}

void gen_sample4(gsl_rng* r, arma::fmat &m, arma::uvec &t, int P,double lambda0,double lambda1, double pmu){
    m.fill(0);
    int i,k;
    int ncells=m.n_rows;
    int ntimes=m.n_cols;
    double lambda[2] = {lambda0,lambda1};
    arma::vec n(P); n.fill(1./P);
    arma::vec p(P); p.fill(pmu);
    int omega[P];

    int label=0;
    for(i=0;i<ncells;++i){
        t(i) = label;
        if(i>0 && i%(ncells/P)==0 && label<P-1) label++;
    }

    for(k=0;k<ntimes;++k){
        for(i=0;i<P;++i) omega[i] = (gsl_rng_uniform(r) < p[i]) ? 1 : 0;
        for(i=0;i<ncells;++i){
            if(omega[t(i)]==1){
                if(gsl_rng_uniform(r)< lambda[1]) m(i,k) = 1;
            } else {
                if(gsl_rng_uniform(r)<lambda[0]) m(i,k) = 1;
            }
        }
    }
}

void gen_sample4(gsl_rng* r, arma::imat &m, arma::uvec &t, int P,double lambda0,double lambda1, double pmu){
    m.fill(0);
    int i,k;
    int ncells=m.n_rows;
    int ntimes=m.n_cols;
    double lambda[2] = {lambda0,lambda1};
    arma::vec n(P); n.fill(1./P);
    arma::vec p(P); p.fill(pmu);
    int omega[P];

    int label=0;
    for(i=0;i<ncells;++i){
        t(i) = label;
        if(i>0 && i%(ncells/P)==0 && label<P-1) label++;
    }

    for(k=0;k<ntimes;++k){
        for(i=0;i<P;++i) omega[i] = (gsl_rng_uniform(r) < p[i]) ? 1 : 0;
        for(i=0;i<ncells;++i){
            if(omega[t(i)]==1){
                if(gsl_rng_uniform(r)< lambda[1]) m(i,k) = 1;
            } else {
                if(gsl_rng_uniform(r)<lambda[0]) m(i,k) = 1;
            }
        }
    }
}

// Here I am adding more variability among size, activity and coherence of the ensembles
void gen_sample5(gsl_rng* r, arma::imat &m, arma::uvec &t, int P){
    m.fill(0);
    int i,k;
	int precision=200;
    int ncells=m.n_rows;
    int ntimes=m.n_cols;
    arma::mat lambda(P,2);
    arma::vec n(P);
    arma::vec sumn(P);
    arma::vec alphaN(P); alphaN.fill(2);
    n = dirichlet(r,alphaN);
    arma::vec p(P);
    for(i=0; i<P; ++i) {
        p(i)=gsl_ran_beta(r,0.01*precision,0.99*precision);
        lambda(i,0)=gsl_ran_beta(r,0.01*precision,0.99*precision);
        lambda(i,1)=gsl_ran_beta(r,0.8*precision,0.2*precision);
    }

    sumn=arma::cumsum(n);

    int omega[P];

    int label=0;
    for(i=0;i<ncells;++i){
        t(i) = label;
        if(i>=sumn(label)*ncells && label<P-1) label++;
    }

    for(k=0;k<ntimes;++k){
        for(i=0;i<P;++i) omega[i] = (gsl_rng_uniform(r) < p[i]) ? 1 : 0;
        for(i=0;i<ncells;++i){
            if(omega[t(i)]==1){
                if(gsl_rng_uniform(r)< lambda(t(i),1)) m(i,k) = 1;
            } else {
                if(gsl_rng_uniform(r)<lambda(t(i),0)) m(i,k) = 1;
            }
        }
    }
}

// Here I am adding more variability among size, activity and coherence of the ensembles but fixing the averages
void gen_sample6(gsl_rng* r, arma::imat &m, arma::uvec &t, int P,double lambda0,double lambda1, double pmu, bool lastIsFree){
    m.fill(0);
    int i,k;
	int precision=200;
    int ncells=m.n_rows;
    int ntimes=m.n_cols;
    arma::mat lambda(P,2);
    arma::vec n(P);
    arma::vec sumn(P);
    arma::vec alphaN(P); alphaN.fill(2);
    n = dirichlet(r,alphaN);
    arma::vec p(P);
    for(i=0; i<P; ++i) {
        p(i)=gsl_ran_beta(r,pmu*precision,(1-pmu)*precision);
        lambda(i,0)=gsl_ran_beta(r,lambda0*precision,(1-lambda0)*precision);
        lambda(i,1)=gsl_ran_beta(r,lambda1*precision,(1-lambda1)*precision);
    }

    if(lastIsFree) lambda(P-1,1)=lambda(P-1,0);

    sumn=arma::cumsum(n);

    int omega[P];

    int label=0;
    for(i=0;i<ncells;++i){
        t(i) = label;
        if(i>=sumn(label)*ncells && label<P-1) label++;
    }

    for(k=0;k<ntimes;++k){
        for(i=0;i<P;++i) omega[i] = (gsl_rng_uniform(r) < p[i]) ? 1 : 0;
        for(i=0;i<ncells;++i){
            if(omega[t(i)]==1){
                if(gsl_rng_uniform(r)< lambda(t(i),1)) m(i,k) = 1;
            } else {
                if(gsl_rng_uniform(r)<lambda(t(i),0)) m(i,k) = 1;
            }
        }
    }
}

// same as gen_sample6 but with multimembership
void gen_sample6_multimem(gsl_rng* r, arma::imat &m, arma::umat &t, int P,double lambda0,double lambda1, double pmu, bool lastIsFree){
    m.fill(0);
    int i,j,k;
    int precision=40;
    int ncells=m.n_rows;
    int ntimes=m.n_cols;
    arma::mat lambda(P,2);
    arma::vec n(P);
    arma::uvec multiplicity(ncells); multiplicity.fill(1);
    arma::vec sumn(P);
    arma::vec alphaN(P); alphaN.fill(2);
    n = dirichlet(r,alphaN);
    arma::vec p(P);
    for(i=0; i<P; ++i) {
        p(i)=gsl_ran_beta(r,pmu*precision,(1-pmu)*precision);
        lambda(i,0)=gsl_ran_beta(r,lambda0*precision,(1-lambda0)*precision);
        lambda(i,1)=gsl_ran_beta(r,lambda1*precision,(1-lambda1)*precision);
    }

    if(lastIsFree) lambda(P-1,1)=lambda(P-1,0);

    sumn=arma::cumsum(n);

    int omega[P];

    int label=0;
    for(i=0;i<ncells;++i){
        t(i,0) = label;
        t(i,1) = label;
        t(i,2) = label;
        if(gsl_rng_uniform(r)<.2){
            multiplicity(i)++;
            while(t(i,1)==label) t(i,1)=gsl_rng_uniform_int(r,P);

            if(gsl_rng_uniform(r)<.2){
                multiplicity(i)++;
                while(t(i,2) == label || t(i,2) == t(i,1)) t(i,2) = gsl_rng_uniform_int(r,P);
            }
        }

        if(i>=sumn(label)*ncells && label<P-1) label++;
    }

    for(k=0;k<ntimes;++k){
        for(i=0;i<P;++i) omega[i] = (gsl_rng_uniform(r) < p[i]) ? 1 : 0;
        for(i=0;i<ncells;++i){
            for(j=0;j<multiplicity(i);j++){
                if(omega[t(i,j)]==1){
                    if(gsl_rng_uniform(r) < lambda(t(i,j),1)) m(i,k) = 1;
                } else {
                    if(gsl_rng_uniform(r) < lambda(t(i,j),0)) m(i,k) = 1;
                }

                if(m(i,k)==1) break;

            }
        }
    }
}

void write_sample(ofstream &file, arma::mat m){
	for(int i = 0; i<m.n_rows;++i){
		for(int j =0 ;j<m.n_cols;++j){
			file<<m(i,j)<< ' ';
		}
		file<<' '<<endl;
	}
}

bool match_prediction(const arma::uvec &u, int* pred, int N, int P,int* res,double* acc){
    unsigned int i,j;
	int occu[P];
    int orig[P];

	for(i=0;i<P;i++){
		orig[i]=0;
		for(j=0;j<P;j++) occu[j]=0;
		for(j=0;j<N;j++){
			if(u(j)==i){
				occu[pred[j]]++;
				orig[i]++;
			}
		}
		int max_match=0;
		int nmatch=0;
		for(j=0;j<P;j++){
			if(occu[j]>nmatch){
				nmatch=occu[j];
				max_match=j;
			}
		}

		res[i]=max_match;
        acc[i]=nmatch*1.0/orig[i];
	}


    for(i=0;i<P;i++){
        unsigned int counter=0;
        for(j=0;j<P;j++) if(res[i]==res[j]) counter++;
        if(counter>1) return(false);
    }

    return(true);

}

bool match_prediction(const arma::uvec &u, int* pred, int N, int P,int* res,double &acc){
    unsigned int i,j;
	int occu[P];
    int orig[P];

    acc=1;
	for(i=0;i<P;i++){
		orig[i]=0;
		for(j=0;j<P;j++) occu[j]=0;
		for(j=0;j<N;j++){
			if(u(j)==i){
				occu[pred[j]]++;
				orig[i]++;
			}
		}
		int max_match=0;
		int nmatch=0;
		for(j=0;j<P;j++){
			if(occu[j]>nmatch){
				nmatch=occu[j];
				max_match=j;
			}
		}

		res[i]=max_match;
        acc=min(acc,nmatch*1.0/orig[i]);
	}


    for(i=0;i<P;i++){
        unsigned int counter=0;
        for(j=0;j<P;j++) if(res[i]==res[j]) counter++;
        if(counter>1) return(false);
    }

    return(true);

}

// This function returns TRUE if different original labels are mapped into different predicted membership
// Otherwise FALSE if two or more original groups are labeled with the same predicted membership.
// NOTE: P is the final membership generated by the sampler, not the original one.
// The original number of ensembles is given as the maximum value of the vector u + 1
bool match_prediction2(const arma::uvec &u, int* pred, int N, int P,int* res,double &acc){

    unsigned int i,j;
    int Ptrue=(int)(u.max())+1;

    int* occu = new int[P]; // Vector of occurrencies of predicted membership i from 1:Ptrue
    int* orig = new int[Ptrue]; // Vector of original group sizes

    acc=1;
    for(i=0;i<Ptrue;i++){
        orig[i]=0;
        // Initialize the number of occurrencies
        for(j=0;j<P;j++) occu[j]=0;

        // search all neurons with original label = i
        for(j=0;j<N;j++){
            if(u(j)==i){
                occu[pred[j]]++;
                orig[i]++;
            }
        }

        // Calculate the predicted membership which corresponds to the original one
        int max_match=0;
        int nmatch=0;
        for(j=0;j<P;j++){
            if(occu[j]>nmatch){
                nmatch=occu[j];
                max_match=j;
            }
        }

        // Assign the matching predicted membership to the vector res
        res[i]=max_match;

        // Set the accuracy of the prediction to the minimum value of the ratio between
        // right assignments and total neurons in the corresponding original group.
        acc=min(acc,nmatch*1.0/orig[i]);

    }


    delete [] occu;    
    delete [] orig;


    // Check for groups assigned to the same predicted membership
    for(i=0;i<Ptrue;i++){
        unsigned int counter=0;
        for(j=0;j<Ptrue;j++) if(res[i]==res[j]) counter++;
        if(counter>1) return(false);
    }

    return(true);

}

// GIBBS
arma::vec dirichlet(gsl_rng *r, const arma::vec &alpha){
    arma::vec theta(alpha.n_elem);
    gsl_ran_dirichlet(r,alpha.n_elem,alpha.memptr(),theta.memptr());
    return(theta);
}

bool fileExists(const char* file) {
    struct stat buf;
    return (stat(file, &buf) == 0);
}

void expand_binary(arma::fmat &s,int size){
    unsigned int i,l,k;
    arma::uword min_ind, max_ind;
    arma::uvec selec;

    for(i=0;i<s.n_rows;++i){

        selec=find(s.row(i));

        for(k=0;k<selec.n_elem;++k){
            min_ind= selec(k) >= size  ? selec(k)-size
                                         : 0;
            max_ind= selec(k)+size < s.n_cols ? selec(k)+size
                                              : s.n_cols-1;

            for(l=min_ind;l<=max_ind;l++) s(i,l)=1;
        }

    }

}

void filter_binary(const arma::fmat &in,arma::fmat &out, arma::uvec &cell_selection, int threshold, int threshold2){
    arma::frowvec tmp;
    arma::uvec selec;
    tmp=arma::sum(in);
    selec=find(tmp>threshold);
    out=in.cols(selec);
    selec=find(arma::sum(out,1)>threshold2);
    out=out.rows(selec);
    cell_selection=selec;
}

void filter_binary(const arma::imat &in,arma::imat &out, arma::uvec &cell_selection, int threshold){
    arma::irowvec tmp;
    arma::uvec selec;
    tmp=arma::sum(in);
    selec=find(tmp>threshold);
    out=in.cols(selec);
    selec=find(arma::sum(out,1));
    out=out.rows(selec);
    cell_selection=selec;
}

void filter_binary2(const arma::imat &in,arma::imat &out, arma::uvec &cell_selection, int threshold, int threshold2){
    arma::irowvec tmp;
    arma::uvec selec_row, selec_col, selec;
    selec_row=find(arma::sum(in,1)>threshold2);
	out=in.rows(selec_row);

    tmp=arma::sum(out);
    selec_col=find(tmp>threshold);
    out=out.cols(selec_col);

    selec=find(arma::sum(out,1));
    out=out.rows(selec);
    cell_selection=selec_row.rows(selec);
}

void SVDrot(arma::mat &input, arma::uvec &cell_labels, unsigned int nguess)
{
    // Extraction of synchronous frames based on the threshold value
    unsigned int i,j,k;

    cell_labels=arma::randi<arma::uvec>(input.n_rows,arma::distr_param(0,nguess-1));

    arma::mat U,V;
    arma::mat Urot;
    arma::vec s;

    arma::vec activity=arma::sum(input,1);
    arma::uvec selec=arma::find(activity);

    arma::svd(U,s,V,input.rows(selec));

    // generate Varimax instance
    arma::mat Ucut=U.cols(0,nguess-1);

    Vari rot(Ucut);
    rot.varimax(1e-5,true);
    Urot = Ucut * rot.rot;

    // Note: the function .max(uword &idx_max) is not mentioned in the armadillo documentation but
    // people use it.
    for(j=0;j<selec.n_elem;++j)
        Urot.row(j).max(cell_labels(selec(j)));


}
			
double gsl_ran_beta_pdflog(double x, double alpha, double beta){
				if(x==0 || x==1){
					std::cout<<"beta out of range"<<endl;
        std::abort();
    } else {
        return( (alpha-1)*log(x) + (beta-1)*log(1-x)) - gsl_sf_lnbeta(alpha,beta);
    }
}
