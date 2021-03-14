#pragma once

#include "Sampler.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace std;

class SigTracer{
    // arguments
    string data_directory; // Experiment name.
    string sample_name; // sample name.
    int iteration; // Replicate index.
    int MaxIter; // # max loop for CVB;
    string ABS_path;

    // Constants
    int K; // # signatures.
    int N; // # mutations for each sample.
    int V; // # mutation types.

    vector<int> MC; // ith mutation type.
    vector<int> D; // total reads.
    vector<int> B; // variant reads.
    vector<int> n_ref; // reference reads.
    vector<int> CN_normal; // copy number of normal cell.
    vector<int> CN_major; // copy number of major allele in cancer cell.
    vector<int> CN_tumor; // copy number of cancer cell.
    double purity; // sample purity.

    string inference_mode; // ms = model_selection or pe = parameter_estimation
    bool with_noise;
    string pe_state;
    string theta_state;

    // initialize
    vector<vector<vector<vector<double> > > > qzUC;
        // responsibility

    int J; // # clones.
    vector<double> pi; // fraction of jth clone.
    vector<double> alpha_pi; // hyper-parameter for clone.
    
    vector<double> lambda; // expected CCF for jth clone.
    vector<double> rho; // hyper-parameter for overdispersion.
    double xi;
    vector<vector<double> > qU;
        // expected responsibility with jth clone in mutation[i].

    vector<vector<double> > theta; // activity for jth clone and kth signature.
    vector<vector<double> > alpha_theta; // hyper-parameter for activity.

    vector<vector<double> > phi; // Mutational distributions of reference signatures.
    vector<vector<double> > reference_phi;
    vector<string> ref_names; // Reference signature names, e.g., SBS1.
    vector<string> ref_sig; // Used reference signature.

    vector<vector<double> > qz;
        // expected responsibility with kth mutation signature in mutation[i].
    vector<vector<vector<double> > > qzU;
    vector<vector<double> > qC;
        // expected responsibility with cth copy number in mutation[i].
    vector<vector<vector<double> > > qUC;

    //expected_values
    vector<vector<double> > Njk;
    vector<double> Nj;

    // run
    double temp_vlb, old_vlb; // Approximated Variational Lower Bound.

public:
    SigTracer(string x, string y, int j, int z, vector<string> r, string path, string m, bool w);
    bool exit_call;
    
    void run();
    void initialize();
    void initialize_clone();
    void initialize_copy_number();
    void initialize_activity();
    void initialize_signature();
    void load_reference();
    void initialize_expected_value();

    void Update_qzUC();
    void Update_clone();
    void Update_activity();
    void Update_BB_parameters();
    void Update_alpha();
    void calc_vlb();

    void load_data();
    void write_result();
    
    void Normalize(vector<double> &vec);
    void Normalize_tensor(vector<vector<vector<double> > > &vec);
    double my_digamma(double def);
    double my_lgamma(double def);
};


SigTracer::SigTracer(string x, string y, int j, int z, vector<string> r, string path, string m, bool w){
    data_directory = x;
    sample_name = y;
    J = j;
    iteration = z;
    for (int i=0; i<r.size(); i++){
        ref_sig.push_back(r[i]);
    }
    ABS_path = path;
    inference_mode = m;
    with_noise = w;
    
    exit_call = false;
    MaxIter = 1000;
}


void SigTracer::run(){
    if(inference_mode == "ms"){
        pe_state = "vb"; theta_state = "sym";
    }
    else{
        if(with_noise){
            pe_state = "cvb"; theta_state = "asy";
        }
        else{
            pe_state = "cvb"; theta_state = "sym";
        }
    }

    initialize(); old_vlb = -1e100;
    for (int iter=0; iter<MaxIter; iter++){
        Update_qzUC();
        Update_clone();
        Update_activity();
        calc_vlb();
        if(isfinite(temp_vlb) == false) exit_call = true;
        if(((fabs((temp_vlb-old_vlb)/old_vlb) < 1e-10) && (iter >= 500))|| exit_call == true) break;
        else old_vlb = temp_vlb;
        Update_BB_parameters();
        Update_alpha();
    }
}


void SigTracer::initialize(){
    initialize_signature();

    qzUC.resize(N);
    for (int i=0; i<N; i++){
        qzUC[i].resize(K);
        for (int k=0; k<K; k++){
            qzUC[i][k].resize(J);
            for (int j=0; j<J; j++){
                qzUC[i][k][j].resize(CN_major[i], -1.0);
                for (int c=0; c<CN_major[i]; c++){
                    qzUC[i][k][j][c] = Sampler::uniform(1.0, 5.0);
                }
            }
        }
        Normalize_tensor(qzUC[i]);
    }

    initialize_clone();
    initialize_copy_number();
    initialize_activity();
    initialize_expected_value();
    Update_clone();
    Update_activity();
    Update_BB_parameters();
    Update_alpha();
}


void SigTracer::initialize_clone(){
    alpha_pi.resize(J, 1.0);
    if(inference_mode == "pe"){
        for (int j=0; j<J; j++) alpha_pi[j] = pow(10.0, Sampler::uniform(-1.0, 0.0));
    }
    
    pi.resize(J, 1.0);
    
    qU.resize(N);
    for (int i=0; i<N; i++){
        qU[i].resize(J, 0.0);
        for (int j=0; j<J; j++){
            for (int k=0; k<K; k++){
                for (int c=0; c<CN_major[i]; c++){
                    qU[i][j] += qzUC[i][k][j][c];
                }
            }
        }
    }
    qzU.resize(N);
    for (int i=0; i<N; i++){
        qzU[i].resize(K);
        for (int k=0; k<K; k++){
            qzU[i][k].resize(J, 0.0);
            for (int j=0; j<J; j++){
                for (int c=0; c<CN_major[i]; c++){
                    qzU[i][k][j] += qzUC[i][k][j][c];
                }
            }
        }
    }
    
    lambda.resize(J, -1);
    for (int j=0; j<J; j++){
        lambda[j] = Sampler::uniform(0.001, 0.999);
    }
    rho.resize(J, -1);
    for (int j=0; j<J; j++) rho[j] = pow(10.0, Sampler::uniform(1.0, 2.0));
    xi = pow(10.0, Sampler::uniform(-1.0,1.0));
}


void SigTracer::initialize_copy_number(){
    qC.resize(N);
    for (int i=0; i<N; i++){
        qC[i].resize(CN_major[i], 0.0);
        for (int c=0; c<CN_major[i]; c++){
            for (int k=0; k<K; k++){
                for (int j=0; j<J; j++){
                    qC[i][c] += qzUC[i][k][j][c];
                }
            }
        }
    }
    qUC.resize(N);
    for (int i=0; i<N; i++){
        qUC[i].resize(J);
        for (int j=0; j<J; j++){
            qUC[i][j].resize(CN_major[i], 0.0);
            for (int c=0; c<CN_major[i]; c++){
                for (int k=0; k<K; k++){
                    qUC[i][j][c] += qzUC[i][k][j][c];
                }
            }
        }
    }

    D.resize(N, 0.0);
    for (int i=0; i<N; i++){
        D[i] = B[i] + n_ref[i];
    }
}


void SigTracer::initialize_activity(){
    alpha_theta.resize(J);
    double temp_alpha = pow(10.0, Sampler::uniform(-1.0, 0.0));
    for (int j=0; j<J; j++){
        alpha_theta[j].resize(K, 1.0);
        for (int k=0; k<K; k++){
            if(theta_state == "sym") alpha_theta[j][k] = temp_alpha;
            else if(theta_state == "asy") alpha_theta[j][k] = pow(10.0, Sampler::uniform(-1.0, 0.0));
        }
    }

   
    theta.resize(J);
    for (int j=0; j<J; j++){
        theta[j].resize(K, -1);
        for (int k=0; k<K; k++){
            theta[j][k] = Sampler::uniform(1.0, 5.0);
        }
    }

    qz.resize(N);
    for (int i=0; i<N; i++){
        qz[i].resize(K, 0.0);
        for (int k=0; k<K; k++){
            for(int j=0; j<J; j++){
                for (int c=0; c<CN_major[i]; c++){
                    qz[i][k] += qzUC[i][k][j][c];
                }
            }
        }
    }
}


void SigTracer::initialize_signature(){ 
    load_reference();
    bool all_flag = false;
    for (int k=0; k<ref_sig.size(); k++){
        if(ref_sig[k] == "all" || ref_sig[k] == "All") all_flag = true;
    }
    if(all_flag == true) K = ref_names.size();
    else K = ref_sig.size();
    
    phi.resize(K);
    for (int k=0; k<K; k++) phi[k].resize(V, -1);
    
    for (int k=0; k<K; k++){
        if(all_flag == true){
            for (int v=0; v<V; v++){
                phi[k][v] = reference_phi[k][v];
            }
        }
        else{
            string temp_ref_sig = ref_sig[k];
            int target_k = -1;
            for (int i=0; i<ref_names.size(); i++){
                if(temp_ref_sig == ref_names[i]){
                    target_k = i;
                    break;
                }
            }
            if(target_k == -1){
                cout << "Cannot find " + temp_ref_sig + " in reference file." << endl;
                exit(1);
            }
            for (int v=0; v<V; v++){
                phi[k][v] = reference_phi[target_k][v];
            }
        }
    }
    for (int k=0; k<K; k++){
        for (int v=0; v<V; v++){
            if(phi[k][v] < 1e-100) phi[k][v] = 1e-100;
        }
        Normalize(phi[k]);
    }
}


void SigTracer::initialize_expected_value(){
    Njk.resize(J);
    for (int j=0; j<J; j++){
        Njk[j].resize(K, 0.0);
        for (int k=0; k<K; k++){
            for (int i=0; i<N; i++){
                for (int c=0; c<CN_major[i]; c++){
                    Njk[j][k] += qzUC[i][k][j][c];
                }
            }
        }
    }
    Nj.resize(J, 0.0);
    for (int j=0; j<J; j++){
        for (int i=0; i<N; i++){
            for (int k=0; k<K; k++){
                for (int c=0; c<CN_major[i]; c++){
                    Nj[j] += qzUC[i][k][j][c];
                }
            }
        }
    }
}


void SigTracer::Update_qzUC(){
    for (int i=0; i<N; i++){
        vector<vector<vector<double> > > old_qzUC; old_qzUC.resize(K);
        for (int k=0; k<K; k++){
            old_qzUC[k].resize(J);
            for (int j=0; j<J; j++){
                old_qzUC[k][j].resize(CN_major[i], -1.0);
                for (int c=0; c<CN_major[i]; c++){
                    old_qzUC[k][j][c] = qzUC[i][k][j][c];
                }
            }
        }
        double binom = my_lgamma(D[i]+1) - my_lgamma(B[i]+1) - my_lgamma(D[i]-B[i]+1);
        for (int k=0; k<K; k++){
            for (int j=0; j<J; j++){
                if(pe_state == "vb"){
                    double sum_theta = 0.0;
                    for (int kk=0; kk<K; kk++) sum_theta += theta[j][kk];
                    for (int c=0; c<CN_major[i]; c++){
                        qzUC[i][k][j][c] = phi[k][MC[i]];
                        double eta = (purity * (c+1)) / (purity * CN_tumor[i] + (1-purity) * CN_normal[i]);
                        double AA = rho[j] * lambda[j] * eta;
                        double BB = rho[j] * (1-lambda[j] * eta);
                        qzUC[i][k][j][c] *= exp(binom + my_lgamma(rho[j]) - my_lgamma(D[i] + rho[j]) + my_lgamma(B[i]+AA) 
                                         + my_lgamma(D[i]-B[i]+BB) - my_lgamma(AA) - my_lgamma(BB));
                        double temp = my_digamma(theta[j][k]) + my_digamma(pi[j]) - my_digamma(sum_theta);
                        qzUC[i][k][j][c] *= exp(temp);
                    }
                }
                else if(pe_state == "cvb"){
                    double temp_Njk = Njk[j][k] - qzU[i][k][j];
                    double temp_Nj = Nj[j] - qU[i][j];
                    double sum_alpha = 0.0;
                    for (int kk=0; kk<K; kk++) sum_alpha += alpha_theta[j][kk];
                    for (int c=0; c<CN_major[i]; c++){
                        qzUC[i][k][j][c] = phi[k][MC[i]];
                        double eta = (purity * (c+1)) / (purity * CN_tumor[i] + (1-purity) * CN_normal[i]);
                        double AA = rho[j] * lambda[j] * eta;
                        double BB = rho[j] * (1-lambda[j] * eta);
                        qzUC[i][k][j][c] *= exp(binom + my_lgamma(rho[j]) - my_lgamma(D[i] + rho[j]) + my_lgamma(B[i]+AA) 
                                         + my_lgamma(D[i]-B[i]+BB) - my_lgamma(AA) - my_lgamma(BB));
                        qzUC[i][k][j][c] *= (temp_Njk + alpha_theta[j][k]);
                        qzUC[i][k][j][c] /= (temp_Nj + sum_alpha);
                        qzUC[i][k][j][c] *= (temp_Nj + alpha_pi[j]);
                    }
                }
            }
        }
        Normalize_tensor(qzUC[i]);
        for (int k=0; k<K; k++){
            for (int j=0; j<J; j++){
                for (int c=0; c<CN_major[i]; c++){
                    if(qzUC[i][k][j][c] < 1e-100) qzUC[i][k][j][c] = 1e-100;
                }
            }
        }
        Normalize_tensor(qzUC[i]);
        for (int k=0; k<K; k++){
            for (int j=0; j<J; j++){
                for (int c=0; c<CN_major[i]; c++){
                    double diff = qzUC[i][k][j][c] - old_qzUC[k][j][c];
                    Njk[j][k] += diff;
                    Nj[j] += diff; qz[i][k] += diff;
                    qU[i][j] += diff; qC[i][c] += diff;
                    qUC[i][j][c] += diff; qzU[i][k][j] += diff;
                }
            }
        }
    }
}


void SigTracer::Update_clone(){
    for(int j=0; j<J; j++){
        pi[j] = Nj[j] + alpha_pi[j];
    }
}


void SigTracer::Update_BB_parameters(){
    // Update dispersion parameter.
    for (int j=0; j<J; j++){
        for (int count=0; count<10; count++){
            double old_rho = rho[j];
            double numerator = 0.0; double denominator = xi;
            for (int i=0; i<N; i++){
                denominator += qU[i][j] * (my_digamma(D[i] + old_rho) - my_digamma(old_rho));
                for (int c=0; c<CN_major[i]; c++){
                    double eta = (purity * (c+1))/(purity * CN_tumor[i] + (1-purity) * CN_normal[i]);
                    double AA = lambda[j] * eta * (my_digamma(B[i] + old_rho * lambda[j] * eta) - my_digamma(old_rho * lambda[j] * eta));
                    double BB = (1 - lambda[j] * eta) * (my_digamma(D[i]-B[i]+old_rho * (1-lambda[j] * eta)) - my_digamma(old_rho * (1-lambda[j] * eta)));
                    numerator += qUC[i][j][c] * (AA + BB);
                }
            }
            rho[j] *= numerator/denominator;
            if((fabs((rho[j] - old_rho)/old_rho) < 1e-3) || rho[j] < min(1.0/lambda[j], 1.0/(1.0-lambda[j]))){
                if(rho[j] < min(1.0/lambda[j], 1.0/(1.0-lambda[j]))) rho[j] = min(1.0/lambda[j], 1.0/(1.0-lambda[j]));
                break;
            }
        }
    }
    double denominator = 0.0;
    for (int j=0; j<J; j++) denominator += rho[j];
    xi = (double)J/denominator;
    
    // Update CCF, lambda[j];
    for (int j=0; j<J; j++){
        for (int count=0; count<10; count++){
            double old_lambda = lambda[j];
            double numerator = 0.0; double denominator = 0.0;
            for (int i=0; i<N; i++){
                for (int c=0; c<CN_major[i]; c++){
                    double eta = (purity * (c+1))/(purity * CN_tumor[i] + (1-purity) * CN_normal[i]);
                    double AA = rho[j] * old_lambda * eta * (my_digamma(B[i]+rho[j]*old_lambda*eta) - my_digamma(rho[j]*old_lambda*eta));
                    double BB = rho[j] * (1-old_lambda*eta) * (my_digamma(D[i]-B[i]+rho[j]*(1-old_lambda*eta)) - my_digamma(rho[j]*(1-old_lambda*eta)));
                    numerator += qUC[i][j][c] * AA;
                    denominator += qUC[i][j][c] * eta * (AA + BB);
                }
            }
            lambda[j] = numerator/denominator;
            if(fabs(lambda[j]-old_lambda)/old_lambda < 1e-3 || lambda[j] > 1.0){
                if(lambda[j] > 1.0) lambda[j] = 1 - 1e-5;
                break;
            }
        }
    }
}


void SigTracer::Update_activity(){
    for (int j=0; j<J; j++){
        for (int k=0; k<K; k++){
            theta[j][k] = Njk[j][k] + alpha_theta[j][k];
        }
    }
}


void SigTracer::Update_alpha(){
    // Update_alpha_theta
    // If all parameters take the same values;
    if(theta_state == "sym"){
        for (int count=0; count<10; count++){
            double old_alpha = alpha_theta[0][0];
            double numerator = 0.0; double denominator = 0.0;
            for (int j=0; j<J; j++){
                for (int k=0; k<K; k++){
                    numerator += my_digamma(Njk[j][k] + old_alpha) - my_digamma(old_alpha);
                }
                denominator += my_digamma(Nj[j] + K * old_alpha) - my_digamma(K * old_alpha);
            }
            denominator *= K;
            if(numerator != 0.0 && denominator != 0.0){
                double temp_alpha = old_alpha * numerator / denominator;
                for (int j=0; j<J; j++){
                    for (int k=0; k<K; k++) alpha_theta[j][k] = temp_alpha;
                }
            }
            if(fabs(alpha_theta[0][0] - old_alpha)/old_alpha < 1e-3) break;
        }
    }
    else if(theta_state == "asy"){ // If all parameters takes the different values;
        for (int j=0; j<J; j++){
            for (int k=0; k<K; k++){
                for (int count=0; count<10; count++){
                    double old_alpha = alpha_theta[j][k];
                    double sum_alpha = 0.0;
                    for(int kk=0; kk<K; kk++) sum_alpha += alpha_theta[j][kk];
                    double numerator = my_digamma(Njk[j][k] + old_alpha) - my_digamma(old_alpha);
                    double denominator = my_digamma(Nj[j] + sum_alpha) - my_digamma(sum_alpha);
                    if(numerator != 0.0 && denominator != 0.0) alpha_theta[j][k] *= numerator/denominator;
                    if(fabs(alpha_theta[j][k] - old_alpha)/old_alpha < 1e-3) break;
                }
            }
        }
    }
    
  
    // Update_alpha_pi
    for (int j=0; j<J; j++){
        for (int count=0; count<10; count++){
            double old_alpha = alpha_pi[j];
            double sum_alpha = 0.0;
            for (int jj=0; jj<J; jj++) sum_alpha += alpha_pi[jj];
            double numerator = my_digamma(Nj[j] + old_alpha) - my_digamma(old_alpha);
            double denominator = my_digamma(N + sum_alpha) - my_digamma(sum_alpha);
            if(numerator != 0.0 && denominator != 0.0) alpha_pi[j] *= numerator/denominator;
            if(fabs(alpha_pi[j]-old_alpha)/old_alpha < 1e-3) break;
        }
    }
}


void SigTracer::calc_vlb(){
    temp_vlb = 0.0;
    for (int i=0; i<N; i++){
        for (int k=0; k<K; k++){
            temp_vlb += qz[i][k] * log(phi[k][MC[i]]);
            for (int j=0; j<J; j++){
                for (int c=0; c<CN_major[i]; c++){
                    temp_vlb -= qzUC[i][k][j][c] * log(qzUC[i][k][j][c]);
                }
            }
        }
        double binom = my_lgamma(D[i]+1) - my_lgamma(B[i]+1) - my_lgamma(D[i]-B[i]+1);
        for (int c=0; c<CN_major[i]; c++){
            double eta = (purity * (c+1)) / (purity * CN_tumor[i] + (1-purity) * CN_normal[i]);
            for (int j=0; j<J; j++) {
                double AA = rho[j] * lambda[j] * eta;
                double BB = rho[j] * (1-lambda[j] * eta);
                temp_vlb += qUC[i][j][c] * (binom + my_lgamma(rho[j]) - my_lgamma(D[i] + rho[j]) 
                         + my_lgamma(B[i]+AA) + my_lgamma(D[i]-B[i]+BB) - my_lgamma(AA) - my_lgamma(BB));
            }
        }
    }
    double sum_alpha_pi = 0.0; double sum_pi = 0.0;
    for (int j=0; j<J; j++){
        sum_alpha_pi += alpha_pi[j]; sum_pi += pi[j];
    }
    temp_vlb += my_lgamma(sum_alpha_pi) - my_lgamma(sum_pi);
    for (int j=0; j<J; j++){
        temp_vlb -= my_lgamma(alpha_pi[j]) - my_lgamma(pi[j]);
    }
    for (int j=0; j<J; j++){
        double sum_alpha_theta = 0.0; double sum_theta = 0.0;
        for (int k=0; k<K; k++){
            sum_alpha_theta += alpha_theta[j][k];
            sum_theta += theta[j][k];
        }
        temp_vlb += my_lgamma(sum_alpha_theta) - my_lgamma(sum_theta);
        for (int k=0; k<K; k++){
            temp_vlb -= my_lgamma(alpha_theta[j][k]) - my_lgamma(theta[j][k]);
        }
    }
    for (int j=0; j<J; j++){
        temp_vlb += log(xi) - (xi * rho[j]);
    }
}


void SigTracer::load_data(){
    ifstream ifs;
    string file_name = ABS_path + "data/" + data_directory;
    file_name += "/" + sample_name + ".csv";
    ifs.open(file_name.c_str(), ios::in);
    if(!ifs){
        cout << "Cannot open " + file_name  << endl;
        exit(1);
    }
    
    static char buf[10000000];
    char *temp;

    MC.reserve(1000000);
    n_ref.reserve(1000000);
    B.reserve(1000000);
    CN_normal.reserve(1000000);
    CN_major.reserve(1000000);
    CN_tumor.reserve(1000000);

    ifs.getline(buf, 10000000);
    while(ifs.getline(buf, 10000000)){
        temp = strtok(buf, ","); temp = strtok(NULL, ",");
        temp = strtok(NULL, ","); temp = strtok(NULL, ",");
        n_ref.push_back(atoi(temp));
        temp = strtok(NULL, ","); B.push_back(atoi(temp));
        temp = strtok(NULL, ","); CN_normal.push_back(atoi(temp));
        temp = strtok(NULL, ","); temp = strtok(NULL, ",");
        temp = strtok(NULL, ","); CN_major.push_back(atoi(temp));
        temp = strtok(NULL, ","); CN_tumor.push_back(atoi(temp));
        temp = strtok(NULL, ","); MC.push_back(atoi(temp));
    }
    N = MC.size();
    ifs.close();
    
    file_name = ABS_path + "data/" + data_directory + "/purity.csv";
    ifs.open(file_name.c_str(), ios::in);
    if(!ifs){
        cout << "Cannot open " + file_name << endl;
        exit(1);
    }
    ifs.getline(buf, 10000000);
    while(ifs.getline(buf, 10000000)){
        temp = strtok(buf, ",");
        if(temp == sample_name){
            temp = strtok(NULL, ",");
            purity = atof(temp);
        }
    }
    ifs.close();
}


void SigTracer::load_reference(){
    ifstream ifs;
    string file_name = ABS_path + "data/ref/signature_probability_v3.1.csv";
    ifs.open(file_name.c_str(), ios::in);
    if(!ifs){
        cout << "Cannot open " + file_name << endl;
        exit(1);
    }
    char buf[1000000]; char *temp;
    ifs.getline(buf, 1000000); 
    temp = strtok(buf, ","); temp = strtok(NULL, ",");
    while(true){
        temp = strtok(NULL, ",");
        if(temp == NULL) break;
        ref_names.push_back(temp);
    }
    ref_names[ref_names.size()-1].pop_back();
    int ref_len = ref_names.size();
    reference_phi.resize(ref_len);
    for (int k=0; k<ref_len; k++) reference_phi[k].reserve(2000);
    while(ifs.getline(buf, 1000000)){
        temp = strtok(buf, ","); temp = strtok(NULL, ",");
        for (int k=0; k<ref_len; k++){
            temp = strtok(NULL, ",");
            reference_phi[k].push_back(atof(temp));
        }
    }
    V = reference_phi[0].size();
}


void SigTracer::write_result(){
    ofstream ofs;
    string file_name = ABS_path + "result/" + data_directory;
    file_name += "/" + sample_name + "/J" + to_string(J) + "_run"
              + to_string(iteration) + "_" + inference_mode;
    
    // write ELBO
    ofs.open(file_name + "_elbo.txt", ios::out);
    ofs << to_string(temp_vlb) << "\n";
    ofs.close();
    
    // write activity
    ofs.open(file_name + "_activity.txt", ios::out);
    for (int j=0; j<J; j++){
        Normalize(theta[j]);
        for (int k=0; k<K; k++){
            ofs << to_string(theta[j][k]);
            if(k != K-1) ofs << " ";
        }
        ofs << "\n";
    }
    ofs.close();
    
    // write pi
    ofs.open(file_name + "_pi.txt", ios::out);
    Normalize(pi);
    for (int j=0; j<J; j++){
        ofs << to_string(pi[j]);
        if(j != J-1) ofs << " ";
    }
    ofs << "\n";
    ofs.close();
    
    // write_BB_parameters
    ofs.open(file_name + "_BB.txt", ios::out);
    
    for (int j=0; j<J; j++){
        ofs << to_string(rho[j]);
        if(j != J-1) ofs << " ";
    }
    ofs << endl;
    for (int j=0; j<J; j++){
        ofs << to_string(lambda[j]);
        if(j != J-1) ofs << " ";
    }
    ofs << endl;
    ofs.close();

    // write signatures
    ofs.open(file_name + "_signature.txt", ios::out);
    for (int k=0; k<K; k++){
        for (int v=0; v<V; v++){
            ofs << to_string(phi[k][v]);
            if(v != V-1) ofs << " ";
        }
        ofs << "\n";
    }
    ofs.close();
    
    // write alpha
    ofs.open(file_name + "_alpha.txt", ios::out);
    for (int j=0; j<J; j++){
        for (int k=0; k<K; k++){
            ofs << to_string(alpha_theta[j][k]);
            if(k != K-1) ofs << " ";
        }
        ofs << "\n";
    }
    for (int j=0; j<J; j++){
        ofs << to_string(alpha_pi[j]);
        if(j != J-1) ofs << " ";
    }
    ofs << "\n";
    ofs.close();

    // write_qU
    ofs.open(file_name + "_qU.txt", ios::out);
    for (int i=0; i<N; i++){
        for (int j=0; j<J; j++){
            ofs << to_string(qU[i][j]);
            if(j != J-1) ofs << " ";
        }
        ofs << "\n";
    }
    ofs.close();

    // write_qz
    ofs.open(file_name + "_qz.txt", ios::out);
    for (int i=0; i<N; i++){
        for (int k=0; k<K; k++){
            ofs << to_string(qz[i][k]);
            if(k != K-1) ofs << " ";
        }
        ofs << "\n";
    }
    ofs.close();

    // write_qC
    ofs.open(file_name + "_qC.txt", ios::out);
    for (int i=0; i<N; i++){
        for (int c=0; c<CN_major[i]; c++){
            ofs << to_string(qC[i][c]);
            if(c != CN_major[i]-1) ofs << " ";
        }
        ofs << "\n";
    }
    ofs.close();

    // write_qzU
    ofs.open(file_name + "_qzU.txt", ios::out);
    ofs << to_string(N) << " " << to_string(K) << " " << to_string(J) << "\n";
    for (int i=0; i<N; i++){
        for (int k=0; k<K; k++){
            for (int j=0; j<J; j++){
                ofs << to_string(qzU[i][k][j]);
                if(j != J-1) ofs << " ";
            }
            ofs << "\n";
        }
    }
    ofs.close();

    // write_qUC
    ofs.open(file_name + "_qUC.txt", ios::out);
    ofs << to_string(N) << " " << to_string(J) << "\n";
    for (int i=0; i<N; i++){
        ofs << to_string(CN_major[i]);
        if(i != N-1) ofs << " ";
    }
    ofs << "\n";
    for (int i=0; i<N; i++){
        for (int j=0; j<J; j++){
            for (int c=0; c<CN_major[i]; c++){
                ofs << to_string(qUC[i][j][c]);
                if(c != CN_major[i] - 1) ofs << " ";
            }
            ofs << "\n";
        }
    }
    ofs.close();
}


void SigTracer::Normalize(vector<double> &vec){
    int length = vec.size();
    double sum = 0.0;
    for (int i=0; i<length; i++) sum += vec[i];
    for (int i=0; i<length; i++) vec[i] /= sum;
}


void SigTracer::Normalize_tensor(vector<vector<vector<double> > > &vec){
    int dim_one = vec.size();
    double sum = 0.0;
    for (int i=0; i<dim_one; i++){
        int dim_two = vec[i].size();
        for (int j=0; j<dim_two; j++){
            int dim_three = vec[i][j].size();
            for (int k=0; k<dim_three; k++){
                sum += vec[i][j][k];
            }
        }
    }
    for (int i=0; i<dim_one; i++){
        int dim_two = vec[i].size();
        for (int j=0; j<dim_two; j++){
            int dim_three = vec[i][j].size();
            for (int k=0; k<dim_three; k++){
                vec[i][j][k] /= sum;
            }
        }
    }
}


double SigTracer::my_digamma(double def){
    double res;
    if(def < 1e-300){
        res = boost::math::digamma(1e-300);
    }
    else if(def > 1e300){
        res = boost::math::digamma(1e300);
    }
    else{
        res = boost::math::digamma(def);
    }
    return(res);
}


double SigTracer::my_lgamma(double def){
    double res;
    if(def < 1e-300){
        res = boost::math::lgamma(1e-300);
    }
    else{
        res = boost::math::lgamma(def);
    }
    return(res);
}


