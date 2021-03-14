import sys
import os
import utility as util
import random
import numpy as np
import math
from scipy.spatial.distance import cosine
import pickle

def main():
    purity = make_purity()
    true_CCF = make_true_CCF()
    
    CN_tumour_total = make_CN_tumour_total()
    CN_major, CN_minor = make_CN_tumour(CN_tumour_total)
    CN_mut = make_CN_mut(CN_major)
    
    n_total = make_n_total()
    rho = make_rho()
    pi = make_pi()
    U = make_U(pi)
    n_var, n_ref = make_n_var(true_CCF, n_total, U, rho, CN_tumour_total, CN_mut, purity)
    
    K, ref_names, mutational_distributions, activity = make_md()
    z = make_z(K, U, activity)
    M = make_M(mutational_distributions, z)

    sample_names = []
    for s in range(S):
        sample_names.append("sample-" + str(s+1))
        write_mutations(sample_names[s], CN_tumour_total[s], CN_major[s], CN_minor[s], CN_mut[s],
                        n_total[s], U[s], n_var[s], n_ref[s], z[s], M[s], ref_names[s])
    write_purity(sample_names, purity)
    write_ref_sig(sample_names, ref_names)
    save_ground_truth(sample_names, true_CCF, pi, mutational_distributions, activity, rho, U,
                      CN_mut, z, ref_names)


def make_purity():
    purity = np.zeros([S])
    for s in range(S):
        log_purity = np.random.uniform(-1.0, 0.0)
        purity[s] = math.exp(log_purity)
    return purity


def make_true_CCF():
    true_CCF = np.zeros([S, J])
    for s in range(S):
        while(True):
            for j in range(J):
                true_CCF[s,j] = np.random.uniform(0.0, 1.0)
            true_CCF[s] = np.sort(true_CCF[s])[::-1]
            flag = True
            for j in range(1, J):
                if(true_CCF[s,j-1]-true_CCF[s,j] < 0.10): flag = False 
            if(flag): break
    return true_CCF


def make_CN_tumour_total():
    CN_tumour_total = np.zeros([S,N])
    for s in range(S):
        for n in range(N):
            CN_tumour_total[s,n] = int(np.random.lognormal(1,0.3))
            if(CN_tumour_total[s,n] == 0): CN_tumour_total[s,n] = 1
    return CN_tumour_total


def make_CN_tumour(CN_tumour_total):
    CN_major = np.zeros([S,N]); CN_minor = np.zeros([S,N])
    for s in range(S):
        for n in range(N):
            temp_total = int(CN_tumour_total[s,n])
            temp = int(np.random.beta(5, 3) * temp_total)
            if(temp >= CN_tumour_total[s,n]/2):
                CN_major[s,n] = temp; CN_minor[s,n] = temp_total - temp
            else:
                CN_minor[s,n] = temp; CN_major[s,n] = temp_total - temp;
    return CN_major, CN_minor


def make_CN_mut(CN_major):
    CN_mut = np.zeros([S,N])
    for s in range(S):
        for n in range(N):
            temp_major = int(CN_major[s,n])
            prob_var = list(range(1, temp_major+1))
            prob = [1.0/temp_major for cn in range(temp_major)]
            CN_mut[s,n] = np.random.choice(prob_var, p=prob)
    return CN_mut


def make_n_total():
    n_total = np.zeros([S,N])
    prob_var = list(range(read_counts_min, read_counts_max+1))
    prob = [1.0/(read_counts_max-read_counts_min+1) for count in range(read_counts_max-read_counts_min+1)]
    for s in range(S):
        for n in range(N):
            n_total[s,n] = np.random.choice(prob_var, p=prob)
    return n_total


def make_rho():
    rho = np.zeros([S,J])
    for s in range(S):
        if(rho_state == "same"):
            temp = np.random.normal(60.0, 5.0)
            if (temp<=0.0): temp = 0.1
            for j in range(J): rho[s,j] = temp
        elif(rho_state == "different"):
            while(True):
                for j in range(J):
                    rho[s,j] = np.random.uniform(5.0, 100.0)
                temp = np.sort(rho[s])[::-1]
                flag = True
                for j in range(1, J):
                    if(temp[j-1]-temp[j] < 5.0): flag = False
                if(flag): break
    return rho


def make_pi():
    pi = np.zeros([S,J])
    alpha = [1.0 for j in range(J)]
    for s in range(S):
        while(True):
            pi_s = np.random.dirichlet(alpha)
            if(min(pi_s) > 0.05): break
        for j in range(J):
            pi[s,j] = pi_s[j]
    return pi


def make_U(pi):
    U = np.zeros([S,N])
    for s in range(S):
        for n in range(N):
            U[s,n] = np.random.choice(range(J), p=pi[s])
    return U


def make_n_var(true_CCF, n_total, U, rho, CN_tumour_total, CN_mut, purity):
    n_var = np.zeros([S,N]); n_ref = np.zeros([S,N])
    for s in range(S):
        temp_rho = rho[s]
        for n in range(N):
            temp_U = int(U[s,n])
            temp_total = n_total[s,n]
            temp_CCF = true_CCF[s,temp_U]
            eta = (purity[s] * CN_mut[s,n])/((purity[s] * CN_tumour_total[s,n]) + ((1-purity[s]) * CN_normal))
            temp_prob = np.random.beta(temp_rho[temp_U] * temp_CCF * eta, temp_rho[temp_U] * (1-temp_CCF * eta))
            temp_n_var = int(np.random.binomial(temp_total, temp_prob))
            if(temp_n_var == 0): n_var[s,n] = 1
            elif(temp_n_var == temp_total): n_var[s,n] = temp_total-1
            else: n_var[s,n] = temp_n_var
            n_ref[s,n] = temp_total - n_var[s,n]
    return n_var, n_ref


def make_md():
    K = [0 for s in range(S)]
    for s in range(S):
        K[s] = int(np.random.poisson(mean_K)+2)
    
    ref_names = [["" for k in range(K[s])] for s in range(S)]
    mutational_distributions = [np.zeros([K[s],V]) for s in range(S)]
    activity = [[[0.0 for k in range(K[s])] for j in range(J)] for s in range(S)]
    for s in range(S):
        alpha = [1.0 for k in range(K[s])]
        while(True):
            K_known = list(range(len(known_names)))
            sig_selected = random.sample(K_known, K[s])
            sig_selected = np.sort(sig_selected)
            for k,x in enumerate(sig_selected):
                ref_names[s][k] = known_names[x]
            for k,ref_name in enumerate(ref_names[s]):
                target_k = known_names.index(ref_name)
                for v in range(V):
                    mutational_distributions[s][k,v] = signatures[target_k,v]
            for k in range(K[s]):
                temp_sum = 0.0
                for v in range(V):
                    temp_sum += mutational_distributions[s][k,v]
                for v in range(V):
                    mutational_distributions[s][k,v] /= temp_sum
            for j in range(J):
                temp_activity = np.random.dirichlet(alpha)
                for k in range(K[s]):
                    activity[s][j][k] = temp_activity[k]
            flag = True
            mix = [[0.0 for v in range(V)] for j in range(J)]
            for j in range(J):
                for k in range(K[s]):
                    for v in range(V):
                        mix[j][v] += activity[s][j][k] * mutational_distributions[s][k][v]
            for j in range(J):
                for jj in range(j+1, J):
                    if(cosine(mix[j], mix[jj]) < 0.10 or 
                       cosine(activity[s][j], activity[s][jj]) < 0.10): flag = False
            if(flag == True): break
    return K, ref_names, mutational_distributions, activity


def make_z(K, U, activity):
    z = np.zeros([S,N])
    for s in range(S):
        for n in range(N):
            temp_U = int(U[s,n])
            temp_activity = activity[s][temp_U]
            z[s,n] = np.random.choice(range(K[s]), p=temp_activity)
    return z


def make_M(mutational_distributions, z):
    M = np.zeros([S,N])
    for s in range(S):
        for n in range(N):
            temp_z = int(z[s,n])
            temp_signature = mutational_distributions[s][temp_z]
            M[s,n] = np.random.choice(range(V), p=temp_signature)
    return M


def write_mutations(sample_name, CN_tumour_total, CN_major, CN_minor, CN_mut,
                    n_total, U, n_var, n_ref, z, M, ref_names):
    out = open("data/" + exp_name + "/" + sample_name + ".csv", "w")
    out.write("mutation_id,chromosome,position,ref_counts,var_counts,normal_cn,minor_cn,mut_cn,major_cn,total_cn,trinucleotide,signature,clone\n")
    for n in range(N):
        out.write("mut_" + str(n+1) + ",")
        out.write("1," + str(int(5*(n+1))) + ",")
        out.write(str(int(n_ref[n])) + ",")
        out.write(str(int(n_var[n])) + ",")
        out.write(str(int(CN_normal)) + ",")
        out.write(str(int(CN_minor[n])) + ",")
        out.write(str(int(CN_mut[n])) + ",")
        out.write(str(int(CN_major[n])) + ",")
        out.write(str(int(CN_tumour_total[n])) + ",")
        out.write(str(int(M[n])) + ",")
        out.write(ref_names[int(z[n])] + ",")
        out.write(str(int(U[n])) + "\n")
    out.close()


def write_purity(sample_names, purity):
    out = open("data/" + exp_name + "/purity.csv", "w")
    out.write("Sample_name,purity\n")
    for s in range(S):
        out.write(sample_names[s] + "," + str(purity[s]) + "\n")
    out.close()


def write_ref_sig(sample_names, ref_names):
    out = open("data/" + exp_name + "/ref_sig.csv", "w")
    for s, sample_name in enumerate(sample_names):
        out.write(sample_names[s])
        for k,x in enumerate(ref_names[s]):
            out.write("," + x)
        out.write("\n")
    out.close()


def save_ground_truth(sample_names, true_CCF, pi, mutational_distributions, activity, rho, U,
                      CN_mut, z, ref_names):
    ground_truth = {}
    ground_truth.update({"sample_names":sample_names})
    ground_truth.update({"CCF":true_CCF})
    ground_truth.update({"pi":pi})
    ground_truth.update({"signature":mutational_distributions})
    ground_truth.update({"activity":activity})
    ground_truth.update({"rho":rho})
    ground_truth.update({"U":U})
    ground_truth.update({"CN_mut":CN_mut})
    ground_truth.update({"z":z})
    ground_truth.update({"ref_sig":ref_names})
    file_name = "data/" + exp_name + "/ground_truth.pkl"
    out = open(file_name, "wb")
    pickle.dump(ground_truth, out)


def make_CCF(purity, CN_tumour_total, CN_major, n_var, n_ref):
    CCF = np.zeros([S,N])
    for s in range(S):
        temp_purity = purity[s]
        for n in range(N):
            temp_CN_total = int(CN_tumour_total[s,n])
            temp_CN_mut = np.random.choice(range(int(CN_major[s,n]))) + 1
            temp_n_var = int(n_var[s,n]); temp_n_ref = int(n_ref[s,n])
            temp_VAF = np.random.beta(temp_n_var, temp_n_ref)
            CCF[s,n] = (temp_VAF * (temp_purity * temp_CN_total + (1-temp_purity) * CN_normal)) / (temp_purity * temp_CN_mut)
    return CCF


def make_bootstrap_samples(M, CCF, sample_names, i):
    for s in range(S):
        out = open("data/" + exp_name + "/" + sample_names[s] + "/" + sample_names[s] + "_" + str(i) + ".csv", "w")
        out.write("Mutation_type,CCF,mut_name\n")
        
        temp_M = M[s]; temp_CCF = CCF[s]
        n_sorted = np.argsort(temp_CCF)[::-1]
        for n in range(N):
            out.write(str(int(temp_M[n_sorted[n]])) + "," + str(temp_CCF[n_sorted[n]]) 
                      + ",mut_" + str(int(n_sorted[n]+1)) + "\n")
        out.close()
            


if __name__ == "__main__":
    exp_name = sys.argv[1]
    ABS_path = os.getcwd() + "/"
    S, mean_K, N, J, CN_normal, read_counts_min, read_counts_max,\
        rho_state = util.load_simulation_setting(ABS_path, exp_name)
    known_names, signatures, mutation_names = util.load_known(ABS_path)
    V = len(mutation_names)

    main()

