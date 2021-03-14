import os
import utility as util
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
import multiprocessing as mp
import numpy as np
from scipy.stats import binom_test


def test(ABS_path, exp_name, best_j, best_r, ref_names):
    sample_names, purity = util.load_sample_names_and_purity(ABS_path, exp_name)
    fig_dir, tab_dir = make_result_dirs(ABS_path, exp_name, sample_names)
    
    all_sigs = make_all_sigs(ABS_path, ref_names)

    sig_clone = []; mutations = []; primary_clones = []
    for s,sample_name in enumerate(sample_names):
        temp_sig_clone, temp_mutations, temp_primary_clones =\
            load_clone_info(ABS_path, exp_name, sample_name, purity[s],\
                            best_j[s], best_r[s], ref_names[sample_name])
        sig_clone.append(temp_sig_clone)
        mutations.append(temp_mutations)
        primary_clones.append(temp_primary_clones)
    
    all_mutations = make_all_mutations(mutations)
    positives, negatives, null_p = make_set(sig_clone, mutations, primary_clones, all_sigs, all_mutations,\
                                    sample_names)
    
    importance = make_test(positives, negatives, null_p, all_sigs, all_mutations)
    write_test(importance, all_mutations, all_sigs, tab_dir)


def make_result_dirs(ABS_path, exp_name, sample_names):
    fig_dir = ABS_path + "result/" + exp_name + "/figures/"
    tab_dir = ABS_path + "result/" + exp_name + "/tables/"
    for sample_name in sample_names:
        os.makedirs(fig_dir + sample_name, exist_ok=True)
        os.makedirs(tab_dir + sample_name, exist_ok=True)
    os.makedirs(fig_dir + "summary", exist_ok=True)
    os.makedirs(tab_dir + "summary", exist_ok=True)
    return fig_dir, tab_dir


def make_all_sigs(ABS_path, ref_names):
    all_sigs = []
    for temp_key in list(ref_names.keys()):
        temp_ref_names = ref_names[temp_key]
        for x in temp_ref_names:
            if(x not in all_sigs): all_sigs.append(x)
    all_sigs.sort()
    all_sig_ref = util.load_all_sig_names(ABS_path)
    for x in all_sig_ref[::-1]:
        if(x not in all_sigs): all_sig_ref.remove(x)
    all_sigs = all_sig_ref
    return all_sigs


def load_clone_info(ABS_path, exp_name, sample_name, purity,\
                    J, R, ref_sig):
    mut_catalog = util.load_mut_catalog(ABS_path, exp_name, sample_name)
    qU = util.load_qU(ABS_path, exp_name, sample_name, J, R, "pe")
    activity = util.load_activity(ABS_path, exp_name, sample_name, J, R, "pe")
    pi = util.load_pi(ABS_path, exp_name, sample_name, J, R, "pe")
    rho,CCF = util.load_BB(ABS_path, exp_name, sample_name, J, R, "pe")
    qC = util.load_qC(ABS_path, exp_name, sample_name, J, R, "pe")
    B = mut_catalog["var_counts"]; D = mut_catalog["ref_counts"]
    CN_total = mut_catalog["total_cn"]; annotation = mut_catalog["annotation"]

    N = len(qU); K = len(ref_sig)
    M = [0.0 for n in range(N)]
    for n in range(N):
        for c in range(len(qC[n])):
            M[n] += (1+c) * qC[n][c]
    eta = [0.0 for n in range(N)]
    for n in range(N):
        eta[n] = (purity * M[n]) / (purity * CN_total[n] + (1.0-purity) * 2.0)
    expected_CCF = [0.0 for n in range(N)]
    for n in range(N):
        expected_CCF[n] = (B[n]/(B[n]+D[n])) / eta[n]
    sorted_n = np.argsort(expected_CCF)[::-1]
    
    mutations = [[] for j in range(J)]
    for n in sorted_n:
        temp_J = np.argmax(qU[n])
        mutations[temp_J].append(annotation[n])

    sig_clone = [[] for j in range(J)]
    for j in range(J):
        for k in range(K):
            if(N * pi[j] * activity[j][k] > 100): sig_clone[j].append(ref_sig[k])
    
    sorted_j = list(np.argsort(CCF)[::-1])
    primary_j = sorted_j[0]
    for j in sorted_j:
        if(pi[j] > pi[primary_j] and CCF[j] > 0.95): primary_j = j
    primary_clones = [False for j in range(J)]
    primary_clones[primary_j] = True

    return sig_clone, mutations, primary_clones


def make_all_mutations(mutations):
    all_mutations = []
    for x in mutations:
        for y in x:
            for z in y:
                if(z not in all_mutations): all_mutations.append(z)
    all_mutations.sort()
    return all_mutations


def make_set(sig_clone, mutations, primary_clones, all_sigs, all_mutations, sample_names,\
             initiate_rate = 1.0, clone_threshold = 100):
    K = len(all_sigs)
    positive_set = [{} for k in range(K)]
    negative_set = [{} for k in range(K)]
    null_p = [[0,0] for k in range(K)]
    for s, sample_name in enumerate(sample_names):
        for j,temp_clone_mutations in enumerate(mutations[s]):
            temp_num_mutation = len(temp_clone_mutations)
            max_initiate_index = int(len(temp_clone_mutations) * initiate_rate)
            for m, temp_mutation in enumerate(temp_clone_mutations):
                for k, sig in enumerate(all_sigs):
                    if(m < max_initiate_index and temp_num_mutation > clone_threshold and sig in
                       sig_clone[s][j] and primary_clones[s][j]):
                        null_p[k][0] += 1
                        if(temp_mutation not in positive_set[k].keys()):
                            positive_set[k].update({temp_mutation:1})
                        else:
                            temp = positive_set[k][temp_mutation]
                            positive_set[k].update({temp_mutation:temp+1})
                    else:
                        null_p[k][1] += 1
                        if(temp_mutation not in negative_set[k].keys()):
                            negative_set[k].update({temp_mutation:1})
                        else:
                            temp = negative_set[k][temp_mutation]
                            negative_set[k].update({temp_mutation:temp+1})
    return positive_set, negative_set, null_p


def make_test(positives, negatives, null_p, all_sigs, all_mutations):
    M = len(all_mutations); K = len(all_sigs)
    importance = np.zeros([M,K])
    for k in range(K):
        temp_null_p = null_p[k][0]/(null_p[k][0]+null_p[k][1])
        for m, mutation in enumerate(all_mutations):
            positive = 0
            if(mutation in positives[k].keys()):
                positive = positives[k][mutation]
            trials = positive
            if(mutation in negatives[k].keys()):
                trials += negatives[k][mutation]
            importance[m,k] = binom_test(positive, trials, p=temp_null_p, alternative="greater")
    return importance


def write_test(mean_test_result, all_mutations, all_sigs, tab_dir):
    out = open(tab_dir + "summary/test_mutation.tsv", "w")
    out.write("Mutation")
    for sig in all_sigs: out.write("\t" + sig)
    out.write("\n")
    for m, mut in enumerate(all_mutations):
        out.write(mut)
        for s,sig in enumerate(all_sigs):
            out.write("\t" + str(mean_test_result[m,s]))
        out.write("\n")
    out.close()

