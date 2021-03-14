import os
import sys
ABS_path = os.getcwd() + "/"

import numpy as np
import utility as util
import select_best as sb
import evaluate_simulation as es
import pickle
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()


def main():
    best_j = sb.select_best_j(ABS_path, exp_name, max_J, max_run, min_J = min_J)
    best_rs = []
    for j in range(min_J, max_J+1):
        temp_best_r = []
        for s, sample_name in enumerate(sample_names):
            ELBO = []
            for r in range(1, max_run+1):
                ELBO.append(util.load_elbo(ABS_path, exp_name, sample_name,\
                                           j, r, "ms"))
            r_selected = np.argmax(ELBO) + 1
            temp_best_r.append(r_selected)
        best_rs.append(temp_best_r)

    color_palette = sns.color_palette("hls", max_J-min_J+1)
    colors = [color_palette[best_j[s]-1] for s in range(S)]

    LL = [[0.0 for s in range(S)] for j in range(max_J-min_J+2)]

    for j in range(min_J, max_J+1):
        for s,sample_name in enumerate(sample_names):
            mut_catalog = util.load_mut_catalog(ABS_path, exp_name, sample_name)
            qU = util.load_qU(ABS_path, exp_name, sample_name, j, best_rs[j-min_J][s], "ms")
            qzU = util.load_qzU(ABS_path, exp_name, sample_name, j, best_rs[j-min_J][s], "ms")
            qUz = es.convert_qzU_to_qUz(qzU)
            qUC = util.load_qUC(ABS_path, exp_name, sample_name, j, best_rs[j-min_J][s], "ms")
            rho,CCF = util.load_BB(ABS_path, exp_name, sample_name, j, best_rs[j-min_J][s], "ms")
            signature = util.load_signature(ABS_path, exp_name, sample_name, j, best_rs[j-min_J][s], "ms")
            temp_LL = es.calc_log_likelihood(mut_catalog["var_counts"], mut_catalog["ref_counts"],\
                                      mut_catalog["total_cn"], mut_catalog["trinucleotide"],\
                                      mut_catalog["major_cn"], signature, purity[s], rho, CCF, qU,\
                                      qUz, qUC, j)
            LL[j-min_J][s] = temp_LL

    gt = pickle.load(open(ABS_path + "data/" + exp_name + "/ground_truth.pkl", "rb"))
    gt_pi = gt["pi"]
    gt_CCF = gt["CCF"]; gt_rho = gt["rho"]
    gt_signature = gt["signature"]
    for s, sample_name in enumerate(sample_names):
        gt_z = util.load_mut_catalog(ABS_path, exp_name, sample_name)["signature"]
        qz = convert_gt_z(gt_z, ref_names[s])
        gt_U = list(map(int, util.load_mut_catalog(ABS_path, exp_name, sample_name)["clone"]))
        gt_C = list(map(int, util.load_mut_catalog(ABS_path, exp_name, sample_name)["mut_cn"]))
        gt_major = list(map(int, util.load_mut_catalog(ABS_path, exp_name,
                            sample_name)["major_cn"]))
        signature = gt_signature[s]
        rho = gt_rho[s]; CCF = gt_CCF[s]
        qUC, qU = convert_gt_UC(gt_U, gt_C, gt_major, J)
        qUz = make_qUz(qz, qU)
        mut_catalog = util.load_mut_catalog(ABS_path, exp_name, sample_name)
        temp_LL = es.calc_log_likelihood(mut_catalog["var_counts"], mut_catalog["ref_counts"],\
                                      mut_catalog["total_cn"], mut_catalog["trinucleotide"],\
                                      mut_catalog["major_cn"], signature, purity[s], rho, CCF, qU, qUz, qUC, J)
        LL[-1][s] = temp_LL


    fig = plt.figure()
    left = list(range(max_J-min_J+2))
    settings = ["J=" + str(j) for j in range(min_J, max_J+1)]
    settings.append("ground-truth")
    
    ax_ll = fig.add_subplot(111)
    sns.boxplot(data=LL, showfliers=False, ax=ax_ll)
    #sns.stripplot(data=RR_mut, jitter=True, color=colors, size=3, ax=ax_ll)
    sns.stripplot(data=LL, jitter=True, color=".3", size=3, ax=ax_ll)
    ax_ll.set_xticks(left)
    ax_ll.set_xticklabels(settings)
    ax_ll.set_ylabel("Log-likelihood")
    ax_ll.set_xlabel("# clones and ground-truth")
    fig.tight_layout()
    plt.savefig(ABS_path + "result/" + exp_name + "/figures/summary/likelihood_by_J.png", dpi=300)
    plt.close(1)
    

def convert_gt_z(gt_z, ref_name):
    qz = [[0.0 for k in range(len(ref_name))] for n in range(len(gt_z))]
    for n,x in enumerate(gt_z):
        qz[n][ref_name.index(x)] = 1.0
    return qz


def convert_gt_UC(gt_U, gt_C, gt_major, temp_j):
    qUC = [[[0.0 for c in range(gt_major[n])] for j in range(temp_j)] for n in range(len(gt_U))]
    qU = [[0.0 for j in range(temp_j)] for n in range(len(gt_U))]
    for n,x in enumerate(gt_U):
        qU[n][x] = 1.0
        qUC[n][x][gt_C[n]-1] = 1.0
    return qUC, qU


def make_qUz(qz, qU):
    N = len(qz); K = len(qz[0]); temp_J = len(qU[0])
    qUz = [[[0.0 for k in range(K)] for j in range(temp_J)] for n in range(N)]
    for n in range(N):
        for j in range(temp_J):
            qUz[n][j][qz[n].index(1.0)] = 1.0
    return qUz


if __name__ == "__main__":
    args = sys.argv
    exp_name = args[1]
    min_J = 1
    max_J = 4
    max_run = 2

    S,mean_K, N,J, CN_normal, read_counts_min, read_counts_max, rho_state =\
        util.load_simulation_setting(ABS_path, exp_name)
    ref_names = util.load_ref_sig_simulation(ABS_path, exp_name)
    sample_names, purity = util.load_sample_names_and_purity(ABS_path, exp_name)
    
    main()
