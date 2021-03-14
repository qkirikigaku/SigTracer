import os
ABS_path = os.getcwd() + "/"
import sys
import utility as util
import select_best as sb
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()


def main():
    best_j = sb.select_best_j(ABS_path, exp_name, max_J, max_run, min_J = min_J)
    best_r = sb.select_best_r(ABS_path, exp_name, best_j, max_run_pe) 
    
    all_signature = []
    for s, sample_name in enumerate(sample_names):
        ref_name = ref_names[sample_name]
        for x in ref_name:
            if(x not in all_signature):
                all_signature.append(x)
    all_signature = list(np.sort(all_signature))
    all_signatures = util.load_all_sig_names(ABS_path)
    for x in all_signatures[::-1]:
        if(x not in all_signature): all_signatures.remove(x)
    all_signature = all_signatures
    K = len(all_signature)
    
    sig_by_clone = [[0.0 for k in range(K)] for i in range(2)]
    for s, sample_name in enumerate(sample_names):
        J = best_j[s]; R = best_r[s]
        rho, CCF = util.load_BB(ABS_path, exp_name, sample_name, J, R, "pe")
        pi = util.load_pi(ABS_path, exp_name, sample_name, J, R, "pe")
        sorted_j = np.argsort(CCF)[::-1]
        primary_j = sorted_j[0]
        for j in sorted_j:
            if(pi[j] > pi[primary_j] and CCF[j] > 0.95): primary_j = j
        qzU = util.load_qzU(ABS_path, exp_name, sample_name, J, R, "pe")
        N = len(qzU)
        activity = util.load_activity(ABS_path, exp_name, sample_name, J, R, "pe")
        temp_ref_names = ref_names[sample_name]
        for k,x in enumerate(temp_ref_names):
            target_k = all_signature.index(x)
            for j in range(J):
                if(j == primary_j): sig_by_clone[0][target_k] += pi[j] * activity[j][k] * N
                else: sig_by_clone[1][target_k] += pi[j] * activity[j][k] * N
    
    sum_ = [0.0 for k in range(K)]
    for i in range(2):
        for k in range(K):
            sum_[k] += sig_by_clone[i][k]
    for i in range(2):
        for k in range(K):
            sig_by_clone[i][k] /= sum_[k]

    labels = ["Primary clone", "Subclones"]
    signatures = np.arange(K)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    bottom = [0.0 for k in range(K)]
    ax.barh(signatures, sig_by_clone[0], left=bottom, label="Primary clone")
    bottom = sig_by_clone[0]
    ax.barh(signatures, sig_by_clone[1], left=bottom, label="Subclones")
    ax.set_xlabel("Proportion of clones to which the mutation belongs for each signature")
    ax.set_ylabel("Signature")
    ax.set_xlim(0.0,1.0)
    ax.set_yticks(signatures)
    ax.set_yticklabels(all_signature)
    ax.invert_yaxis()
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.savefig(ABS_path + "result/" + exp_name + "/PostAnalysis/sig_by_clone.png", dpi=300,
                bbox_inches= "tight")
    plt.close(1)


if __name__ == '__main__':
    args = sys.argv
    exp_name = args[1]
    min_J = 1
    max_J = 5
    max_run = 3
    max_run_pe = 3

    sample_names, purity = util.load_sample_names_and_purity(ABS_path, exp_name)
    ref_names = util.load_ref_sig(ABS_path, exp_name)

    main()
