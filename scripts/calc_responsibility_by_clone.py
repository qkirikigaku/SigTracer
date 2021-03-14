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
    
    all_signature = [];
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
    
    sig_by_clone = [[0.0 for i in range(2)] for k in range(K)]
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
                if(j == primary_j): sig_by_clone[target_k][0] += pi[j] * activity[j][k] * N
                else: sig_by_clone[target_k][1] += pi[j] * activity[j][k] * N
    
    sum_ = [0.0 for i in range(2)]
    for i in range(2):
        for k in range(K):
            sum_[i] += sig_by_clone[k][i]
    for i in range(2):
        for k in range(K):
            sig_by_clone[k][i] /= sum_[i]

    labels = ["Primary clone", "Subclones"]
    signatures = np.arange(2)

    raw_palette = sns.color_palette("hls", K)
    sig_palette = []
    for k in range(0, K, 2): sig_palette.append(raw_palette[k])
    for k in range(1, K, 2): sig_palette.append(raw_palette[k])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    bottom = [0.0 for i in range(2)]
    for k in range(K):
        ax.barh(signatures, sig_by_clone[k], left=bottom, color=sig_palette[k], label=all_signature[k])
        for i in range(2): bottom[i] += sig_by_clone[k][i]
    ax.set_xlabel("Proportion of signatures")
    ax.set_ylabel("Clone type")
    ax.set_xlim(0.0,1.0)
    ax.set_yticks(signatures)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.savefig(ABS_path + "result/" + exp_name + "/PostAnalysis/clone_by_sig.png", dpi=300,
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
