import os
import sys
import numpy as np
import utility as util
ABS_path = os.getcwd() + "/"

def main():
    ref_sig = util.load_all_sig_names(ABS_path)
    sig_prob,mut_names = make_sig_prob(ref_sig)
    sample_names, purity = util.load_sample_names_and_purity(ABS_path, exp_name)

    for s, sample_name in enumerate(sample_names):
        lines = open("data/"+ exp_name + "/" + sample_name + ".csv", "r").readlines()[1:]
        out = open("data/" + exp_name + "/" + sample_name + ".csv", "w")
        out.write("mutation_id,chromosome,position,ref_counts,var_counts,normal_cn,minor_cn,mut_cn,major_cn,total_cn,trinucleotide,annotation\n")
        for line in lines[1:]:
            out.write(line[:line.index("SBS")])
            temp_sig = line[line.index("SBS"):]
            temp_k = ref_sig.index(temp_sig[:temp_sig.index(",")])
            temp_mut = mut_names[list(np.random.multinomial(1,sig_prob[temp_k])).index(1)]
            out.write(temp_mut + "\n")


def make_sig_prob(ref_sig):
    K = len(ref_sig)
    N = 81
    sig_prob = [[0.0 for n in range(N)] for k in range(K)]
    for k in range(K):
        for n in range(N):
            if(k==n): sig_prob[k][n] = 0.2
            else: sig_prob[k][n] = 0.01
    mut_names = ["" for n in range(N)]
    for i,sig in enumerate(ref_sig):
        mut_names[i] = sig + "-target"
    for i in range(K, N, 1):
        mut_names[i] = "other-" + str(i+1-K)
    return sig_prob, mut_names


if __name__ == "__main__":
    args = sys.argv
    exp_name = args[1]
    main()
