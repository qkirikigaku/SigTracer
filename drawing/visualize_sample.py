import utility as util
import numpy as np
from scipy.special import betainc
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()


def visualize_sample(ABS_path, exp_name, fig_dir, tab_dir, sample_name, purity, J, R, ref_sig):
    rho,CCF = util.load_BB(ABS_path, exp_name, sample_name, J, R, "pe")
    activity = util.load_activity(ABS_path, exp_name, sample_name, J, R, "pe")
    signature = util.load_signature(ABS_path, exp_name, sample_name, J, R, "pe")
    pi = util.load_pi(ABS_path, exp_name, sample_name, J, R, "pe")
    mut_catalog = util.load_mut_catalog(ABS_path, exp_name, sample_name)
    qzU = util.load_qzU(ABS_path, exp_name, sample_name, J, R, "pe")
    qUC = util.load_qUC(ABS_path, exp_name, sample_name, J, R, "pe")
    qz = util.load_qz(ABS_path, exp_name, sample_name, J, R, "pe")
    qU = util.load_qU(ABS_path, exp_name, sample_name, J, R, "pe")
    qC = util.load_qC(ABS_path, exp_name, sample_name, J, R, "pe")

    RR_mut = calc_RR_mut(mut_catalog["trinucleotide"], qz, signature, J)
    RR_vaf = calc_RR_vaf(mut_catalog["var_counts"], mut_catalog["ref_counts"],\
                         mut_catalog["total_cn"], purity, pi, rho, CCF, qUC, J,\
                         fig_dir, sample_name, qU)
    
    draw_clone(J, CCF, activity, pi, qU, qzU, ref_sig, fig_dir, sample_name,\
               mut_catalog["var_counts"], mut_catalog["ref_counts"], purity,\
               mut_catalog["total_cn"], qC)
    write_stats(tab_dir, sample_name, J, RR_mut, RR_vaf)


def calc_RR_mut(mutations, qz, signature, J):
    N = len(mutations); K = len(signature); V = len(signature[0])
    mut = [0 for v in range(V)]
    for x in mutations: mut[x] += 1
    rec = [0 for v in range(V)]
    for n in range(N):
        temp = np.dot(qz[n], signature)
        for v in range(V): rec[v] += temp[v]
    RR = 0.0
    for v in range(V): RR += min(mut[v], rec[v])
    RR /= N
    return RR


def calc_RR_vaf(B, D, CN_tumor, purity, pi, rho, CCF, qUC, J, fig_dir, sample_name, qU, div_num=30):
    N = len(B)
    M = [[[] for j in range(J)] for n in range(N)]
    for n in range(N):
        for j in range(J):
            for c in range(len(qUC[n][j])):
                M[n][j].append((1+c) * qUC[n][j][c])
    eta = [[[]for j in range(J)] for n in range(N)]
    for n in range(N):
        for j in range(J):
            for c in range(len(M[n][j])):
                eta[n][j].append((purity * M[n][j][c]) / ((purity * CN_tumor[n]) + ((1-purity) * 2.0)))

    left = [0 for div in range(div_num+1)]
    for div in range(1, div_num+1):
        left[div] = div * 1.0/div_num
    
    rec = [0.0 for div in range(div_num)]
    rec_wc = [[0.0 for div in range(div_num)] for j in range(J)]
    for n in range(N):
        for div in range(div_num):
            for j in range(J):
                temp = 0.0
                for c in range(len(eta[n][j])):
                    if(eta[n][j][c] != 0.0):
                        temp += qUC[n][j][c] * (betainc(rho[j]*CCF[j]*eta[n][j][c],rho[j]*(1.0-CCF[j]*eta[n][j][c]),left[div+1])\
                                - betainc(rho[j]*CCF[j]*eta[n][j][c], rho[j]*(1.0-CCF[j]*eta[n][j][c]), left[div]))
                rec[div] += temp
                rec_wc[j][div] += temp

    mut = [0 for div in range(div_num)]
    for n in range(N):
        temp_VAF = float(B[n])/float(B[n] + D[n])
        temp_div = int(temp_VAF/(1.0/div_num))
        if(temp_div == div_num): temp_div = div_num-1
        mut[temp_div] += 1
    RR = 0.0
    for div in range(div_num): RR += min(mut[div], rec[div])
    RR /= N
    
    sorted_j = np.argsort(CCF)[::-1]
    fig = plt.figure()
    ax = fig.add_subplot(211)
    bottom = [0 for div in range(div_num)]
    for j in range(J):
        ax.bar(left[1:], rec_wc[sorted_j[j]], bottom=bottom, width=1.0/div_num)
        for div in range(div_num): bottom[div] += rec_wc[sorted_j[j]][div]
    ax.set_title("Reconstructed")
    ax.set_ylabel("Mutation burden")
    ax.invert_xaxis()
    ax = fig.add_subplot(212)
    ax.bar(left[1:], mut, width=1.0/div_num)
    ax.set_title("Observed")
    ax.set_xlabel("VAF")
    ax.set_ylabel("Mutation burden")
    ax.invert_xaxis()
    fig.tight_layout()
    fig.savefig(fig_dir + sample_name + "/" + sample_name + "_VAF_Reconstructed.png")
    plt.close(1)
    
    return RR


def draw_clone(J, CCF, activity, pi, qU, qzU, ref_sig, fig_dir, sample_name,\
               var, ref, purity, CN_total, qC):
    N = len(qU); K = len(ref_sig)
    div_num = N // 100
    if(div_num > 100): div_num=100
    if(div_num < 10): div_num = 10

    M = [0.0 for n in range(N)]
    for n in range(N):
        for c in range(len(qC[n])):
            M[n] += (1+c) * qC[n][c]
    eta = [0.0 for n in range(N)]
    for n in range(N):
        eta[n] = (purity * M[n]) / (purity * CN_total[n] + (1.0-purity) * 2.0)

    expected_CCF = [0.0 for n in range(N)]
    for n in range(N):
        expected_CCF[n] = (var[n] / (var[n] + ref[n])) / eta[n]
    max_CCF = max(expected_CCF)

    left = np.linspace(0.0, max_CCF, div_num)

    rec = [[0.0 for div in range(div_num)] for j in range(J)]
    rec_with_K = [[[0.0 for div in range(div_num)] for k in range(K)] for j in range(J)]
    for n in range(N):
        temp_div = int(expected_CCF[n]/(max_CCF/div_num))
        if(temp_div >= div_num): temp_div = div_num-1
        for j in range(J):
            rec[j][temp_div] += qU[n][j]
            for k in range(K):
                rec_with_K[j][k][temp_div] += qzU[n][k][j]

    raw_palette = sns.color_palette("hls", K)
    sig_palette = []
    for k in range(0, K, 2): sig_palette.append(raw_palette[k])
    for k in range(1, K, 2): sig_palette.append(raw_palette[k])

    """
    sig_left = [[] for j in range(J)]; sig_height = [[[] for k in range(K)] for j in range(J)]

    sorted_n = np.argsort(expected_CCF)[::-1]
    temp_CCF = [[] for j in range(J)]; temp_qU = [[] for j in range(J)]
    temp_qzU = [[[] for k in range(K)] for j in range(J)]
    pop_flag = [False for j in range(J)]
    for n in sorted_n:
        for j in range(J):
            temp_CCF[j].append(expected_CCF[n])
            temp_qU[j].append(qU[n][j])
            for k in range(K): temp_qzU[j][k].append(qzU[n][k][j])
        for j in range(J):
            if(sum(temp_qU[j]) > 50):
                temp_CCF_sum = 0.0
                for i,x in enumerate(temp_CCF[j]):
                    temp_CCF_sum += temp_qU[j][i] * x
                temp_CCF_sum /= sum(temp_qU[j])
                sig_left[j].append(temp_CCF_sum)
                for k in range(K): sig_height[j][k].append(np.mean(temp_qzU[j][k]))
                pop_flag[j] = True
            else:
                pop_flag[j] = False
            if(pop_flag[j]):
                temp_CCF[j].pop(0)
                temp_qU[j].pop(0)
                for k in range(K): temp_qzU[j][k].pop(0)
    for j in range(J):
        if(len(sig_left[j]) == 0):
            sig_left[j].append(0.0)
            for k in range(K): sig_height[j][k].append(activity[j][k])
        else:
            sig_left[j].append(0.0)
            for k in range(K): sig_height[j][k].append(sig_height[j][k][-1])
        sig_left[j].insert(0, max_CCF)
        for k in range(K): sig_height[j][k].insert(0, sig_height[j][k][0])
        for div in range(len(sig_left[j])):
            temp_sum = 0.0;
            for k in range(K): temp_sum += sig_height[j][k][div]
            for k in range(K): sig_height[j][k][div] /= temp_sum
    """

    sorted_j = np.argsort(CCF)[::-1]
    
    fig = plt.figure()
    ax = fig.add_subplot(1+J,1,1)
    bottom = [0 for div in range(div_num)]
    for j in range(J):
        ax.bar(left, rec[sorted_j[j]], bottom=bottom, width=max_CCF/div_num, label="Clone "+ str(j+1))
        for div in range(div_num): bottom[div] += rec[sorted_j[j]][div]
    ax.set_title("Clone composition (" + str(N) + " mutations)")
    ax.set_ylabel("Mutation")
    ax.set_xlim(0.0, max_CCF)
    ax.invert_xaxis()
    ax.set_xticklabels([])
    ax.legend(bbox_to_anchor=(1.01,1), loc="upper left")
    plt.subplots_adjust(hspace=0.5)
    for j in range(J):
        ax = fig.add_subplot(1+J,1,2+j)
        left = np.linspace(0.0, max_CCF, div_num)
        bottom = [0.0 for div in range(div_num)]
        for k in range(K):
            ax.bar(left, rec_with_K[sorted_j[j]][k], bottom=bottom, width=max_CCF/div_num,
                   color=sig_palette[k], label=ref_sig[k])
            for div in range(div_num): bottom[div] += rec_with_K[sorted_j[j]][k][div]
        ax.set_ylabel("Exposure")
        ax.set_xlim(0.0, max_CCF)
        """
        left = [0.0, max_CCF]
        height = [[activity[sorted_j[j]][k] for x in range(2)] for k in range(K)]
        for k in range(K):
            ax.plot(left, height[k], linestyle = "dashed", alpha=0.5, color=sig_palette[k])
        for k in range(K):
            ax.plot(sig_left[sorted_j[j]], sig_height[sorted_j[j]][k], color=sig_palette[k], label=ref_sig[k])
        ax.set_ylabel("Activity")
        ax.set_xlim(0.0,max_CCF); ax.set_ylim(0.0, 1.05)
        """
        if(j == J-1):
            ax.set_xlabel("Expected cancer cell fraction (CCF)")
        ax.invert_xaxis()
        if((j==0 and J in [1,2,3]) or (j==1 and J>=4)): ax.legend(bbox_to_anchor=(1.01,1), loc="upper left")
        if(j != J-1):
            ax.set_xticklabels([])
        ax.set_title("Clone " + str(j+1) + " (CCF : " + str(round(CCF[sorted_j[j]], 3)) + ", proportion : " + str(round(pi[sorted_j[j]], 3)) + ")")
    fig.savefig(fig_dir + sample_name + "/" + sample_name + "_clone.png", bbox_inches="tight")
    plt.close(1)


def write_stats(tab_dir, sample_name, J, RR_mut, RR_vaf):
    out = open(tab_dir + sample_name + "/stats.tsv", "w")
    out.write("The_number_of_clones\t" + str(J) + "\n")
    out.write("Reconstruction_rate_for_mutational_type\t" + str(RR_mut) + "\n")
    out.write("Reconstruction_rate_for_VAF\t" + str(RR_vaf) + "\n")
    out.close()
