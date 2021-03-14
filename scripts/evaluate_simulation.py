import pickle
import utility as util
import itertools
from scipy.spatial.distance import cosine
import numpy as np
from scipy.special import betainc
from scipy.special import loggamma


def evaluate_model_selection(True_J, best_j, min_J, max_J, sample_names, result_path):
    out_all = open(result_path + "/tables/summary/model_selection.tsv", "w")
    out_all.write("Sample_name\tPredicted_J\n")
    summary_table = [0 for j in range(max_J-min_J+1)]
    for s,sample_name in enumerate(sample_names):
        out_all.write(sample_name + "\t" + str(best_j[s]) + "\n")
        summary_table[best_j[s]-min_J] += 1

    out = open(result_path + "/tables/summary/model_selection_summary.tsv", "w")
    out.write("Method")
    for j in range(min_J, max_J+1):
        out.write("\tJ=" + str(j))
    out.write("\nSigTracer")
    for j in range(max_J-min_J+1):
        out.write("\t" + str(summary_table[j]))
    out.write("\n")


def evaluate_accuracy(ABS_path, exp_name, J, best_r, ref_names, sample_names, purity, result_path):
    # Evaluate CCF accuracy
    gt_CCF = pickle.load(open(ABS_path + "data/" + exp_name + "/ground_truth.pkl", "rb"))["CCF"]
    out = open(result_path + "/tables/summary/delta_CCF.tsv", "w")
    out.write("Sample_name\tΔCCF\n")
    for s,sample_name in enumerate(sample_names):
        rho, CCF = util.load_BB(ABS_path, exp_name, sample_name, J, best_r[s], "pe")
        temp_delta_CCF = calc_delta_CCF(CCF, gt_CCF[s])
        out.write(sample_name + "\t" + str(temp_delta_CCF) + "\n")

    # Evaluate activity accuracy
    gt_activity = pickle.load(open(ABS_path + "data/" + exp_name + "/ground_truth.pkl",
                                   "rb"))["activity"]
    out = open(result_path + "/tables/summary/cos_activity.tsv", "w")
    out.write("Sample_name\tCosine_distance_of_activity\n")
    for s,sample_name in enumerate(sample_names):
        activity = util.load_activity(ABS_path, exp_name, sample_name, J, best_r[s], "pe")
        temp_cos = calc_cos_activity(activity, gt_activity[s])
        out.write(sample_name + "\t" + str(temp_cos) + "\n")

    # Evaluate Reconstruction-Rate for mutation-type
    out = open(result_path + "/tables/summary/RR_mut.tsv", "w")
    out.write("Sample_name\tRR_mutaiotn-type\n")
    for s, sample_name in enumerate(sample_names):
        qz = util.load_qz(ABS_path, exp_name, sample_name, J, best_r[s], "pe")
        mutations = util.load_mut_catalog(ABS_path, exp_name, sample_name)["trinucleotide"]
        signature = util.load_signature(ABS_path, exp_name, sample_name, J, best_r[s], "pe")
        temp_RR = calc_RR_mut(mutations, qz, signature, J)
        out.write(sample_name + "\t" + str(temp_RR) + "\n")

    # Evaluate Reconstruction-Rate for VAF
    out = open(result_path + "/tables/summary/RR_VAF.tsv", "w")
    out.write("Sample_name\tRR_VAF\n")
    for s, sample_name in enumerate(sample_names):
        mut_catalog = util.load_mut_catalog(ABS_path, exp_name, sample_name)
        pi = util.load_pi(ABS_path, exp_name, sample_name, J, best_r[s], "pe")
        rho,CCF = util.load_BB(ABS_path, exp_name, sample_name, J, best_r[s], "pe")
        qUC = util.load_qUC(ABS_path, exp_name, sample_name, J, best_r[s], "pe")
        qU = util.load_qU(ABS_path, exp_name, sample_name, J, best_r[s], "pe")
        temp_RR = calc_RR_VAF(mut_catalog["var_counts"], mut_catalog["ref_counts"],\
                              mut_catalog["total_cn"], purity[s], pi, rho, CCF, qUC,\
                              J, qU)
        out.write(sample_name + "\t" + str(temp_RR) + "\n")
    
    # Evaluate log-likelihood
    out = open(result_path + "/tables/summary/log_likelihood.tsv", "w")
    out.write("Sample_name\tlog_likelihood\n")
    for s, sample_name in enumerate(sample_names):
        mut_catalog = util.load_mut_catalog(ABS_path, exp_name, sample_name)
        qU = util.load_qU(ABS_path, exp_name, sample_name, J, best_r[s], "pe")
        qzU = util.load_qzU(ABS_path, exp_name, sample_name, J, best_r[s], "pe")
        qUz = convert_qzU_to_qUz(qzU)
        qUC = util.load_qUC(ABS_path, exp_name, sample_name, J, best_r[s], "pe")
        rho,CCF = util.load_BB(ABS_path, exp_name, sample_name, J, best_r[s], "pe")
        signature = util.load_signature(ABS_path, exp_name, sample_name, J, best_r[s], "pe")
        temp_LL = calc_log_likelihood(mut_catalog["var_counts"], mut_catalog["ref_counts"],\
                                      mut_catalog["total_cn"], mut_catalog["trinucleotide"],\
                                      mut_catalog["major_cn"], signature, purity[s], rho, CCF, qU, qUz, qUC, J)
        out.write(sample_name + "\t" + str(temp_LL) + "\n")


def calc_delta_CCF(predicted, correct):
    J = len(predicted)
    j_candidates = list(itertools.permutations(list(range(J)), J))
    total_distance = [0.0 for l in range(len(j_candidates))]
    for l, j_candidate in enumerate(j_candidates):
        for j, sorted_j in enumerate(j_candidate):
            total_distance[l] += abs(correct[j] - predicted[sorted_j]) / J
    min_distance = min(total_distance)
    return min_distance


def calc_cos_activity(predicted, correct):
    J = len(predicted)
    j_candidates = list(itertools.permutations(list(range(J)), J))
    total_distance = [0.0 for l in range(len(j_candidates))]
    for l, j_candidate in enumerate(j_candidates):
        for j, sorted_j in enumerate(j_candidate):
            total_distance[l] += cosine(correct[j], predicted[sorted_j]) / J
    min_distance = min(total_distance)
    return min_distance


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


def calc_RR_VAF(B, D, CN_tumor, purity, pi, rho, CCF, qUC, J, qU, div_num=30):
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
    for n in range(N):
        for div in range(div_num):
            for j in range(J):
                temp = 0.0
                for c in range(len(eta[n][j])):
                    if(eta[n][j][c] != 0.0):
                        temp += qUC[n][j][c] * (betainc(rho[j]*CCF[j]*eta[n][j][c],rho[j]*(1.0-CCF[j]*eta[n][j][c]),left[div+1])\
                                - betainc(rho[j]*CCF[j]*eta[n][j][c], rho[j]*(1.0-CCF[j]*eta[n][j][c]), left[div]))
                rec[div] += temp


    mut = [0 for div in range(div_num)]
    for n in range(N):
        temp_VAF = float(B[n])/float(B[n] + D[n])
        temp_div = int(temp_VAF/(1.0/div_num))
        if(temp_div == div_num): temp_div = div_num-1
        mut[temp_div] += 1
    
    RR = 0.0
    for div in range(div_num): RR += min(mut[div], rec[div])
    RR /= N
    return RR


def convert_qzU_to_qUz(qzU):
    N = len(qzU); K = len(qzU[0]); J = len(qzU[0][0])
    qUz = [[[0.0 for k in range(K)] for j in range(J)] for n in range(N)]
    for n in range(N):
        for j in range(J):
            for k in range(K):
                qUz[n][j][k] = qzU[n][k][j]
    return qUz


def calc_log_likelihood(B, D, CN_tumor, mutations, CN_major, signature, purity, rho, CCF, qU, qUz, qUC, J):
    LL = 0.0
    N = len(B); K = len(signature)
    for n in range(N):
        for j in range(J):
            temp_sum = sum(qUz[n][j])
            if(temp_sum != 0.0):
                for k in range(K): qUz[n][j][k] /= temp_sum
            else:
                for k in range(K): qUz[n][j][k] = 1.0/K
            temp_sum = sum(qUC[n][j])
            if(temp_sum != 0.0):
                for c in range(len(qUC[n][j])): qUC[n][j][c] /= temp_sum
            else: 
                for c in range(len(qUC[n][j])): qUC[n][j][c] = 1.0/len(qUC[n][j])

    for n in range(N):
        temp = 0.0
        for j in range(J):
            for k in range(K):
                if(signature[k][mutations[n]] != 0.0): temp += qU[n][j] * qUz[n][j][k] * np.log(signature[k][mutations[n]])
        LL += temp
    
    eta = [[[]for j in range(J)] for n in range(N)]
    for n in range(N):
        for j in range(J):
            for c in range(CN_major[n]):
                eta[n][j].append((purity * (1+c)) / ((purity * CN_tumor[n]) + ((1-purity) * 2.0)))

    for n in range(N):
        temp = 0.0
        for j in range(J):
            for c in range(len(eta[n][j])):
                x = BetaBinomial(B[n], B[n]+D[n], rho[j], CCF[j] * eta[n][j][c])
                temp += qU[n][j] * qUC[n][j][c] * x
        LL += temp
    return LL


def BetaBinomial(B, D, M, mu):
    log_likelihood = loggamma(D+1) - loggamma(B+1) - loggamma(D-B+1)
    log_likelihood += loggamma(M) - loggamma(D+M) + loggamma(B+M*mu) + loggamma(D-B+M*(1-mu))
    log_likelihood -= loggamma(M*mu) + loggamma(M*(1-mu))
    return log_likelihood

