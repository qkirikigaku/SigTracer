import numpy as np
import csv
import os
import pickle

def load_known(ABS_path, ver="3.1"):
    input_file = ABS_path + "data/ref/signature_probability_v" + ver + ".csv"
    lines = open(input_file, "r").readlines()
    known_names = lines[0].split(",")[2:]
    known_names[-1] = known_names[-1][:-1]
    K = len(known_names); V = len(lines)-1
    signatures = np.zeros([K,V]); mutation_names = []
    for v, line in enumerate(lines[1:]):
        temp_list = line.split(",")
        temp_list[-1] = temp_list[-1][:-1]
        mutation_name = temp_list[1][0] + "[" + temp_list[0] + "]" + temp_list[1][2]
        mutation_names.append(mutation_name)
        for k, temp in enumerate(temp_list[2:]):
            signatures[k,v] = float(temp)
    return known_names, signatures, mutation_names
    

def load_simulation_setting(ABS_path, exp_name):
    lines = open(ABS_path + "data/" + exp_name + "/setting.txt", "r").readlines()
    S = int(lines[0].split()[1])
    mean_K = int(lines[1].split()[1])
    N = int(lines[2].split()[1])
    J = int(lines[3].split()[1])
    CN_normal = int(lines[4].split()[1])
    read_counts_min = int(lines[5].split()[1])
    read_counts_max = int(lines[6].split()[1])
    rho_state = lines[7].split()[1]
    return S, mean_K, N, J, CN_normal, read_counts_min,\
           read_counts_max, rho_state


def load_MT(ABS_path, mt_type = "SBS"):
    lines = open(ABS_path + "data/ref/MT_" + mt_type + ".txt", "r").readlines()
    mt_names = []
    for line in lines:
        mt_names.append(line[:-1])
    return mt_names


def load_sex(ABS_path, aliquot_ids, sample_names):
    lines = open(ABS_path + "raw_data/pcawg_sample_sheet.tsv", "r").readlines()
    ali_to_donor = {}
    for line in lines[1:]:
        temp = line.split("\t")
        ali_to_donor.update({temp[5]:temp[0]})
    lines = open(ABS_path + "raw_data/pcawg_donor_clinical_August2016_v9.csv", "r").readlines()
    donor_to_sex = {}
    for line in lines[1:]:
        temp = line.split(",")
        donor_to_sex.update({temp[0]:temp[5]})

    sex_dic = {}; no_sex_count = 0
    for sample_name,aliquot_id in zip(sample_names, aliquot_ids):
        temp_donor = ali_to_donor[aliquot_id]
        if(temp_donor in donor_to_sex.keys()):
            sex_dic.update({sample_name:donor_to_sex[temp_donor]})
        else:
            no_sex_count += 1
            sex_dic.update({sample_name:"undefined"})
    print("no_sex_count : " + str(no_sex_count))
    return sex_dic


def load_age(ABS_path, aliquot_ids, sample_names):
    lines = open(ABS_path + "raw_data/pcawg_sample_sheet.tsv", "r").readlines()
    ali_to_donor = {}
    for line in lines[1:]:
        temp = line.split("\t")
        ali_to_donor.update({temp[5]:temp[0]})
    lines = open(ABS_path + "raw_data/pcawg_donor_clinical_August2016_v9.csv", "r").readlines()
    donor_to_age = {}
    for line in lines[1:]:
        temp = line.split(",")
        if(temp[10] != ""): donor_to_age.update({temp[0]:int(temp[10])})

    age_dic = {}; no_age_count = 0
    for sample_name,aliquot_id in zip(sample_names, aliquot_ids):
        temp_donor = ali_to_donor[aliquot_id]
        if(temp_donor in donor_to_age.keys()):
            age_dic.update({sample_name:donor_to_age[temp_donor]})
        else:
            no_age_count += 1
            age_dic.update({sample_name:"undefined"})
    print("no_age_count : " + str(no_age_count))
    return age_dic


def load_ref_sig(ABS_path, exp_name):
    reader = csv.reader(open(ABS_path + "data/" + exp_name + "/ref_sig.csv"))
    lines = [line for line in reader]
    ref = {}
    for line in lines:
        temp_sample = line[0]
        ref.update({temp_sample:line[1:]})
    return ref


def load_ref_sig_simulation(ABS_path, exp_name):
    ground_truth = pickle.load(open(ABS_path + "data/" + exp_name + "/ground_truth.pkl", "rb"))
    ref_names = ground_truth["ref_sig"]
    return ref_names


def load_all_sig_names(ABS_path):
    reader = csv.reader(open(ABS_path + "raw_data/PCAWG_sigProfiler_SBS_signatures_in_samples.csv",
                "r"))
    temp_list = [line for line in reader][0][3:]
    return temp_list


def load_sample_names_and_purity(ABS_path, exp_name):
    reader = csv.reader(open(ABS_path + "data/" + exp_name + "/purity.csv"))
    lines = [line for line in reader]
    sample_names = []; purity = [];
    for line in lines[1:]:
        sample_names.append(line[0])
        purity.append(float(line[1]))
    return sample_names, purity


def define_color_list(mt_type = "SBS"):
    if(mt_type == "SBS"):
        color_list = ["_" for v in range(96)]
        for i in range(96):
            if(i < 16): color_list[i] = "r"
            elif(i < 32): color_list[i] = "g"
            elif(i < 48): color_list[i] = "b"
            elif(i < 64): color_list[i] = "c"
            elif(i < 80): color_list[i] = "m"
            else: color_list[i] = "y"
        return color_list


def load_sample_info(ABS_path, tumor_type):
    lines = open(ABS_path + "raw_data/" + tumor_type + "/sample_info.txt", "r").readlines()
    sample_names = []; aliquot_ids = []; purity = []
    for line in lines[1:]:
        temp_list = line.split("\t")
        sample_names.append(temp_list[0])
        aliquot_ids.append(temp_list[1])
        purity.append(float(temp_list[2]))
    return sample_names, aliquot_ids, purity


def load_elbo(ABS_path, exp_name, sample_name, J, R, mode, tag=""):
    file_name = ABS_path + "result/" + exp_name + "/" + sample_name + "/J" +\
                str(J) + "_run" + str(R) + "_" + mode + tag +\
                "_elbo.txt"
    elbo = -1e100
    if(os.path.exists(file_name)):
        lines = open(file_name, "r").readlines()
        elbo = float(lines[0].split()[0])
    return elbo


def load_pi(ABS_path, exp_name, sample_name, J, R, mode, tag=""):
    lines = open(ABS_path + "result/" + exp_name + "/" + sample_name + "/J" +\
                 str(J) + "_run" + str(R) + "_" + mode + tag + \
                 "_pi.txt", "r").readlines()
    pi = list(map(float, lines[0].split()))
    return pi


def load_BB(ABS_path, exp_name, sample_name, J, R, mode, tag=""):
    lines = open(ABS_path + "result/" + exp_name + "/" + sample_name + "/J" +\
                 str(J) + "_run" + str(R) + "_" + mode + tag +\
                 "_BB.txt", "r").readlines()
    rho = list(map(float, lines[0].split()))
    CCF = list(map(float, lines[1].split()))
    for j,x in enumerate(CCF):
        if(x == 1.0): CCF = 1.0-1e-10
    return rho, CCF


def load_activity(ABS_path, exp_name, sample_name, J, R, mode, tag=""):
    lines = open(ABS_path + "result/" + exp_name + "/" + sample_name + "/J" +\
                 str(J) + "_run" + str(R) + "_" + mode + tag +\
                 "_activity.txt", "r").readlines()
    activity = []
    for line in lines:
        activity.append(list(map(float, line.split())))
    return activity


def load_signature(ABS_path, exp_name, sample_name, J, R, mode, tag=""):
    lines = open(ABS_path + "result/" + exp_name + "/" + sample_name + "/J" +\
                 str(J) + "_run" + str(R) + "_" + mode + tag +\
                 "_signature.txt", "r").readlines()
    signature = []
    for line in lines:
        signature.append(list(map(float, line.split())))
    return signature


def load_mut_catalog(ABS_path, exp_name, sample_name):
    reader = csv.reader(open(ABS_path + "data/" + exp_name + "/" + sample_name + ".csv"))
    lines = [line for line in reader]
    header = lines[0]
    mut_catalog = {}
    for head in header: mut_catalog.update({head:[]})
    for line in lines[1:]:
        for i,temp in enumerate(line):
            if(header[i] in ["ref_counts", "var_counts", "normal_cn", "minor_cn",\
                             "major_cn", "total_cn", "trinucleotide"]):
                mut_catalog[header[i]].append(int(temp))
            else: mut_catalog[header[i]].append(temp)
    return mut_catalog


def load_qzU(ABS_path, exp_name, sample_name, J, R, mode, tag=""):
    lines = open(ABS_path + "result/" + exp_name + "/" + sample_name + "/J" +\
                 str(J) + "_run" + str(R) + "_" + mode + tag +\
                 "_qzU.txt", "r").readlines()
    header = lines[0].split()
    temp_N = int(header[0]); temp_K = int(header[1])
    qzU = [[] for n in range(temp_N)]; temp_n = 0
    for i,line in enumerate(lines[1:]):
        temp_k = i % temp_K
        temp_list = list(map(float, line.split()))
        qzU[temp_n].append(temp_list)
        if(temp_k == temp_K-1): temp_n += 1
    return qzU


def load_qUC(ABS_path, exp_name, sample_name, J, R, mode, tag=""):
    lines = open(ABS_path + "result/" + exp_name + "/" + sample_name + "/J" +\
                 str(J) + "_run" + str(R) + "_" + mode + tag +\
                 "_qUC.txt", "r").readlines()
    header = lines[0].split()
    temp_N = int(header[0])
    qUC = [[] for n in range(temp_N)]; temp_n = 0
    for i,line in enumerate(lines[2:]):
        temp_j = i % J
        temp_list = list(map(float, line.split()))
        qUC[temp_n].append(temp_list)
        if(temp_j == J-1): temp_n += 1
    return qUC


def load_qz(ABS_path, exp_name, sample_name, J, R, mode, tag=""):
    lines = open(ABS_path + "result/" + exp_name + "/" + sample_name + "/J" +\
                 str(J) + "_run" + str(R) + "_" + mode + tag +\
                 "_qz.txt", "r").readlines()
    qz = []
    for line in lines: qz.append(list(map(float, line.split())))
    return qz


def load_qU(ABS_path, exp_name, sample_name, J, R, mode, tag=""):
    lines = open(ABS_path + "result/" + exp_name + "/" + sample_name + "/J" +\
                 str(J) + "_run" + str(R) + "_" + mode + tag +\
                 "_qU.txt", "r").readlines()
    qU = []
    for line in lines:
        qU.append(list(map(float, line.split())))
    return qU


def load_qC(ABS_path, exp_name, sample_name, J, R, mode, tag=""):
    lines = open(ABS_path + "result/" + exp_name + "/" + sample_name + "/J" +\
                 str(J) + "_run" + str(R) + "_" + mode + tag +\
                 "_qC.txt", "r").readlines()
    qC = []
    for line in lines: qC.append(list(map(float, line.split())))
    return qC
