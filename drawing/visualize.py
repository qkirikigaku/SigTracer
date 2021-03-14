import os
import utility as util
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
import multiprocessing as mp
import visualize_sample as vs
import numpy as np


def visualize_result(ABS_path, exp_name, best_j, best_r, ref_names, num_processes):
    sample_names, purity = util.load_sample_names_and_purity(ABS_path, exp_name)
    fig_dir, tab_dir = make_result_dirs(ABS_path, exp_name, sample_names)
    draw_num_clones(fig_dir, best_j, sample_names)
    
    arguments = []
    for s,sample_name in enumerate(sample_names):
        vs.visualize_sample(ABS_path, exp_name, fig_dir, tab_dir, sample_name,\
                            purity[s], best_j[s], best_r[s], ref_names[sample_name])

    summarize_RR(ABS_path, exp_name, tab_dir, sample_names, best_j)


def make_result_dirs(ABS_path, exp_name, sample_names):
    fig_dir = ABS_path + "result/" + exp_name + "/figures/"
    tab_dir = ABS_path + "result/" + exp_name + "/tables/"
    for sample_name in sample_names:
        os.makedirs(fig_dir + sample_name, exist_ok=True)
        os.makedirs(tab_dir + sample_name, exist_ok=True)
    os.makedirs(fig_dir + "summary", exist_ok=True)
    os.makedirs(tab_dir + "summary", exist_ok=True)
    return fig_dir, tab_dir


def draw_num_clones(fig_dir, best_j, sample_names):
    min_j = min(best_j); max_j = max(best_j)
    left = range(min_j, max_j+1)
    height = [0 for i in range(max_j - min_j + 1)]
    for x in best_j: height[x-min_j] += 1

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(left, height)
    ax.set_xlabel("The predicted number of clones")
    ax.set_ylabel("The number of samples")
    ax.set_title("Model selection (total : " + str(len(sample_names)) + " samples)")
    ax.set_xticks(left); ax.set_xticklabels(left)
    plt.savefig(fig_dir + "summary/num_clone.png", dpi=300)
    plt.close(1)


def summarize_RR(ABS_path, exp_name, tab_dir, sample_names, best_j):
    RR_mut = []; RR_vaf = []
    for sample_name in sample_names:
        lines = open(tab_dir + sample_name + "/stats.tsv", "r").readlines()
        RR_mut.append(float(lines[1].split("\t")[1]))
        RR_vaf.append(float(lines[2].split("\t")[1]))
    out = open(tab_dir + "summary/summary.tsv", "w")
    out.write("Sample_name\tThe_number_of_clones\tReconstruction_rate_for_mutational_type\t" +\
              "Reconstruction_rate_for_VAF\n")
    out.write("Mean\t" + str(np.mean(best_j)) + "\t" + str(np.mean(RR_mut)) + "\t" + str(np.mean(RR_vaf)) + "\n")
    for s,sample_name in enumerate(sample_names):
        out.write(sample_name + "\t" + str(best_j[s]) + "\t" + str(RR_mut[s]) + "\t" + str(RR_vaf[s]) + "\n")
    out.close()

