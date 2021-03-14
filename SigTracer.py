import os
ABS_path = os.getcwd() + "/" # FIXME according to your environment.
num_processes=8 # FIXME according to your environment.
Simulation_flag = True # FIXME to change "False" when applying an actual data.

import sys
sys.path.append(ABS_path + "scripts/")
sys.path.append(ABS_path + "drawing/")
from multiprocessing import Pool
import subprocess
import utility as util
import select_best as sb
import visualize
import evaluate_simulation as es

def main():
    result_path = ABS_path + "result/" + exp_name
    for sample_name in sample_names:
        os.makedirs(result_path + "/" + sample_name, exist_ok=True)
    
    arguments = []
    for sample_name in sample_names:
        for j in range(min_J, max_J + 1):
            for temp_r in range(1, max_run+1):
                arguments.append((exp_name, sample_name, ABS_path, j, temp_r,\
                              ref_names[sample_name], "ms"))
    pool = Pool(processes=num_processes)
    _ = pool.starmap(execute, arguments)
    best_j = sb.select_best_j(ABS_path, exp_name, max_J, max_run, min_J = min_J)
 
    print("------------------------------")
    print("Complete to select num_clones.")
    print("------------------------------\n")
    
    arguments = []
    for s,sample_name in enumerate(sample_names):
        for r in range(1, max_run_pe+1):
            arguments.append((exp_name, sample_name, ABS_path, best_j[s],\
                              r, ref_names[sample_name],"pe"))
    pool = Pool(processes=num_processes)
    _= pool.starmap(execute, arguments)
    
    print("--------------------------------")
    print("Complete to estimate parameters.")
    print("--------------------------------\n")

    best_r = sb.select_best_r(ABS_path, exp_name, best_j, max_run_pe)
    file_foots = ["BB", "activity", "alpha", "elbo", "pi", "qC", "qU", "qUC", "qz", "qZU", "signature"]
    for s, sample_name in enumerate(sample_names):
        best_path = ABS_path + "result/" + exp_name + "/" + sample_name + "/suggested_solution"
        os.makedirs(best_path, exist_ok=True)
        file_head = ABS_path + "result/" + exp_name + "/" + sample_name + "/J" + str(best_j[s])\
                    + "_run" + str(best_r[s]) + "_pe_"
        for file_foot in file_foots:
            file_name = file_head + file_foot + ".txt"
            cmd = "cp " + file_name + " " + ABS_path + "result/" + exp_name + "/" + sample_name\
                    + "/suggested_solution/J" + str(best_j[s]) + "_" + file_foot + ".txt"
            subprocess.call(cmd.split())


    visualize.visualize_result(ABS_path, exp_name, best_j, best_r, ref_names, num_processes)

    print("-----------------------")
    print("Complete visualization.")
    print("------------------------\n")
    
    if(Simulation_flag == True):
        lines = open(ABS_path + "data/" + exp_name + "/setting.txt").readlines()
        true_J = int(lines[3].split()[1])
        es.evaluate_model_selection(true_J, best_j, min_J, max_J, sample_names, result_path)
        run_with_true_J(exp_name, sample_names, true_J, ABS_path, ref_names, result_path)
        print("--------------------------------")
        print("Complete to evaluate simulation.")
        print("--------------------------------")


def execute(exp_name, sample_name, ABS_path, j, i, ref_sig, mode):
    cmd = ABS_path + "bin/ST -d " + exp_name\
        + " -s " + sample_name\
        + " -j " + str(j)\
        + " -i " + str(i)\
        + " -p " + ABS_path\
        + " -m " + mode
    for x in ref_sig:
        cmd += " -r " + x
    subprocess.call(cmd.split())


def run_with_true_J(exp_name, sample_name, true_J, ABS_path, ref_names, result_path):
    arguments = []
    for s, sample_name in enumerate(sample_names):
        for temp_r in range(1, max_run_pe+1):
            arguments.append((exp_name, sample_name, ABS_path, true_J, temp_r, ref_names[sample_name], "pe"))
    pool = Pool(processes=num_processes)
    _ = pool.starmap(execute, arguments)
    best_r = sb.select_best_r(ABS_path, exp_name, [true_J for s in range(len(sample_names))], max_run_pe)
    es.evaluate_accuracy(ABS_path, exp_name, true_J, best_r, ref_names, sample_names,
                         purity, result_path)


if __name__ == '__main__':
    args = sys.argv
    exp_name = args[1]
    min_J = 1
    max_J = 4
    max_run = 2
    max_run_pe = 2

    sample_names, purity = util.load_sample_names_and_purity(ABS_path, exp_name)
    ref_names = util.load_ref_sig(ABS_path, exp_name)

    main()
