import os
ABS_path = os.getcwd() + "/" #FIXME according to your environment

import sys
sys.path.append(ABS_path + "scripts/")
sys.path.append(ABS_path + "drawing/")
from multiprocessing import Pool
import subprocess
import utility as util
import select_best as sb
import test_result as tr

def main():
    result_path = ABS_path + "result/" + exp_name
    for sample_name in sample_names:
        os.makedirs(result_path + "/" + sample_name, exist_ok=True)
    
    best_j = sb.select_best_j(ABS_path, exp_name, max_J, max_run, min_J = min_J)
    best_r = sb.select_best_r(ABS_path, exp_name, best_j, max_run_pe)
    tr.test(ABS_path, exp_name, best_j, best_r, ref_names)

    print("--------------")
    print("Complete test.")
    print("--------------\n")


if __name__ == '__main__':
    args = sys.argv
    exp_name = args[1]
    min_J = 1
    max_J = 4
    max_run = 3
    max_run_pe = 3

    sample_names, purity = util.load_sample_names_and_purity(ABS_path, exp_name)
    ref_names = util.load_ref_sig(ABS_path, exp_name)

    main()
