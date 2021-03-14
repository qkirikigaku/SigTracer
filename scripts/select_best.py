import numpy as np
import utility as util

def select_best_j(ABS_path, exp_name, max_J, max_run, min_J = 1):
    sample_names,purity = util.load_sample_names_and_purity(ABS_path, exp_name)
    best_j = []
    for s, sample_name in enumerate(sample_names):
        ELBO = []; max_r = []
        for j in range(min_J, max_J+1):
            temp_elbo = []
            for r in range(1, max_run+1):
                temp_elbo.append(util.load_elbo(ABS_path, exp_name, sample_name,\
                                 j, r, "ms"))
            ELBO.append(max(temp_elbo))
            max_r.append(np.argmax(temp_elbo)+1)
        J_selected = np.argmax(ELBO) + min_J
        r_selected = max_r[np.argmax(ELBO)]
        modified_J = remove_minor_clone(ABS_path, exp_name, sample_name, J_selected, r_selected)
        best_j.append(modified_J)
    return best_j


def select_best_r(ABS_path, exp_name, best_j, max_run, tag=""):
    sample_names, purity = util.load_sample_names_and_purity(ABS_path, exp_name)
    best_r = []
    for s, sample_name in enumerate(sample_names):
        ELBO = []
        for r in range(1, max_run+1):
            ELBO.append(util.load_elbo(ABS_path, exp_name, sample_name,\
                        best_j[s], r, "pe", tag=tag))
        r_selected = np.argmax(ELBO) + 1
        best_r.append(r_selected)
    return best_r


def remove_minor_clone(ABS_path, exp_name, sample_name, J_selected, r_selected):
    pi = util.load_pi(ABS_path, exp_name, sample_name, J_selected, r_selected, "ms")
    num_to_remove = 0
    for temp_pi in pi:
        if(temp_pi < 0.01): num_to_remove += 1
    return J_selected - num_to_remove
