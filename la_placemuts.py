from numpy import log2, zeros, log10
import numpy as np
import pandas as pd

def enumerate_mutations(node, all_mutstr, all_cells):
    """Enumerate possible mutation placements on inferred tree

    Args:
        node (TreeNode): inferred subtree root
        all_mutstr ([pandas.Series]): all possible genotypes
        all_cells (set[string]): all cells in global inferred tree
    """
    if not node or len(node.cells) == 0:
        return
    curr_mut = pd.Series(zeros(len(all_cells)), all_cells, dtype=np.int32)
    
    # set mut to 1 for all cells of that node
    curr_mut.loc[list(node.cells)] = 1
    
    # add to all_str
    all_mutstr.append(curr_mut)
    enumerate_mutations(node.left_child, all_mutstr, all_cells)
    enumerate_mutations(node.right_child, all_mutstr, all_cells)

def compute_loglikelihood(place, actual, conf_mat):
    """Computes log-likelihood score for a potential mutation placement

    Args:
        place (pandas.Series): new placement genotype
        actual (pandas.Series): current genotype
        conf_mat (np.ndarray): confusion matrix with false positive/negative rates

    Returns:
        float: log-likelihood of placement
    """
    like = 0
    for d_i, e_i in zip(place, actual):
        like += log10(conf_mat[d_i, e_i])
    return like

def correct_geno_matrix(geno_mat, root_node, conf_mat):
    """Places mutations according to greatest log-likelihood score

    Args:
        geno_mat (pandas.DataFrame): current inferred genotype matrix
        root_node (TreeNode): node of inferred lineage tree
        conf_mat (np.ndarray): confusion matrix with false positive/negative rates
    """
    # build possible strings
    all_place = []
    all_cells = geno_mat.index.values
    enumerate_mutations(root_node, all_place, all_cells)
    
    # for each mutation
    for curr_mut in geno_mat.columns.values:
        orig_place = geno_mat.loc[:, curr_mut]
        
        best_like = compute_loglikelihood(all_place[0], orig_place, conf_mat)
        best_place = all_place[0]
                
        for placement in all_place[1:]:
            curr_like = compute_loglikelihood(placement, orig_place, conf_mat)
            if curr_like >  best_like:
                best_like = curr_like
                best_place = placement
        geno_mat.loc[:, curr_mut] = best_place