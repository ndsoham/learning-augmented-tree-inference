from utils import TreeNode, convert_to_probs
from functools import partial
from ete3 import Tree
import scistree2 as s2

def construct_tree_recursively(node, ete3_tree):
    """subroutine to attached inferred scistree2 subtree to current node

    Args:
        node (TreeNode): root of tree to resolve
        ete3 (ete3.TreeNode): scistree2 resolved tree
    """
    if len(node.cells) == 1:
        return
    left_ete3 = ete3_tree.children[0]
    right_ete3 = ete3_tree.children[1]
    # assign left child
    node.left_child = TreeNode(left_ete3.get_leaf_names(), None, None)
    # assign right child
    node.right_child = TreeNode(right_ete3.get_leaf_names(), None, None)
    # recur left
    construct_tree_recursively(node.left_child, left_ete3)
    # recur right
    construct_tree_recursively(node.right_child, right_ete3)
    

def run_scis2 (cell_matrix_to_mutate, node, rem_muts, fp_rate, fn_rate):
    """runs scistree2 to resolve small subtree

    Args:
        cell_matrix_to_mutate (pd.DataFrame): cell genotype matrix
        node (TreeNode): root of tree to resolve
        rem_muts ([string]): remaining mutations to place
        fp_rate (float): false positive rate
        fn_rate (float): false negative rate
    """
    
    # only give scistree2 probably coordinates
    cells_to_use = list(node.cells)
    all_sums = cell_matrix_to_mutate.loc[cells_to_use, rem_muts].sum(axis="index")
    muts_to_use = all_sums.loc[all_sums > 0].index.values
    place_muts = True
    
    if len(muts_to_use) == 0 and len(node.cells) > 1:
        # if no mutations left, relax condiation and allow all muts
        muts_to_use = cell_matrix_to_mutate.columns.values
        place_muts = False  
    elif len(muts_to_use) == 0 or len(node.cells) == 1:
        return
    
    # create the prob matrix
    prob_matrix = partial(convert_to_probs, fp_rate=fp_rate, fn_rate=fn_rate)
    cells_df = cell_matrix_to_mutate.loc[cells_to_use, muts_to_use]
    cells_df = cells_df.map(prob_matrix)
    probs = cells_df.values.transpose()
    
    cell_names = cells_df.index.to_list()
    mut_names = cells_df.columns.to_list()
    
    # submit prob matrix to s2
    gp_matrix = s2.probability.from_probs(probs, cell_names, mut_names)
    caller = s2.ScisTree2(threads=10)
    s2_tree, imputed_genotype, _ = caller.infer(gp_matrix)
    ete_tree = Tree(s2_tree.output(), format=9)
    
    if place_muts:
        # change cell_matrix_to_mutate
        cell_matrix_to_mutate.loc[cells_to_use, muts_to_use] = imputed_genotype.transpose()
        
        for mut in muts_to_use:
            rem_muts.remove(mut)

    # recursively construct subtree
    construct_tree_recursively(node, ete_tree)
    