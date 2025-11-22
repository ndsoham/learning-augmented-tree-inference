from la_oracle import xgb_oracle_predict
from numpy import log2
from utils import TreeNode
from la_bruteforce import run_scis2
    
def counterpart_tester_strong(allowed_cells, u, v, thresh1, thresh2, cell_df, m_i, n_i, fp_i, fn_i):
    """Tests if u is a counterpart of v

    Args:
        allowed_cells ([string]): which cells can be tested
        u (string): baseline vertex
        v (string): test vertex
        thresh1 (float): threshold for counter 1
        thresh2 (float): threshold for counter 2
        cell_df (pandas.DataFrame): genotype matrix dataframe
        m_i (int): mutation count
        n_i (int): cell count
        fp_i (float): false positive rate
        fn_i (float): false negative rate
    Returns:
        bool: whether v is a counterpart of u
    """
    c1 = 0
    c2 = 0
    
    test_cells = allowed_cells - {u, v}
    
    for t_i in test_cells:
        split_cell = xgb_oracle_predict([u, v, t_i], cell_df, m_i, n_i, fp_i, fn_i)
        
        if (split_cell == v):
            c1+=1
        elif (split_cell == t_i):
            c2+=1
            
    if c1 <= thresh1 and c2 <= thresh2:
        return True
    return False


def strong_partial_tree(cell_matrix_to_mutate, node, m_i, n_i, fp_i, fn_i, rem_muts, epsilon=0.01):
    """Constructs strong partial tree
    Args:
        cell_matrix_to_mutate (pd.DataFrame): cell genotype matrix
        node (TreeNode): treenode to split
        m_i (int): mutation count
        n_i (int): cell count
        fp_i (int): false positive rate
        fn_i (int): false negative rate
        rem_muts set[string]: mutations left to place
        epsilon (float): sufficiently small constant to control counterpart threshold
    """
    
    sp_tree_array = dict()
    
    # if too few cells left, use bruteforce
    if len(node.cells) < log2(n_i):
        run_scis2(cell_matrix_to_mutate, node, rem_muts, fp_i, fn_i)
        return

    # find the strongest counterpart
    for base_cell in node.cells:
        sp_tree = TreeNode([], None, None)
        
        sample_verts = (node.cells - set({base_cell}))
        
        thresh1 = (3/5 - epsilon) * len(sample_verts)
        thresh2 = (1/6 - epsilon) * len(sample_verts)
        
        for test_cell in (node.cells - set({base_cell})):
            is_counterpart = counterpart_tester_strong(node.cells, base_cell, test_cell, 
                                                       thresh1, thresh2, cell_matrix_to_mutate, 
                                                       m_i, n_i, fp_i, fn_i)
            if is_counterpart:
                sp_tree.cells.add(test_cell)
                
        sp_tree_array[base_cell] = sp_tree
        
    largest_size = lambda base_cell: len(sp_tree_array[base_cell].cells)
    best_split_cell = max(sp_tree_array, key=largest_size)
    left_tree = sp_tree_array[best_split_cell]
    
    # if not strong counterpart was found, use bruteforce
    if len(left_tree.cells) == 0:
        run_scis2(cell_matrix_to_mutate, node, rem_muts, fp_i, fn_i)
        return
    
    right_cells = (node.cells - left_tree.cells)
    right_tree = TreeNode(right_cells, None, None)
    
    node.left_child = left_tree
    node.right_child = right_tree
    
    # recursively construct left and right child
    strong_partial_tree(cell_matrix_to_mutate, left_tree, m_i, n_i, fp_i, fn_i, rem_muts)
    strong_partial_tree(cell_matrix_to_mutate, right_tree, m_i, n_i, fp_i, fn_i, rem_muts)

