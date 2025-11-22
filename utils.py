class TreeNode():
    """represents internal tree node
    """
    def __init__(self, cells, left_child, right_child):
        """constructes TreeNode object

        Args:
            cells ([string]): cell list for node
            left_child (TreeNode): node's left subtree
            right_child (TreeNode): node's right subtree
        """
        self.cells = set(cells) 
        self.left_child = left_child
        self.right_child = right_child
        
def build_newick_string(current_node):
    """Build newick string representation for inferred tree

    Args:
        current_node (TreeNode): root node of subtree

    Returns:
        string: newick string representation for inferred tree
    """
    no_left_child = not current_node.left_child or len(current_node.left_child.cells) == 0
    no_right_child = not current_node.right_child or len(current_node.right_child.cells) == 0
    
    if no_left_child and no_right_child:
        base_str = ""
        for cell in current_node.cells:
            base_str += f"{cell}"
        return base_str
    else:
        return ("(" + build_newick_string(current_node.left_child) + ","
                    + build_newick_string(current_node.right_child) + ")")
    
def convert_to_probs(mut, fp_rate, fn_rate):
    """converts genotype entry to posterier probabilities for scistree2

    Args:
        mut (int): mutation present or not
        fp_rate (float): false positive rate
        fn_rate (float): false negative rate

    Returns:
        float: posterier probability
    """
    if mut==0:
        return (1-fp_rate)/(1-fp_rate+fn_rate)
    else:
        return (fp_rate)/(fp_rate+1-fn_rate)