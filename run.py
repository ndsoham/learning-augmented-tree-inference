from utils import TreeNode, build_newick_string
from la_augmented import strong_partial_tree
from la_placemuts import correct_geno_matrix
import pandas as pd
import numpy as np
from ete3 import Tree
import argparse

parser = argparse.ArgumentParser(description="Run learning-augmented algorithm"
                                             "for inferring evolutionary history")
parser.add_argument("n", help="cell count")
parser.add_argument("m", help="mutation count")
parser.add_argument("fp", help="false positive rate")
parser.add_argument("fn", help="false negative rate")
parser.add_argument("opt_tree", help="tree render file")
parser.add_argument("opt_mat", help="noise corrected matrix file")

args = parser.parse_args()

n_i=int(args.n); m_i=int(args.m); fp_i=float(args.fp); fn_i=float(args.fn)

geno_fp = f"sample-data/simNo_{1}-n_{n_i}-m_{m_i}-fp_{fp_i}-fn_{fn_i}-na_0.SC.after_noise"

noisy_genotype = pd.read_csv(geno_fp, sep="\t").set_index("cell_id_x_mut_id")

all_muts = noisy_genotype.columns.to_list()
all_cells = noisy_genotype.index.values
conf_mat = np.array([[1-fp_i, fp_i], [fn_i, 1-fn_i]])
root = TreeNode(all_cells, None, None)

strong_partial_tree(noisy_genotype, root, m_i, n_i, fp_i, fn_i, all_muts)
correct_geno_matrix(noisy_genotype, root, conf_mat)

inf_newick = build_newick_string(root)+";"

inferred_tree = Tree(inf_newick, format=9)
inferred_tree.show()
inferred_tree.render(f"output-trees/{args.opt_tree}")
noisy_genotype.to_csv(f"output-matrices/{args.opt_mat}", sep="\t")