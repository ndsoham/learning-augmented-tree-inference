## Learning-Augmented Evolutionary Tree Inference Algorithm
Contains sample implementation for learning-augmented algorithm.

### Running

_Note_: see requirements.txt for necessary prerequisites

`usage: run.py [-h] n m fp fn opt_tree opt_mat`

Run learning-augmented algorithmfor inferring evolutionary history

positional arguments:<br>
-  n:           cell count
-  m:           mutation count
-  fp:          false positive rate
-  fn:          false negative rate
-  opt_tree:    tree render file
-  opt_mat:     noise corrected matrix file

_Note_: `sample-data/` contains example input genotype matrices for $n=20$ and $n=50$ cells. `sample-oracles/` contains oracles used for inference for the sample data. The following commands can be used to test the implementation:

- `python run.py 20 20 0.01 0.05 test-n_20.png test-n_20.tsv`
- `python run.py 50 500 0.0001 0.2 test-n_50.png test-n_500.tsv`
