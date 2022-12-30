import os
import glob
import pickle
import pandas as pd
import numpy as np
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

# read text file into pandas DataFrame
ex_matrix = pd.read_csv("1.1_exprMatrix_filtered_t.txt", sep='\t', header=0, index_col=0)
ex_matrix.shape
tf_names = load_tf_names("1.1_inputTFs.txt")
adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)
adjacencies.to_csv("adjacencies.tsv", sep='\t', index=False, header=False)
