import os, pickle
import numpy as np
import pandas as pd

from common_functions import get_indeces_identical_mols, calc_similarity_matrix


root = os.path.dirname(os.getcwd())

D = pd.read_csv(root + "/data/data.csv")
cpds_inchikey_smiles = D[["inchi_key", "smiles"]].drop_duplicates()
print("data =", D.shape)

# remove identical compounds
if os.path.exists('data_cpd_similarity_matrix.npy'):
    similarity_matrix = np.load('data_cpd_similarity_matrix.npy')
else:
    similarity_matrix = calc_similarity_matrix(cpds_inchikey_smiles["smiles"].values)

# [61, 234, 682, 733 ...]
rm_smile_indeces = get_indeces_identical_mols(similarity_matrix) 

rm_smiles = [cpds_inchikey_smiles["smiles"].values[i] for i in rm_smile_indeces]

to_rm_mols = D[D['smiles'].isin(rm_smiles)].index

D = D.drop(to_rm_mols)

# remove duplicated smiles and inchi keys
cpds_inchikey_smiles = D[["inchi_key", "smiles"]].drop_duplicates()

duplicated_inchikeys = cpds_inchikey_smiles[cpds_inchikey_smiles['inchi_key'].duplicated()]
duplicated_smiles = cpds_inchikey_smiles[cpds_inchikey_smiles['smiles'].duplicated()]

X = D[~D['inchi_key'].isin(duplicated_inchikeys['inchi_key'].values)]
X = X[~X['smiles'].isin(duplicated_smiles['smiles'].values)]
X.to_csv(root + "/data/data_ml.csv", index=None)
count_experiments = X["epxeriment"].shape

cpd_count = D.groupby(["inchi_key"])["inchi_key"].count()

# generate the train and test sets
selected_cpds = cpd_count[cpd_count >=60].index.values
print("data =", X.shape)

test = X[X["inchi_key"].isin(selected_cpds)]
train = X[~X["inchi_key"].isin(selected_cpds)]

test.to_csv(root + "/data/dense_test.csv", index=False)
train.to_csv(root + "/data/dense_train.csv", index=False) 

print(test.shape)
print(train.shape)