import os
import pandas as pd

root = os.path.dirname(os.getcwd())
D = pd.read_csv(root + "/data/data", 
                names=["epxeriment", "uniprot", "kinase", "group", "inchi_key", "value", "smiles", "mw"])

cpd_count = D.groupby(["inchi_key"])["inchi_key"].count()
cpds = cpd_count[cpd_count >=60].index.values
print(len(cpds))

test = D[D["inchi_key"].isin(cpds)]
train = D[~D["inchi_key"].isin(cpds)]

test.to_csv(root + "/data/dense_test.csv", index=False)
train.to_csv(root + "/data/dense_train.csv", index=False) 

print(test.shape)
print(train.shape)