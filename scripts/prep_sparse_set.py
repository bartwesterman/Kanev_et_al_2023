import pandas as pd
from sklearn.model_selection import train_test_split

root = os.path.dirname(os.getcwd())
D = pd.read_csv(root + "/data/data.csv", 
                names=["epxeriment", "uniprot", "kinase", "group", "inchi_key", "value", "smiles", "mw"])

train, test = train_test_split(D, stratify=D["kinase"], test_size=0.15)

test.to_csv(root + "/data/sparse_test.csv", index=False)
train.to_csv(root + "/data/sparse_train.csv", index=False) 