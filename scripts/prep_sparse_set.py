import os
import pandas as pd
from sklearn.model_selection import train_test_split

root = os.path.dirname(os.getcwd())
D = pd.read_csv(root + "/data/data_ml.csv")
print(D.shape)

train, test = train_test_split(D, stratify=D["kinase"], test_size=0.3, random_state=42)

cpds_train_set = train["inchi_key"].values

to_move = test[test['inchi_key'].isin(cpds_train_set)]

print(test.shape[0], to_move.shape[0], test.shape[0]-to_move.shape[0])

train = pd.concat([train, to_move])
test = test.drop(to_move.index)
print(train.shape)
print(test.shape)

test.to_csv(root + "/data/sparse_test.csv", index=False)
train.to_csv(root + "/data/sparse_train.csv", index=False) 