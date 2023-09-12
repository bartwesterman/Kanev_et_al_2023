import pickle, os, sys
import pandas as pd
import numpy as np
from os import path
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import accuracy_score
from sklearn.metrics import mean_squared_error 
from sklearn.preprocessing import OneHotEncoder

from common_functions import morgan_fp, one_hot_encoding, zscale_sequence, zscale_residue, protvec, cnn

cwd = os.getcwd()
data_path = os.path.dirname(cwd) + "/data"
features_path = os.path.dirname(cwd) + "/features"

full_data_path = data_path + "/data_ml.csv"
#train_path = data_path + "/dense_train.csv"
#test_path = data_path + "/dense_test.csv"
train_path = data_path + "/sparse_train.csv"
test_path = data_path + "/sparse_test.csv"

full_data = pd.read_csv(full_data_path)
train = pd.read_csv(train_path)
test = pd.read_csv(test_path)

index = 0
onehotdict = {}
for kinase in full_data["kinase"].unique():
    onehotdict[kinase] = index
    index += 1

X_train_inhibitor = np.array(morgan_fp(train))
X_test_inhibitor = np.array(morgan_fp(test))

y_train = train["value"].values
y_test = test["value"].values

for model in ["model_1", "model_2", "model_3", "model_4", "model_5", "model_6"]:

    info = ""
    if model == "model_1":
        info = "morganfp"
    elif model == "model_2":
        info = "onehot_morganfp"
    elif model == "model_3":
        info = "zscale_seq_morganfp"
    elif model == "model_4":
        info = "zscale_res_morganfp"
    elif model == "model_5":
        info = "protvec_morganfp"
    elif model == "model_6":
        info = "cnn_morganfp"
    else:
        print('model not recognized ...', file=sys.stderr)

    if "dense" in train_path:
        info += "_dense"
        print("training " + model + " on dense data")
    else:
        info += "_sparse"
        print("training " + model + " on sparse data")

    X_train = None 
    X_test = None 
    X_train_protein = None 
    X_test_protein = None 

    reg = None

    print("preparing data ...")
    if model == "model_1":
        X_train = X_train_inhibitor
        X_test = X_test_inhibitor
    if model == "model_2":
        X_train_protein = one_hot_encoding(train, onehotdict)
        X_test_protein = one_hot_encoding(test, onehotdict)
        if X_train_protein.shape[0] == X_train_inhibitor.shape[0]:
            X_train = np.hstack((X_train_protein, X_train_inhibitor))
            X_test = np.hstack((X_test_protein, X_test_inhibitor))
        else:
            print('shape input features does not match ...', file=sys.stderr)
    elif model == "model_3":
        X_train_protein = zscale_sequence(train,features_path)
        X_test_protein = zscale_sequence(test,features_path)  
        if X_train_protein.shape[0] == X_train_inhibitor.shape[0]:
            X_train = np.hstack((X_train_protein, X_train_inhibitor))
            X_test = np.hstack((X_test_protein, X_test_inhibitor))
        else:
            print('shape input features does not match ...', file=sys.stderr)
    elif model == "model_4":
        X_train_protein = zscale_residue(train,features_path)
        X_test_protein = zscale_residue(test,features_path)
        if X_train_protein.shape[0] == X_train_inhibitor.shape[0]:
            X_train = np.hstack((X_train_protein, X_train_inhibitor))
            X_test = np.hstack((X_test_protein, X_test_inhibitor))
        else:
            print('shape input features does not match ...', file=sys.stderr)
    elif model == "model_5":
        X_train_protein = protvec(train,features_path)
        X_test_protein = protvec(test,features_path)
        if X_train_protein.shape[0] == X_train_inhibitor.shape[0]:
            X_train = np.hstack((X_train_protein, X_train_inhibitor))
            X_test = np.hstack((X_test_protein, X_test_inhibitor))
        else:
            print('shape input features does not match ...', file=sys.stderr)
    elif model == "model_6":
        X_train_protein = cnn(train,features_path)
        X_test_protein = cnn(test,features_path) 
        if X_train_protein.shape[0] == X_train_inhibitor.shape[0]:
            X_train = np.hstack((X_train_protein, X_train_inhibitor))
            X_test = np.hstack((X_test_protein, X_test_inhibitor))
        else:
            print('shape input features does not match ...', file=sys.stderr)

    print(X_train.shape)
    print(X_test.shape)

    print("training ...")
    reg = RandomForestRegressor(n_estimators=100, n_jobs=-1, random_state=42)
    reg = reg.fit(X_train, y_train)

    print("saving model ...")
    sm_feat = pickle.dump(reg, open("evaluate_rf_models/saved_models/" + model +"_"+info+".pl","wb" ))

    y_pred = reg.predict(X_test)
    with open("evaluate_rf_models/saved_models/" + model + "_" + info+'_pred.npy', 'wb') as f:
        np.save(f, y_pred)

    with open("evaluate_rf_models/saved_models/" + model + "_" + info + ".txt", "w") as f_out:
        f_out.write("model = " + str(model) + "\n")
        f_out.write("info = " + str(info) + "\n")
        f_out.write("train size = " + str(X_train.shape) + "\n")
        f_out.write("test size = " + str(X_test.shape) + "\n")
        f_out.write("RMSE = "+ str(mean_squared_error(y_test, y_pred, squared=False)) + "\n")

    print()
