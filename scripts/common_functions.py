import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

def calc_similarity_matrix(smiles_list):

    mols = [Chem.MolFromSmiles(smile) for smile in smiles_list]
    fgrp = [AllChem.GetMorganFingerprintAsBitVect(mol,2,nBits=1024) for mol in mols]
        
    similarities = np.zeros((len(fgrp), len(fgrp)), dtype=np.float16)

    for i in range(len(fgrp)):

        similarity = DataStructs.BulkTanimotoSimilarity(fgrp[i], fgrp[:i])
        similarities[i, :i] = similarity
        similarities[:i, i] = similarity

    if not os.path.exists('data_cpd_similarity_matrix.npy'):
        np.save('data_cpd_similarity_matrix.npy', similarities)

    return(similarities)


def get_indeces_identical_mols(similarities_matrix):

    max_vals = np.max(similarities_matrix, axis=0) # columns
    identical_mols = np.where(max_vals > 0.99)[0]
    indeces_identical_mols = identical_mols 

    for index in identical_mols:
        for i in np.where(similarities_matrix[:,index] == 1)[0]:
            if i not in indeces_identical_mols:
                indeces_identical_mols.append(i)

    return(indeces_identical_mols)


def morgan_fp(df):
    features = []
    for i, r in df.iterrows():
        smile = r["smiles"]
        mol = Chem.MolFromSmiles(smile)
        bits = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        features.append(list(bits))
    return(features)

def one_hot_encoding(df,onehotdict):
    onehot = np.zeros((df.shape[0], 69))
    row_index = 0
    for i, r in df.iterrows():
        onehot[row_index, onehotdict[r["kinase"]]] = 1
        row_index += 1

    return(onehot)

def zscale_sequence(df,features_path):
    kinases = df['kinase']
    zscales = pd.read_csv(features_path + "/zscales_seq.txt", delimiter="\t")
    zscales = zscales.drop(["Pocket"], axis=1)
    merged = pd.merge(left=kinases, right=zscales, how='left', left_on='kinase', right_on='HGNC/MGI')
    merged = merged.drop(['kinase'], axis=1)
    merged = merged.drop(['HGNC/MGI'], axis=1)
    return(np.array(merged))

def zscale_residue(df,features_path):
    kinases = df['kinase']
    zscales = pd.read_csv(features_path +  "/zscales_res.txt", delimiter="\t")
    zscales = zscales.drop(["PocketSequence"], axis=1)
    merged = pd.merge(left=kinases, right=zscales, how='left', left_on='kinase', right_on='HGNC/MGI')
    merged = merged.drop(['kinase'], axis=1)
    merged = merged.drop(['HGNC/MGI'], axis=1)
    return(np.array(merged))

def protvec(df,features_path):
    kinases = df['kinase']
    protvec = pd.read_csv(features_path + "/protvec.txt", delimiter=" ")
    protvec = protvec.drop(["pocket"], axis=1)
    merged = pd.merge(left=kinases, right=protvec, how='left', left_on='kinase', right_on='HGNC/MGI')
    merged = merged.drop(['kinase'], axis=1)
    merged = merged.drop(['HGNC/MGI'], axis=1)
    return(np.array(merged))

def cnn(df,features_path):
    kinases = df['kinase']
    cnn = pd.read_csv(features_path + "/model_3_127/fp_2.txt", delimiter=" ")
    merged = pd.merge(left=kinases, right=cnn, how='left', left_on='kinase', right_on='HGNC/MGI')
    merged = merged.drop(['kinase'], axis=1)
    merged = merged.drop(['HGNC/MGI'], axis=1)
    return(np.array(merged))