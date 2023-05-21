import os
import numpy as np
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