import os, sys
from rdkit import Chem
from rdkit.Chem import Descriptors

root = os.path.dirname(os.getcwd())
with open(root + "/data/db_sets.csv", "r") as f_in:

    with open(root + "/data/db_sets_mw.csv", "w") as f_out:

        for line in f_in:

            line = line.rstrip()
            items = line.split(",")

            smile = items[-1]

            try:
                mol = Chem.MolFromSmiles(smile)
                mol_weight = Descriptors.MolWt(mol)
                if mol_weight >= 180 and mol_weight <= 700:
                    f_out.write(line + "," + str(mol_weight) + "\n")
            except:
                continue

