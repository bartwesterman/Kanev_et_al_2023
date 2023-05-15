import os, sys
import pandas as pd

# 328943 db_sets.csv
root = os.path.dirname(os.getcwd())
chembl = pd.read_csv(root + "/data/chembl/prep_chembl_v31.csv",header=None)

# 447518 db_sets.csv
with open (root + "/data/vanwesten/prep_vanwesten.csv", "r") as f_in:
    
    with open (root + "/data/db_sets.csv", "a") as f_out:
        for line in f_in:

            line = line.rstrip()
            items = line.split(",")

            kinase = items[1]
            inchi_key = items[4]
            experiment = kinase + "." + inchi_key

            check = chembl[chembl[0] == experiment].shape[0]
            if check == 0:
                f_out.write(line + "\n")