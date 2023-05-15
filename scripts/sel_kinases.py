import os
import pandas as pd

root = os.path.dirname(os.getcwd())
kinase_mapping = pd.read_csv(root + "/data/kinases.csv")

sel_kinases = [ "LCK", "KDR", "EGFR", "SRC", "MAPK14", "FGFR1", "MAP2K1", "CDK2", "GSK3B", "MAPK10", "ABL1", "CHEK1", "ITK", "PDPK1", "NTRK1", "SYK", "INSR", "IGF1R", 
                "MAPK1", "HCK", "PTK2", "MAPKAPK2", "IRAK4", "MAPK8", "AKT1", "FGFR2", "KIT", "MET", "PAK4", "RET", "TGFBR1", "BTK", "CHEK2", "PIK3CA", "WEE1", "CSNK1D", 
                "CSNK2A1", "BRAF", "ROCK1", "RPS6KB1", "DYRK1A", "AURKA", "EPHA3", "TTK", "PIM1", "ALK", "MAP2K7", "PIK3CG", "CDK9", "NEK2", "EPHA2", "PAK1", "PRKACA", 
                "CLK1", "MAP3K7", "MERTK", "DAPK1", "CSNK2A2", "MAP3K5", "FGFR4", "DDR1", "DYRK2", "RIPK2", "TBK1", "CDK8", "MELK", "ACVR1", "ERN1", "STK24"]


db_sets_mw = pd.read_csv(root + "/data/db_sets_mw.csv", header=None)

data = db_sets_mw[db_sets_mw[2].isin(sel_kinases)]
data.to_csv(root + "/data/data.csv", header=None, index=None)
print(data.shape)