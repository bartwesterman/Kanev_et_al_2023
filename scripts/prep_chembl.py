import os, sys, math, pickle
import pandas as pd

root = os.path.dirname(os.getcwd())
db_file = root + "/data/chembl/chembl_v31.csv"
cpds = root + "/data/chembl/chembl_tautomers.csv"
kinases = root + "/data/kinases.csv"

db_header = ["kinase_name","hgnc_mgi","accession","kinase_group","chembl_id","standard_inchi_key","canonical_smiles",
             "standard_type","standard_relation","standard_value","standard_units","potential_duplicate","pchembl_value",
             "activity_id","assay_id","doi"]

cpds = pd.read_csv(cpds, delimiter="\t")
kinases = pd.read_csv(kinases)

single_point_inactive = []
dose_response = {}

header = True
with open(db_file, "r") as f_in:

    
    for line in f_in:

        line = line.rstrip()
        items = line.split(",")

        if header == False:

            chembl_id = items[db_header.index("chembl_id")]

            if cpds[cpds["chembl_id"] == chembl_id].shape[0] > 0:

                std_type = items[db_header.index("standard_type")]
                std_relation = items[db_header.index("standard_relation")]
                std_value = float(items[db_header.index("standard_value")])

                accession = items[db_header.index("accession")]
                inchi_key = cpds[cpds["chembl_id"] == chembl_id]["inchi_key"].values[0]
                
                experiment = accession + "." + inchi_key

                # insert only inactive
                if std_type == "Activity":

                    # Pct_Ctrl = 100 × Mean_on_compound/Mean_on_controls
                    if (std_relation == ">" or std_relation == "=") and std_value >= 90:

                        if experiment not in single_point_inactive:
                            single_point_inactive.append(experiment)

                # insert only inactive
                elif std_type == "Inhibition":

                    # Pct_Inhib is defined as 1 – (100 × Mean_on_sample/Mean_on_controls)
                    if (std_relation == "<" or  std_relation == "=") and std_value <= 10:

                        if experiment not in single_point_inactive:
                            single_point_inactive.append(experiment)
                    
                # insert both active and inactive
                else:

                    if not std_relation == "<":

                        if not (std_relation == ">" and std_value < 10000):

                            if experiment not in dose_response.keys():
                                dose_response[experiment] = []

                            if (std_relation == ">" or  std_relation == "=") and std_value >= 10000:
                                dose_response[experiment].append(5)

                            elif (std_relation == "=") and std_value < 10000:

                                try:
                                    log_activity = round(-1*math.log10(std_value*10**-9),3)

                                    if log_activity < 5:
                                        log_activity = 5

                                    dose_response[experiment].append(log_activity)
                                except:
                                    print("log_activity exception occurred: " + str(std_value))
        header = False


chembl = {}
# calc median of dose_response
for experiment in dose_response.keys():
    activities = dose_response[experiment]
    if len(activities) > 0:
        activity = sorted(activities)[int(len(activities)/2)]
        chembl[experiment] = activity


# if single point experiment not in dose_response --> insert
for experiment in single_point_inactive:
    if experiment not in chembl.keys():
        chembl[experiment] = 5

# print to file
with open(root + "/data/chembl/prep_chembl_v31.csv", "w") as f_out:
    for experiment, value in chembl.items():
        
        inchi_key = experiment.split(".")[1]
        smiles = cpds[cpds["inchi_key"] == inchi_key]["tautomer_smiles"].values[0]

        accession = experiment.split(".")[0]
        kinase = kinases[kinases["accession"] == accession]["hgnc_mgi"].values[0]
        group = kinases[kinases["accession"] == accession]["kinase_group"].values[0]

        f_out.write(experiment + "," + accession + "," + kinase + "," + group + "," + inchi_key + "," + str(value) + "," + smiles + "\n")