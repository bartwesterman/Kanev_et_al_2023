import os, sys, math, pickle
import pandas as pd

root = os.path.dirname(os.getcwd())
db_file = root + "/data/vanwesten/si4_complete_dataset.txt"
cpds = root + "/data/vanwesten/vanwesten_tautomers.csv"
kinases = root + "/data/kinases.csv"

db_header = [   "Organism", "HGNC", "Mutation", "Is_Mutant", "HGNC_Mutant", "Compound", "HGNC_Mutant_Compound", "ATP_Concentration_uM", "_External_Organization", 
                "Method_Name", "Technology", "Target_Source", "Method_Class", "Method_Type", "Measured_Effect", "Source", "Original_Data_Type", "Original_Data_Type_General", 
                "Compound_Concentration_uM", "Original_Value", "Act_nM", "Act_Type", "Act_Type_General", "pAct", "Article_PMID", "Article_DOI", "Journal", "Article_Publication_Date", 
                "Article_Publication_Date_Formatted", "Article_Authors", "Article_Title", "Article_Reference", "Article_Abstract", "Target_Binding_Site_QualifierOriginal_Target_Sequence_Info", 
                "Target_Sequence_Domain", "Target_Sequence_Start_End", "Target_Sequence_Modification", "Catalytic_Status"]

cpds = pd.read_csv(cpds, delimiter="\t")
kinases = pd.read_csv(kinases)

single_point_inactive = []
dose_response = {}

header = True
with open(db_file, "r") as f_in:

    for line in f_in:

        line = line.rstrip()
        items = line.split("\t")

        if header == False:

            compound = items[db_header.index("Compound")]
            is_mutant = items[db_header.index("Is_Mutant")] 

            if is_mutant == "No":

                if cpds[cpds["chembl_id"] == compound].shape[0] > 0:

                    concentration = items[db_header.index("Compound_Concentration_uM")] # 10 or NA
                    pact = items[db_header.index("pAct")]

                    std_type = items[db_header.index("Original_Data_Type_General")]
                    std_value = items[db_header.index("Original_Value")]

                    if std_value == ">10000":
                        std_value = "1000"

                    if ">" not in std_value and "<" not in std_value:

                        if std_value == "No Effect":
                            std_value = 100

                        std_value = float(std_value)

                        hgnc = items[db_header.index("HGNC")]
                        inchi_key = cpds[cpds["chembl_id"] == compound]["inchi_key"].values[0]
                        
                        experiment = hgnc + "." + inchi_key

                        # insert only inactive
                        if std_type == "Pct_Ctrl" and concentration == 10:

                            # Pct_Ctrl = 100 × Mean_on_compound/Mean_on_controls
                            if std_value >= 90:

                                if experiment not in single_point_inactive:
                                    single_point_inactive.append(experiment)

                        # insert only inactive
                        elif std_type == "Pct_Inhib" and concentration == 10:

                            # Pct_Inhib is defined as 1 – (100 × Mean_on_sample/Mean_on_controls)
                            if std_value <= 10:

                                if experiment not in single_point_inactive:
                                    single_point_inactive.append(experiment)
                            
                        # insert both active and inactive
                        else:

                            if pact:

                                log_activity = float(pact)

                                if log_activity < 5:
                                    log_activity = 5

                                if experiment not in dose_response.keys():
                                    dose_response[experiment] = []

                                dose_response[experiment].append(log_activity)
        header = False


chembl = {}
# calc median of dose_response
for experiment in dose_response.keys():
    activities = dose_response[experiment]
    activity = sorted(activities)[int(len(activities)/2)]
    chembl[experiment] = activity


# if single point experiment not in dose_response --> insert
for experiment in single_point_inactive:
    if experiment not in chembl.keys():
        chembl[experiment] = 5

print(len(chembl))

# print to file
with open(root + "/data/vanwesten/prep_vanwesten.csv", "w") as f_out:
    for experiment, value in chembl.items():
        
        inchi_key = experiment.split(".")[1]
        smiles = cpds[cpds["inchi_key"] == inchi_key]["tautomer_smiles"].values[0]

        kinase = experiment.split(".")[0]

        check_kinase_hgnc = kinases[kinases["hgnc_mgi"] == kinase].shape[0]
        check_kinase_name = kinases[kinases["kinase_name"] == kinase].shape[0]
        if check_kinase_hgnc:
            accession = kinases[kinases["hgnc_mgi"] == kinase]["accession"].values[0]
            group = kinases[kinases["hgnc_mgi"] == kinase]["kinase_group"].values[0]
        elif check_kinase_name:
            accession = kinases[kinases["kinase_name"] == kinase]["accession"].values[0]
            group = kinases[kinases["kinase_name"] == kinase]["kinase_group"].values[0]    
        else:
            skip=1

        f_out.write(experiment + "," + accession + "," + kinase + "," + group + "," + inchi_key + "," + str(value) + "," + smiles + "\n")
