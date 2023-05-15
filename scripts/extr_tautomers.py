import sys

input_file_path = sys.argv[1]
input_file = input_file_path.split("/")[-1]
output_file_path = sys.argv[2]
header = True

mol_index = 0
mol_name = ""
mol_smile = ""

tautomer_index = 0
tautomer_inchi_key = ""
tautomer_smile = ""

input_header = []

index = 0
tautomers = {}
with open(input_file_path, "r") as f_in:

    with open(output_file_path, "w") as f_out:

        for line in f_in:

            line = line.rstrip()
            items = line.split(",")
            if header == True:
                input_header = items
                input_header = [w.replace('chembl_id', 'NAME') for w in input_header]
                input_header = [w.replace('name', 'NAME') for w in input_header]

            if header == False:
                
                if items[input_header.index("TAUTOMER_RANK")] == "Original structure":
                    mol_smile = items[input_header.index("SMILES")]
                    mol_index = int(items[input_header.index("MOLECULE_NO")])
                    mol_name = items[input_header.index("NAME")]

                if items[input_header.index("TAUTOMER_RANK")] != "Original structure":
                    tautomer_smile = items[input_header.index("SMILES")]
                    tautomer_inchi_key = items[input_header.index("InChIKey")]
                    tautomer_index = int(items[input_header.index("TAUTOMER_OF_MOLECULE_NO")])

                    if mol_index == tautomer_index and (mol_index + tautomer_index > 0):
                        f_out.write(input_file +"\t"+ str(tautomer_index) +"\t"+ mol_name +"\t"+ tautomer_inchi_key +"\t"+ tautomer_smile +"\t"+ mol_smile + "\n")
                        

                index += 1

            header = False