### Required python packages:
GE  : pandas, numpy\
CH  : rdkit\
ML  : scikit-learn, pytorch, [libmolgrid](https://github.com/gnina/libmolgrid)\
VIS : seaborn, matplotlib

### Workflow:
:arrow_right: Prepare tautomers (can be skipped):\
(1) bash gen_tautomers.sh \
(2) python extr_tautomers.py

:arrow_right: Unzip (before running steps 3 to 9)\
../data/chembl/chembl_tautomers.csv.zip\
../data/chembl/chembl_v31.csv.zip\
../data/vanwesten/si4_complete_dataset.txt.zip

:arrow_right: Prepare the data\
(3) python prep_chembl.py\
(4) python prep_vanwesten.py\
(5) python merge_db_sets.py\
(6) python rm_cpds.py\
(7) python sel_kinases.py\
(8) python prep_dense_set.py\
(9) python prep_sparse_set.py

:arrow_right: Prepare the 3D-CNN model:

:arrow_right: Analysis:

### Standalone software:
#### Ambit-Tautomer: An Open Source Tool for Tautomer Generation
Nikolay T Kochev, Vesselina H Paskaleva, Nina Jeliazkova\
Mol Inform. 2013 Jun;32(5-6):481-504.\
PMID: [27481667](https://pubmed.ncbi.nlm.nih.gov/27481667/)\
DOI: [10.1002/minf.201200133](https://onlinelibrary.wiley.com/doi/abs/10.1002/minf.201200133)\
Download [link](https://sourceforge.net/projects/ambit/files/Ambit2/AMBIT%20applications/tautomers/)
