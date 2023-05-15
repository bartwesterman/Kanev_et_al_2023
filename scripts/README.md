### Required packages for the python scripts:
GE  : pandas, numpy\
CH  : rdkit\
ML  : sklearn\
VIS : seaborn, matplotlib

### Workflow:
##### Prepare data:
* Prepare the tautomer files ../data/chembl/chembl_tauromers.csv & ../data/vanwesten/vanwesten_tauromer.csv\
(1) bash gen_tautomers.sh \
(2) python extr_tautomers.py

* Unzip the files in ../data/\
(0) unzip ../data/*.zip

* Prepare the data\
(3) python prep_chembl.py\
(4) python prep_vanwesten.py\
(5) python merge_db_sets.py\
(6) python rm_cpds.py\
(7) python sel_kinases.py\
(8) python prep_dense_set.py\
(9) python prep_sparse_set.py

#### Prepare 3D-CNN model:

#### Analysis:

### Software:
#### Ambit-Tautomer: An Open Source Tool for Tautomer Generation
Nikolay T Kochev, Vesselina H Paskaleva, Nina Jeliazkova\
Mol Inform. 2013 Jun;32(5-6):481-504.\
PMID: [27481667](https://pubmed.ncbi.nlm.nih.gov/27481667/)\
DOI: [10.1002/minf.201200133](https://onlinelibrary.wiley.com/doi/abs/10.1002/minf.201200133)\
Download link: https://sourceforge.net/projects/ambit/files/Ambit2/AMBIT%20applications/tautomers/
