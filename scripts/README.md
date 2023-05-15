### Requires:
pandas\
numpy\
rdkit\
sklearn\
seaborn\
matplotlib

### Workflow:
#### Prepare data:
(1) python prep_chembl.py\
(2) python prep_vanwesten.py\
(3) python merge_db_sets.py\
(4) bash gen_tautomers.sh\
(5) python extr_tautomers.py\
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
