## Predicting the target-landscape of kinase-inhibitors using 3D convolutional neural network on densely covered chemogenomic data
Georgi K. Kanev, Yaran Zhang, Albert J. Kooistra, Andreas Bender, Rob Leurs, David Bailey, Thomas Würdinger, Chris de Graaf, Iwan J.P. de Esch, Bart A. Westerman

#### [Data and results](/data)
#### [Scripts](/scripts)

### Abstract
Many therapies in clinical trials are still based on single drug-single target relationships. To further extend this concept to multi-target approaches using multi-targeted drugs, we developed a machine learning pipeline to unravel the target landscape of kinase inhibitors. This pipeline, which we call 3D-KINEssence, uses a new type of protein fingerprints (3D FP) based on the structure of kinases generated through a 3D convolutional neural network (3D-CNN). These 3D-CNN kinase fingerprints were matched to molecular Morgan fingerprints to predict the targets of each respective kinase inhibitor based on available bioactivity data. The performance of the pipeline was evaluated on two test sets: a sparse drug-target set where each drug is matched in most cases to a single target, but also on a densely-covered drug-target set where each drug is matched to most if not all targets. This latter set is more difficult to train given its non-exclusive character. Our model's root-mean-square error (RMSE) based on the two datasets was 0.61 and 0.775, respectively. These results indicate that 3D FP can confidently predict the target landscape of kinase inhibitors at around 0.7 log unit of bioactivity, also when using stringent criteria not frequently used in the field (i.e., the dense dataset). Our strategy can be utilized in proteochemometric or chemogenomic workflows by consolidating the target landscape of kinase inhibitors and their patient-specific targets.
