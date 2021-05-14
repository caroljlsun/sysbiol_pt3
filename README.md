# Investigation of DNA structural motifs, CpG methylation and their role in transcriptional regulation

Welcome to my github repo for my part III project.

The code for this systems biology project is mostly within the folder "Code". Models that were trained, mostly the CNNs, are within the folder "Models"

This project mainly investigates how we can predict the presence of G-quadruplexes within promoters, in the mouse genome. 3 methods of machine learning (ML) were assessed for this purpose: Support Vector Machines (SVMs), Random Forests (RFs), Convolutional Neural Network (CNNs). 
The best predictive model is a convolutional neural network, named _ML_5_model_ , with an AUC of 0.88. It's model architecture is:

![image](https://user-images.githubusercontent.com/61421828/118273050-5b628780-b4bb-11eb-8402-b5cd523a43a9.png) 

The features used in training the models are shown:
![image](https://user-images.githubusercontent.com/61421828/118273087-674e4980-b4bb-11eb-8ccd-bd273654be1b.png)

The second part of the project investigates if there are cell type specfic features associated with G4s. This was investigated by listing the top 100 predicted G4s from SVMs, RFs, & CNNs. Then the gene names were identified from these G4 locations, and investigated using [enRichr](https://maayanlab.cloud/Enrichr/) as well as "goProfiles" in R. Generally, all ML methods were found to score promoters highly in white blood cell related genes. This is quite interesting, considering how many features used in training have B & T cell information. For example, the cell type associated with the genes the SVMs found: 

![image](https://user-images.githubusercontent.com/61421828/118273327-b4322000-b4bb-11eb-9d2f-e6e12d4d7776.png)

Finally, there is a positive correlation between high G4 score and gene expression:

![image](https://user-images.githubusercontent.com/61421828/118273361-bbf1c480-b4bb-11eb-888e-f0968be09ca1.png)


## Guide to code uploaded
Here is a quick description of the code uploaded in this github repo:

_Collapse_arrays_ changes the feature vector into a form that is useable for SVMs and RFs

_Correlations_ and _Correlations_negative_G4s_ explores the correlations of the features used in ML training

_DeSeq2_ and _Expression_ files explore the expression data in the feature vector (FV)

_FV__ * _G4__ * and _G_quadruplex_ _ * _ *  files process raw G4 counts into a form useable in the FV

_FV_ _* _minimal_ _ _1_ files make differently sized FVs, used to assess ML code before using the full size FV

_Feature_masking_ trains a CNN using the one-hot encoding of the promoters only

_Feature_masking_ _2_ trains a CNN using all feature _except_ the one-hot encoding

_G4_scores_CNN_1_ finds the top 100 G4 predictions using the best CNN model

_GC_content_ calculates the GC content of a promoter sequence using a rolling window of 25 bps

_GO_profiles_ compares the gene lists predicted by SVMs, RFs and CNNS

_Histones_ * _ * _ processes raw histone data

_Kmer_PCA_UMAP_copy_ and 
_write_2_fasta_copy_ conducts k-mer analysis of the promoters with and without G4s

_Methylated_positions_ processes raw methylation data to get the positions of methylated bases in promoters

_Methylkit_2_ uses the methylation data to get DNAshaper features for the FV

_R_loops_ obtains R loop positions

_Random_forests_1_ trains multiple RFs

_SVM_1_ trains multiple SVMs

_Saliency_maps__*_ finds the saliency maps of the multiple CNNs trained in this project

_Venn diagram of G quadruplexes_ assesses briefly how many G4s are shared between different experimental conditions.

_compile_sequences_ lists all the sequences of the promoters used
## Guide to models uploaded
Here is a quick description of the models uploaded in this github repo:

_ML_1_ a CNN which uses 2D convolution, AUC = 0.7772

_ML_4_model_ a CNN which uses the full FV and 1D convolution, AUC = 0.84

_ML_5_model_ a CNN which is feature masked, only looks at one-hot encoding, AUC = 0.88

_ML_6_model_ a CNN which is feature masked, looks at all features except one-hot encoding, AUC = 0.81
