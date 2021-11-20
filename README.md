## Alterations in Connectome Dynamics in ASD
This repository provides data and codes used in our paper 'Alterations in Connectome Dynamics in Autism Spectrum Disorder: A Harmonized Mega- and Meta-Analysis Study Using the ABIDE Dataset'.
## Data
1. Participant information, including Participant ID, site, scaning time, TR, eye status, group, age, IQ, Mean FD.
  <br> Participant.mat
2. Random 512 parcellation
  <br> GroupMask_Random512.nii
3. Global and regional module variability before and after Combat harmonization
  <br> MV_Harmonized_Mean.mat
  <br> MV_Harmonized_std.mat
  <br> MV_Harmonized_Zscore.mat
  <br> Mean_MV.mat
  <br> STD_MV.mat
  <br> MV_Zscore.mat
4. Participants with SRS score
  <br> Participant_SRS_RAW_TOTAL.csv
5. Brain mask used for gene association analysis
  <br> LeftHemi222.nii
6. Gene Ontology enrichment analysis
  <br> OrderedGeneList.txt
7. Gene set enrichment analysis
  <br> PLS1_OrdGeneList.csv
  <br> CustomGeneSets.gmt
  ## Codes
1. Multilayer community detection [1]
<br> http://netwiki.amath.unc.edu/GenLouvain 
2. Modular variability [2]
<br> scaled_inclusivity.m
<br> scaled_inclusivity_wei.m
<br> ami.m
3. Case-control comparison analysis
<br> RegionalModel_Mega.r
<br> RegionalModel_Meta.r 
4. Correlation between ASD-related alterations and cognitive terms [3]
<br> https://neurosynth.org/decode/ 
5. Prediction of Social Impairments Using Connectome Dynamics
<br> SVR_LOOCV_Prediction.m
6. Gene expression data preprocessing [4]
<br> https://github.com/BMHLab/AHBAprocessing 
7. Association between alterations in connectome dynamics and gene expression profiles
<br> GeneAssociationAnalysis.m
8. Gene Ontology enrichment analysis [5]
<br> http://cbl-gorilla.cs.technion.ac.il/ 
9. Gene set enrichment analysis
<br> ClusterProfiler_GSEA.r
## Reference
<br>[1] Jeub LGS, Bazzi M, Jutla IS, Mucha PJ. A generalized Louvain method for community detection implemented in MATLAB. Version 2.1 [software]. 2012 [updated 2016 Nov; cited 2021 Oct 3]. Available from: http://netwiki.amath.unc.edu/GenLouvain.
<br>[2] Liao X, Cao M, Xia M, He Y. Individual differences and time-varying features of modular brain architecture. Neuroimage. 2017;152:94-107.
<br>[3] Yarkoni T, Poldrack RA, Nichols TE, Van Essen DC, Wager TD (2011): Large-scale automated synthesis of human functional neuroimaging data. Nat Methods. 8:665-670.
<br>[4] Arnatkeviciute A, Fulcher BD, Fornito A (2019): A practical guide to linking brain-wide gene expression and neuroimaging data. Neuroimage. 189:353-367.
<br>[5] Eden E, Navon R, Steinfeld I, Lipson D, Yakhini Z (2009): Gorilla: A tool for discovery and visualization of enriched go terms in ranked gene lists. BMC Bioinformatics. 10:48.

   
