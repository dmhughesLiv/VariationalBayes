# VariationalBayes
THis repository contains codes to accompany an article in Biostatistics "Fast approximate inference for multivariate longitudinal data".
The codes will allow reproduction of Figure 5 in the Supplementary material and is based on the Primary Biliary Cirhosis data available within the mixAK package in R.

We make use of three R scripts provide by Lee and Wand in support of their paper (Lee and Wand, 2016, Biometrical Journal: Streamlined mean field variational Bayes for longitudinal and multilevel data analysis), specifically "trapint.r", "accVarApp.r" and "summMCMC.r". These codes can be accessed at:https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201500007

The function mglmm_New.R fits the multivariate generalised linear mixed models using the mean field variational Bayes approach described in the paper. 

The script PBCJointResp10Markers_mixAK.R fits the 10 marker multivariate mixed model described in the paper and produces Supplementary Figure 5, along with additional graphics not provided in the paper.

