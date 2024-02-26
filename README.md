# MoE
A mixture of experts is one of semisupervised learning models. 
It consists of a gating network and expert networks.  

# Finite mixture of regressions.R
This file is just based on regmixEM function in mixtools package in R. 

(Benaglia, T., Chauveau, D., Hunter, D. R., & Young, D. S. (2010). mixtools: an R package for analyzing mixture models. Journal of statistical software, 32, 1-29.)

# Mixture of experts.R 
Basic codes for mixture of experts. 
There are 6 functions in this file. 
1. softmax: Calculation of softmax function.
2. IRLS: Iteratively Reweighted Least Squares estimation algorithm for the gating network.
3. MoE: The main function of the mixture of experts. EM algorithm is used for estimating parameters. 
4. MoE_iter: Due to the inherent multimodality of the mixture likelihood, it is advisable to consider
multiple initial values when applying the EM algorithm in practice. For fixed numbers of experts, this function selects the one that yields the highest likelihood.
5. MoE_best: This function selects the best numbers of experts. 
6. MoE_prediction: This function predicts y based on given x. 

