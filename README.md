# MoE
A mixture of experts is one of semisupervised learning models. 
It consists of a gating network and expert networks.  

# Finite mixture of regressions.R
This file is just based on regmixEM function in mixtools package in R. 
This function is utilized for inital values of parameters in the mixture of experts. 
You should save this function if you would like to use functions in Mixture of experts.R.

(Benaglia, T., Chauveau, D., Hunter, D. R., & Young, D. S. (2010). mixtools: an R package for analyzing mixture models. Journal of statistical software, 32, 1-29.)

# Mixture of experts.R 
Basic codes for mixture of experts. 
There are 6 functions in this file. 
1. softmax: Calculation of softmax function.
2. IRLS: Iteratively Reweighted Least Squares estimation algorithm for the gating network (only 1 iteration).
3. MoE: The main function of the mixture of experts. EM algorithm is used for estimating parameters. 
4. MoE_iter: Due to the inherent multimodality of the mixture likelihood, it is advisable to consider
multiple initial values when applying the EM algorithm in practice. For fixed numbers of experts, this function selects the one that yields the highest likelihood.
5. MoE_best: This function selects the best numbers of experts. 
6. MoE_prediction: This function predicts y based on given x. 

Arguments
1. k: the number of experts
2. x: covariates
3. y: response variable
4. z: latent variable indicating membership
5. alpha: parameters in gating networks
6. exit: If the calculation of derivative is impossible, IRLS algorithm is stopped and exit is updated to exit + 1.
7. eps: tolerance
8. maxit: maximum iteration of EM algorithm
9. iter: maximum iteration of MoE function with fixed numbers of experts.
10. range: range of the number of experts
11. model: results of MoE
