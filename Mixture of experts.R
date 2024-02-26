############################################
# Mixture of linear experts #######################

library(mixtools)
library(boot)
library(MASS)

############################################

softmax = function(k = 2, x, alpha){
	## exp(x)/sum_k=1^K exp(x) ##

	X = cbind(1,x)
	tmp = exp(X %*% alpha)
	denominator = apply(tmp, 1, sum)

	result = tmp
	for(j in 1 : k){
		result[,j] = tmp[,j]/denominator
	}
	return(result)
}


IRLS = function(k = 2, x, alpha = NULL, z = NULL, exit = NULL){
	## Iteratively reweighted least squares estimation (IRLS) ###
    ## alpha_new = alpha - (X^T W_j X)^-1 X^T (Z_j - B_j) ###
	X = cbind(1,x)
	p = ncol(as.matrix(x))
	n = nrow(as.matrix(x))

	B = softmax(k, x, alpha)

	first_derivative = alpha
	for(j in 1 : (k-1) ){
		first_derivative[,j] = t(X) %*% (z[,j] - B[,j])
	}

	w = matrix(0, nrow = n, ncol = k)
	tmp = exp(X %*% alpha)
	denominator = apply(tmp, 1, sum)^2

	alpha_new = alpha 
	for(j in 1: (k-1) ){
		if(k > 2){
			w[,j]= - tmp[,j] * apply(tmp[,-j], 1, sum)/denominator
		}else{
			w[,j]= - tmp[,j] /denominator
		}
		if(sum(is.nan(w[,j])) > 0){
			cat("Warnings! We guess there are too many experts", "\n")
			exit <<- exit + 1 
			return(alpha_new)
			break
		}
		if(min(svd(t(X) %*% diag(w[,j]) %*% X)$d) > 0.1){
			alpha_new[,j] = alpha[,j] - solve(t(X) %*% diag(w[,j]) %*% X) %*% first_derivative[,j]
		}else{
			alpha_new[,j] = alpha[,j] - ginv(t(X) %*% diag(w[,j]) %*% X) %*% first_derivative[,j]
		}
	}

	return(alpha_new)
}


MoE = function(k = 2, y, x, eps = 10^-5, maxit = 500){
	### main function ###

	X = cbind(1,x)
	p = ncol(as.matrix(x))
	n = nrow(as.matrix(x))
	alpha = matrix(0, ncol = k , nrow = p+1 )
	iteration = 0
	exit = 0

	### initial ###
	if(k != 1){
		### finite mixture of regressions ###
		int = regmixEM2(y = y, x = x, k = k, epsilon = eps)
  		beta = int$beta	# expert network
 	   sigma = int$sigma # expert network

		z = int$posterior # latent variable

		for(i in 1 : 10){
			alpha = IRLS(k, x, alpha, z, exit)
		}

   prop = softmax(k, x, alpha)
   xbeta <- X %*% beta
   res <- (y - xbeta)^2

   comp <- t(t(prop) / sqrt(2*pi*sigma^2) * t(exp(-t(t(res)/(2 * sigma^2)))))
	obslik <- sum(log(apply(comp, 1, sum)))

	diff = 10000
	
	lik_record = obslik

	### EM algorithm ###
	while(diff > eps && iteration < maxit){
		exit = 0

		### E-step ###
		z = comp/ apply(comp,1,sum)

		### M-step ###

 	   lm_out = lapply(1:k, function(i) lm(y ~ x, weights = z[, i]))
		
		beta_new = sapply(lm_out, coef)
		xbeta_new = X %*% beta_new
	    res_new = (y - xbeta_new)^2
        sigma_new = sqrt(sapply(1:k, function(i) sum(z[, i] * (res_new[,i]))/sum(z[, i])))
		for(i in 1 : 10){
			alpha_new = IRLS(k, x, alpha, z, exit)
			alpha = alpha_new
		}

		if(exit == 10){
			cat("Warnings! We recommend fewer experts than the number you specified.", "\n")
			diff = Inf
			lik_record = -Inf
			break			
		}

		### diff = newlik - obslik ###
		prop_new = softmax(k, x, alpha_new)
	    comp_new <- t(t(prop_new) / sqrt(2*pi*sigma_new^2) * t(exp(-t(t(res_new)/(2 * sigma_new^2)))))
 	    newlik = sum(log(apply(comp_new, 1, sum)))
	 
		diff  = newlik - obslik
        obslik = newlik
		lik_record = append(lik_record, obslik)

		iteration = iteration + 1
		cat("iteration=", iteration, "diff=", diff, "log-likelihood", 
                  obslik, "\n")

		if(diff < 0 | is.nan(diff) == TRUE | is.na(diff) == TRUE){
			cat("Warnings! We recommend fewer experts than the number you specified.", "\n")
			diff = Inf
			lik_record = -Inf
			break		
		}

		comp = comp_new
		beta = beta_new
		sigma = sigma_new
		alpha = alpha_new


	}
	}else{
		### just linear regression ###
		iteration = iteration + 1
		beta = solve(t(X) %*% X) %*% t(X) %*% y
		xbeta = X %*% beta
		res = X %*% beta
		sigma = t(res) %*% (res) / n
	   comp <- dnorm(res, mean = xbeta, sd = sigma)
  		lik_record <- sum(log(comp))
 		z = matrix(1, nrow = n, ncol = 1)
	}


	## Output ##
	cat("number of iterations=", iteration, "\n")

	BIC = -2 * lik_record[length(lik_record)] + ((p+1) * (2*k-1) + k ) * log(n)
	ICL = BIC - sum(apply(z * log(z+10^-10), 2, sum))

	cluster = apply(z, 1, which.max)


   results = list(posterior = z, beta = beta, sigma = sigma, alpha = alpha, 
        				BIC = BIC, ICL = ICL, cluster = cluster, all.loglik = lik_record)	

	#plot(lik_record, type = "l")

	return(results)	

}

MoE_iter = function(iter = 10, k = 2, y, x, eps = 10^-5){
	## find the best MoPLE with the highest likelihood ##
	candidate = list() # MoE
	loglik = NULL
	for(j in 1 : iter){
		candidate[[j]] = MoE(k = k, y, x, eps = eps)
		loglik = append(loglik, candidate[[j]]$all.loglik[length(candidate[[j]]$all.loglik)])

		cat("number of iterations=", j, "\n")
	}
	return(candidate[[which.max(loglik)]])
}


MoE_best = function(iter = 10, range = 1:5, y, x, eps = 10^-5){
	## Find the best number of experts ##
	candidate = list()
	BIC = NULL
	ICL = NULL
	for(j in range){
		cat("number of compoents=", j, "\n")
		candidate[[j]] = MoE_iter(iter = iter, k = j, y, x, eps)
		BIC = append(BIC, candidate[[j]]$BIC)
		ICL = append(ICL, candidate[[j]]$ICL)
	}
	best_BIC = which.min(BIC)
	best_ICL = which.min(ICL)
	
	result = list()
	result[[1]] = list("Best based on BIC" = best_BIC, "Best based on ICL" = best_ICL)
	result[[2]] = candidate[[best_BIC]]
	result[[3]] = candidate[[best_ICL]]
	names(result) = c("# of experts", "Based on BIC", "Based on ICL")

	return(result)
}

MoE_prediction = function(model, x){
	### Prediction using MoE ###
	X = cbind(1,x)
	beta = model$beta 
	alpha = model$alpha

	experts = X %*% beta
	gating = softmax(k = ncol(beta), x, alpha)

	y_hat = apply(gating * experts, 1, sum)
	return(y_hat)
}

###########################################################
############# Examples ####################################
###########################################################
### simulation ### 
# Generation
set.seed(530)
eps1 = rnorm(200, 0, 0.5)
eps2 = rnorm(200, 0, 0.3)

x1 = rnorm(200, mean = 0, sd = 2)  ; x2 = rnorm(200, mean = 10, sd = 1)
y1 = cbind(1, x1) %*% c(-1,1) + eps1; y2 = cbind(1, x2) %*% c(5,-1)+ eps2
x = c(x1, x2); y = c(y1, y2)

plot(y ~ x)

# fitting
set.seed(530)
model = MoE_best(iter = 10, range = 1:3, y, x)
model[[1]]
model[[2]]
model[[3]]

# prediction
y_hat = MoE_prediction(model[[3]], x)

# result
plot(y ~ x, col = model[[3]]$cluster, pch = model[[3]]$cluster, cex = 1.5)
prediction = cbind(x,y_hat)
prediction = prediction[order(x),]

mat = cbind(x, cbind(1,x) %*% model[[3]]$beta)
lines(mat[order(x),c(1,2)], col = "black", lwd = 3 )
lines(mat[order(x),c(1,3)], col = "red", lwd = 3)

lines(prediction, col = "blue", lwd = 1, type = "p")

### Real dataset ###
data(NOdata)
attach(NOdata)
y = Equivalence
x = NO

# fitting
set.seed(530)
model = MoE_best(iter = 10, range = 1:3, y, x)
model[[1]]

# prediction
y_hat = MoE_prediction(model[[3]], x)

# result
plot(y ~ x, col = model[[3]]$cluster, pch = model[[3]]$cluster, cex = 1.5)
prediction = cbind(x,y_hat)
prediction = prediction[order(x),]

mat = cbind(x, cbind(1,x) %*% model[[3]]$beta)
lines(mat[order(x),c(1,2)], col = "black", lwd = 2)
lines(mat[order(x),c(1,3)], col = "red", lwd = 2)

lines(prediction, col = "blue", lwd = 1, type = "p")

#####################################################