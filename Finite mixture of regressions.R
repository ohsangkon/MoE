############################################
# Finite mixture of regressions #######################

library(mixtools)
library(boot)
library(MASS)

############################################
init_value <- function(k = 2, x, y, samp = 0.95){
	### inital values ####
	library(stats)
	library(MASS)

	model <- NULL
	x = as.matrix(x)
	y = as.matrix(y)

	prop <- NULL
	beta <- NULL
	sigma <- NULL

	clust_x <- list()
	clust_y <- list()
	init <- list()

	ind <- sample(nrow(x), round(nrow(x)*samp))
	x_sample <- as.matrix(x[ ind , ])
	y_sample <- as.matrix(y[ ind,  ])

	clust <- kmeans(as.matrix(cbind(x_sample, y_sample)),k)  # split data

	for(j in 1 : k){  # estimate
		clust_x[[j]] <- as.matrix(x_sample[which(clust$cluster == j),])
		clust_y[[j]] <- as.matrix(y_sample[which(clust$cluster == j),])


		design_X <- model.matrix(~., data=data.frame(clust_x[[j]]))
		beta = cbind(beta, ginv(t(design_X) %*% design_X) %*% t(design_X) %*% clust_y[[j]])

		sigma = append(sigma, sum((clust_y[[j]] - design_X %*% beta[,j])^2)/nrow(clust_x[[j]]))

	}

	model$lambda <- clust$size / sum(clust$size)
	model$beta <- beta
	model$s <- sigma
	model$k = k

	return(model)
}


regmixEM2 = function (y, x, lambda = NULL, beta = NULL, sigma = NULL, k = 2, 
    addintercept = TRUE, arbmean = TRUE, arbvar = TRUE, epsilon = 1e-08, 
    maxit = 10000, verb = FALSE) 
{
    if (arbmean == FALSE && arbvar == FALSE) {
        stop(paste("Must change constraints on beta and/or sigma!", 
            "\n"))
    }
    s = sigma
    if (addintercept) {
        x = cbind(1, x)
    }
    n <- length(y)
    p <- ncol(x)
    tmp = init_value(k = k, x[,-1], y, samp = 0.95)
	
    lambda <- tmp$lambda
    beta <- tmp$beta
    s <- tmp$s
    k <- tmp$k
    diff <- 100000
    iter <- 0
    xbeta <- x %*% beta
    res <- (y - xbeta)^2
    if (arbmean == FALSE) {
        res <- sapply(1:k, function(i) res)
    }
    comp <- t((lambda/sqrt(2 * pi * s^2)) * t(exp(-t(t(res)/(2 * 
        s^2)))))
    obsloglik <- sum(log(apply(comp, 1, sum)))
    ll <- obsloglik
    z = matrix(nrow = n, ncol = k)
    restarts <- 0
    while (diff > epsilon && iter < maxit) {
	z = comp/ apply(comp,1,sum)

        lambda.new <- apply(z, 2, mean)
        if (sum(lambda.new < 1e-08) > 0 || is.na(sum(lambda.new))) {
            sing <- 1
        }
        else {
            if (arbmean == FALSE) {
                if (addintercept) {
                  beta.new <- lm(y ~ x[, -1], weights = apply(t(t(z)/(s^2)), 
                    1, sum))$coef
                }
                else beta.new <- lm(y ~ x - 1, weights = apply(t(t(z)/(s^2)), 
                  1, sum))$coef
            }
            else {
                if (addintercept) {
                  lm.out <- lapply(1:k, function(i) lm(y ~ x[, 
                    -1], weights = z[, i]))
                }
                else lm.out <- lapply(1:k, function(i) lm(y ~ 
                  x - 1, weights = z[, i]))
                beta.new <- sapply(lm.out, coef)
            }
            xbeta.new <- x %*% beta.new
            res <- (y - xbeta.new)^2
            if (arbmean == FALSE) {
                res <- sapply(1:k, function(i) res)
            }
            if (arbvar) {
                s.new <- sqrt(sapply(1:k, function(i) sum(z[, 
                  i] * (res[, i]))/sum(z[, i])))
            }
            else s.new <- sqrt(sum(z * res)/n)
            lambda <- lambda.new
            beta <- beta.new
            xbeta <- x %*% beta
            s <- s.new
            sing <- sum(s < 1e-08)
            comp <- lapply(1:k, function(i) lambda[i] * dnorm(y, 
                xbeta[, i * arbmean + (1 - arbmean)], s[i * arbvar + 
                  (1 - arbvar)]))
            comp <- sapply(comp, cbind)
            compsum <- apply(comp, 1, sum)
            newobsloglik <- sum(log(compsum))
        }
        
            diff <- newobsloglik - obsloglik
				if(diff < 0 | is.nan(diff) == TRUE | is.na(diff) == TRUE){
					cat("Warnings! We recommend fewer components than the number you specified.", "\n")
					diff = Inf
					obsloglik = -Inf
					break		
				}
            obsloglik <- newobsloglik
            ll <- c(ll, obsloglik)
            iter <- iter + 1
            if (verb) {
                cat("iteration=", iter, "diff=", diff, "log-likelihood", 
                  obsloglik, "\n")
            
        }
    }
    scale.order = order(s)
    sigma.min = min(s)
    if (iter == maxit | obsloglik == Inf) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
    if (arbmean == FALSE) {
        z = z[, scale.order]
        names(beta) <- c(paste("beta", ".", 0:(p - 1), sep = ""))
        colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
        a = list(x = x, y = y, lambda = lambda[scale.order], 
            beta = beta, sigma = sigma.min, scale = s[scale.order]/sigma.min, 
            loglik = obsloglik, posterior = z[, scale.order], 
            all.loglik = ll, restarts = restarts, ft = "regmixEM")
        class(a) = "mixEM"
        a
    }
    else {
        rownames(beta) <- c(paste("beta", ".", 0:(p - 1), sep = ""))
        colnames(beta) <- c(paste("comp", ".", 1:k, sep = ""))
        colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
        a = list(posterior = z, prop = lambda, beta = beta, 
            sigma = s, loglik = obsloglik, all.loglik = ll, 
            restarts = restarts, ft = "regmixEM")
        class(a) = "mixEM"
        a
    }
}

			   
### Example ###
data(NOdata)
attach(NOdata)
y = Equivalence
x = NO
#x = cbind(NO, rnorm(nrow(NOdata), mean = 0, sd = 1))

regmixEM2(y, x, k = 2)



