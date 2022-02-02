# used to generate beta surface that is used for comparison with functional regression
make_surface = function(basis_time, alpha, beta_s = 1){
	###
	exp_f = function(s, t, alpha){
		exp(-alpha * (t - s)) * (s <= t)
	}
	
	D = length(basis_time)
	s = t = basis_time
	
	exp_surface = t(outer(s, t, exp_f, alpha = alpha) * beta_s)
	colnames(exp_surface) = paste0("s", 1:D); rownames(exp_surface) = paste0("t", 1:D)
	
	attr(exp_surface, "t") = t; attr(exp_surface, "s") = s
	return(exp_surface)
}

## function to get joint log likelihood of model for monitoring convergence
logLik = function(data, N, D, Kt, Dstar, delta_coefs, beta_coefs, xstar_ls, sigma, lambda_d, Pen,
									random_effects){
	
	if(random_effects){
		lik_y_sigma = sigma * diag(D)  + Dstar %*% (lambda_d * solve(Pen)) %*% t(Dstar)
	}else{
		lik_y_sigma = sigma * diag(D)  
	}
	
	lik_y = rep(NA, length.out = N)
	# loop over trials
	trial_ids = unique(data$trial)
	for(i in 1:length(trial_ids)){
		trial_df = filter(data, trial == trial_ids[i])
		mu = (trial_df$y0_star + as.numeric(xstar_ls[[i]] %*% as.vector(beta_coefs)))
		
		lik_y[i] = dmvnorm(trial_df$value, mean = mu, sigma = lik_y_sigma,
											 checkSymmetry = FALSE,
											 log = TRUE) 
	}
	
	list(ylik = sum(lik_y))

}


# this function puts together Xstar for ith trial
# thetaStar is stored value of integrated intercept term at a particular value of alpha
make_Xstar = function(xdata, alpha, spline_basis, basis_time, thetaStar){
	forcing_functions = as.matrix(xdata)
	D = nrow(forcing_functions)
	P = ncol(forcing_functions) # includes intercept term
	Kt = ncol(spline_basis)
	
	xstar_i = matrix(NA, nrow = D, ncol = P * Kt)
	for(p in 1:P){
		p_index = (p-1) * Kt + 1
		
		if(p > 1){
			xstar_ip = integrate_alpha(basis_time, spline_basis, alpha, covar = forcing_functions[,p])
		}else{
			xstar_ip = thetaStar
		}
		xstar_i[, p_index:(p * Kt)] = xstar_ip
	} # end for loop
	
	return(xstar_i)
}

# This function integrates the exponential surface containing alpha
# used to create Dstar and Xstar terms
# returns Gstar, a D x Kt matrix that will be multiplied by beta_coefs or delta_coefs in later steps
# covar is either a length D vector of ones or covariate values for a particular covariate and subjects
integrate_alpha = function(basis_time, spline_basis, alpha, covar = NULL){
	D = length(basis_time)
	if(is.null(covar)){
		covar = rep(1, length.out = D)
	}
	
	exp_f = function(s, t = 1, basis_time, alpha, covar){
		exp(-alpha * (t - s))  *  (s < t) *  covar[which(basis_time == s)] * spline_basis[which(basis_time == s),]
	}
	
	Gstar = sapply(basis_time, FUN = function(t){cumtrapz(basis_time, exp_f(basis_time, t, basis_time, alpha, covar))[D,]})
	return(t(Gstar))
}


# alpha loss for no random effects
alpha_loss = function(data, x_data, beta_coefs, alpha, spline_basis){
	
	data = data %>% mutate(y0_star = Yo * exp(- alpha * time))
	
	basis_time = filter(data, trial == 1)$time
	thetaStar = integrate_alpha(basis_time, spline_basis, alpha, covar = NULL)
	xstar_ls = lapply(x_data, make_Xstar, alpha = alpha, spline_basis = spline_basis, basis_time = basis_time, thetaStar = thetaStar)
	xstar = do.call(rbind, xstar_ls)
	
	#Yo * exp(- alpha0 * time)
	#unique(subject_df$Yo) * exp(-alpha * subject_df$time)
	
	sum((data$value - data$y0_star - as.numeric(xstar %*% as.vector(beta_coefs)))^2)
}


alpha_loss_re = function(data, x_data, spline_basis, basis_time, beta_coefs, alpha, 
												 random_effects, delta_coefs, C){
	
	if(!random_effects){
		alpha_loss(data, x_data, beta_coefs, alpha, spline_basis)
	}else{
		Kt = dim(spline_basis)[2]
		Dstar = integrate_alpha(basis_time, spline_basis, alpha, covar = NULL) # same as thetaStar
		deltaStar = Dstar %*% delta_coefs
		Dstar_squared = crossprod(Dstar)
		trace_Dstar_C = sum(diag(Dstar_squared %*% C))
		trial_ids = unique(data$trial)
		xstar_ls = lapply(x_data, make_Xstar, alpha = alpha, 
											spline_basis = spline_basis, basis_time = basis_time, 
											thetaStar = Dstar)
		#xstar = do.call(rbind, xstar_ls)
		
		#sum((data$value - data$y0_star - as.numeric(xstar %*% as.vector(beta_coefs)) - deltaStar)^2)
		
		rss = rep(NA, length(trial_ids))
		for(i in 1:length(trial_ids)){
			subject_df = filter(data, trial == trial_ids[i])
			ystar_i = unique(subject_df$Yo) * exp(-alpha * subject_df$time)
			deltastar_i = deltaStar[, i]
			# maybe below shouldn't be subject specific or maybe tehre is a plus when there should be a minus
			deltaStar_squared_i = trace_Dstar_C + crossprod(delta_coefs[,i], Dstar_squared %*% delta_coefs[,i])
			Y_centered = subject_df$value - ystar_i - as.numeric(xstar_ls[[i]]  %*% as.vector(beta_coefs))


			rss[i] = crossprod(Y_centered) -
				2 * crossprod(Y_centered, deltastar_i) +
				deltaStar_squared_i

		}
		return(sum(rss))
	}
	
	
}



