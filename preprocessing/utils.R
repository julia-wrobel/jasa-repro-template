
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

# This function makes DStar and xstar matrices by doing integration
# Integration is alway done for a single subject at a time
make_thetaStar = function(surface, spline_basis, xp = NULL){
	
	s = attributes(surface)$s
	#ds = c(diff(s)[1], diff(s))
	ds = c(0, diff(s))
	
	D = nrow(spline_basis);  Kt = ncol(spline_basis)
	one_vector = matrix(1, nrow = Kt, ncol =  1)
	
	if( is.null(xp) ) {xp = rep(1, D)}
	
	thetaStar = matrix(NA, D, Kt)
	for(t_index in 1:D){
		thetaStar[t_index,] = colSums(kronecker(surface[t_index,] * xp, t(one_vector)) * spline_basis *ds)
	}
	return(thetaStar)
}


# this function puts together Xstar for ith trial
make_Xstar = function(xdata, alpha, spline_basis, time){
	forcing_functions = as.matrix(xdata)
	D = nrow(forcing_functions)
	P = ncol(forcing_functions) # includes intercept term
	Kt = ncol(spline_basis)
	
	surface = make_surface(time, alpha)
	xstar_i = matrix(NA, nrow = D, ncol = P * Kt)
	for(p in 1:P){
		p_index = (p-1) * Kt + 1
		xstar_ip = make_thetaStar(surface, spline_basis, xp = forcing_functions[,p])
		
		xstar_i[, p_index:(p * Kt)] = xstar_ip
	} # end for loop
	
	return(xstar_i)
}

make_y0 = function(data, xstar, alpha, beta_coefs,
									 slice_num = 1, mle = FALSE, use_re = FALSE){
	#yi0 = sum(data$value - deltaStar - xstar %*% beta_coefs)/sum(exp(-alpha * data$time))
	#yi0 = t(exp(alpha * data$time)) %*% (data$value - deltaStar - xstar %*% beta_coefs) %>% as.numeric()
	
	#yi0 = (unique(data$Yo) + (xstar %*% beta_coefs)[1]) %>% as.numeric()
	#yi0 = exp(alpha * data$time[1]) * (data$value- deltaStar - xstar %*% beta_coefs)[1]
	#yi0 = yi0/length(data$time) %>% as.numeric()
	#yi0 = as.numeric(yi0)
	
	if(mle){
		if(use_re){
			yi0 = sum(exp(-alpha * data$time) * (data$value - data$deltaStar - xstar %*% beta_coefs)) / sum(exp(- 2 * alpha * data$time))
		}else{
			yi0 = sum(exp(-alpha * data$time) * (data$value  - xstar %*% beta_coefs)) / sum(exp(- 2 * alpha * data$time))
		}
	}else{
		# take mean of first 3 fitted values
		yi0 = data %>% slice(1:slice_num) %>% summarize(yi0 = mean(yhat)) %>% pull(yi0)
	}
	
	data$Yo = yi0
	data
}


## this is the model function for nls()
## here data is the full dataset
make_alpha = function(data, x_data, beta, alpha, lambda_d = NULL, sigma = NULL, spline_basis = NULL, 
											random_intercept = FALSE){
	# to make this work for random effects could I just add half of the weird extra part to Y?
	basis_time = filter(data, trial == 1)$time
	D = length(basis_time)
	exp_surface = make_surface(basis_time, alpha)
	s = attributes(exp_surface)$s
	#ds = c(diff(s)[1], diff(s))
	ds = c(0, diff(s))
	
	trial_ids = unique(data$trial)
	Y = matrix(NA, nrow = length(trial_ids), ncol = D)
	
	for(i in 1:length(trial_ids)){
		subject_df = filter(data, trial == trial_ids[i])
		

		Y[i,] = unique(subject_df$Yo) * exp(-alpha * subject_df$time) + 
			colSums(t(exp_surface) * rowSums(x_data[[i]] * beta) * ds)
		
	} # end for loop
	
	y = as.vector(t(Y))
	return(y)
} # end make_alpha

# this will need to be edited when we add in random effects
alpha_loss = function(data, x_data, beta, alpha, lambda_d = NULL, sigma = NULL, spline_basis = NULL, 
											random_intercept = NULL){
	sum((data$value - make_alpha(data, x_data, beta, alpha))^2)
	#sum((data$value - make_alpha(data, x_data, beta, alpha, lambda_d, sigma, spline_basis, random_intercept))^2)
}




## when delta is zero make sure this equals the other alpha loss
# probably will be faster if you don't call this each time, like when you only need to calculate loss function
alpha_loss_re = function(data, x_data, spline_basis, basis_time, beta, alpha, delta_coefs, C){
	
	exp_surface = make_surface(basis_time = basis_time, alpha)
	s = attributes(exp_surface)$s
	#ds = c(diff(s)[1], diff(s))
	ds = c(0, diff(s))
	Kt = dim(spline_basis)[2]
	
	Dstar = make_thetaStar(make_surface(basis_time = basis_time, alpha), spline_basis)
  deltaStar = Dstar %*% delta_coefs
	trial_ids = unique(data$trial)
	
	rss = rep(NA, length(trial_ids))
	for(i in 1:length(trial_ids)){
		subject_df = filter(data, trial == trial_ids[i])
		ystar_i = unique(subject_df$Yo) * exp(-alpha * subject_df$time)
		xstar_iB = colSums(t(exp_surface) * rowSums(x_data[[i]] * beta) * ds)
		deltastar_i = deltaStar[, i]
		
		rss[i] = crossprod(subject_df$value) - 
			2 * crossprod(subject_df$value, ystar_i + xstar_iB + deltastar_i) +
			crossprod(ystar_i) + 2 * crossprod(ystar_i, xstar_iB + deltastar_i) + 
			crossprod(xstar_iB) + 2 * crossprod(xstar_iB, deltastar_i) +
			sum(diag(crossprod(Dstar) %*% C)) + 
			crossprod(deltastar_i)
		
	}
	return(sum(rss))
}
