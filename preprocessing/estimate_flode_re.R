
# forcing functions are a design matrix, if you want to include an intercept it needs to include one
# forcing functions are a matrix that is part of the input data frame, with 1 for intercept included
estimate_flode_re = function(data, Kt = 10, alpha0 = 6,  
										 forcing_functions = c(""), max_iter = 100,
										 tol =  0.00001, seed = 1988, 
										 update_y0 = FALSE, mle = TRUE, lambda_b = .1, lambda_d = .1){
	set.seed(seed)
	# process data to extract forcing functions
	x_data = data %>% select(trial, forcing_functions)  %>% nest(-trial) %>% pull(data)
	
	basis_time = filter(data, trial == 1)$time # assumes grid is even across subjects
	spline_basis = bs(basis_time, df = Kt, intercept = TRUE)
	
	N = length(unique(data$trial))
	D = length(basis_time)
	P = length(forcing_functions)
	
	# make yiStar term
	data = data %>% mutate(Yo = y0, y0_star = Yo * exp(- alpha0 * time))
	
	##################################################
	# initialize alpha and beta
	alpha_cur = alpha0
	xstar_ls = lapply(x_data, make_Xstar, alpha = alpha_cur, spline_basis = spline_basis, time = basis_time)
	xstar = do.call(rbind, xstar_ls)
	beta_coefs = solve(crossprod(xstar)) %*% (t(xstar) %*% (data$value - data$y0_star))
	beta_coefs = matrix(beta_coefs, nrow = Kt, ncol = P)
	#beta_coefs = matrix(rnorm(Kt * P), nrow = Kt, ncol = P)
	beta_cur = spline_basis %*% beta_coefs
	
	##################################################
	# initialize variance parameters 
	sigma_cur = 1
	lambda_d_cur = lambda_d
	lambda_b_cur = rep(lambda_b, P)
	
	##################################################
	# define penalty matrix
	diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
	Pen = t(spline_basis) %*% t(diff2) %*% diff2 %*% spline_basis
	Pen_b = rep(1/lambda_b_cur, each = Kt) * Matrix::bdiag(rep(list(Pen), P))
	
	
	##################################################
	# initialize coefficients for random intercepts
	delta_coefs = matrix(rnorm(Kt *  N, 0, sqrt(lambda_d_cur)), nrow = Kt, ncol = N)
	#delta_coefs = matrix(0, nrow = Kt, ncol = N)
	
	Dstar = make_thetaStar(make_surface(basis_time = basis_time, alpha_cur), spline_basis)
	deltaStar = as.vector(Dstar %*% delta_coefs)
	data = data %>% mutate(delta =  as.vector(spline_basis %*% delta_coefs))
	C = solve(1/lambda_d_cur * Pen + crossprod(Dstar)/sigma_cur)
	
	##################################################
	# iteratively calculate alpha and beta
	iter = 1
	error_cur = 100
	y0_mat = matrix(NA, nrow = N, ncol = max_iter + 1)
	colnames(y0_mat) = paste0("iter_", 0:max_iter)
	y0_mat[, 1] = unique(data$Yo)
	sse_vec = alpha_vec = lambda_d_vec = rep(NA, max_iter)
	sse_mat = matrix(NA, nrow = max_iter, ncol = 6)
	colnames(sse_mat) = c("alpha", "beta", "lambdas", "sigma", "delta", "y0")
	while(iter < max_iter && error_cur > tol){
		
		##################################################
		# update alpha
		alpha_nls = optim(alpha_cur, alpha_loss_re, data = data, x_data = x_data,
											spline_basis = spline_basis, basis_time = basis_time,
											beta = beta_cur, 
											delta_coefs = delta_coefs, C = C,
											method = "Brent", lower = 1, upper = 15)
		
		error_cur = abs(alpha_cur - alpha_nls$par) 
		alpha_cur = alpha_vec[iter] = alpha_nls$par
		
		sse_mat[iter, 1] = alpha_nls$value
		
		# update terms that contain alpha
		Dstar = make_thetaStar(make_surface(basis_time = basis_time, alpha_cur), spline_basis)
		data = data %>% mutate(y0_star = Yo * exp(- alpha_cur * time))
		C = solve(1/lambda_d_cur * Pen + crossprod(Dstar)/sigma_cur)
		
		##################################################
		# update beta
		xstar_ls = lapply(x_data, make_Xstar, alpha = alpha_cur, spline_basis = spline_basis, time = basis_time)
		xstar = do.call(rbind, xstar_ls)
		beta_coefs = solve(crossprod(xstar) + sigma_cur * Pen_b) %*% (t(xstar) %*% (data$value - data$y0_star - deltaStar))
		beta_coefs  = matrix(beta_coefs, nrow = Kt, ncol = P)
		beta_cur = spline_basis %*% beta_coefs
		
		# sse
		sse_mat[iter, 2] = alpha_loss_re(data, x_data, spline_basis, basis_time, beta_cur, alpha_cur, delta_coefs, C)
		
		# update lambda_b and Pen_b
		lambda_b_cur = diag(crossprod(beta_coefs, Pen %*% beta_coefs))/Kt 
		Pen_b = rep(1/lambda_b_cur, each = Kt) * Matrix::bdiag(rep(list(Pen),P))
		
		# update lambda_d
		trial_ids = unique(data$trial)
		lambda_di = rep(NA, length(trial_ids))
		trace_PC = sum(diag(Pen %*% C))
		for(i in 1:length(trial_ids)){
			lambda_di[i] = trace_PC + crossprod(delta_coefs[,i], Pen %*% delta_coefs[,i])
		}

		lambda_d_cur = lambda_d_vec[iter] = sum(lambda_di)/(N * Kt)
		
		# sse
		sse_mat[iter, 3] = alpha_loss_re(data, x_data, spline_basis, basis_time, beta_cur, alpha_cur, delta_coefs, C)
		
		##################################################
		# update sigma
		RSS = alpha_loss_re(data, x_data, spline_basis, basis_time, beta_cur, alpha_cur, delta_coefs, C)
		sigma_cur = as.numeric(RSS)/(N * D)
		
		C = solve(1/lambda_d_cur * Pen + crossprod(Dstar)/sigma_cur)
		
		# sse
		sse_mat[iter, 4] = alpha_loss_re(data, x_data, spline_basis, basis_time, beta_cur, alpha_cur, delta_coefs, C)
		
		##################################################
		# update delta
		for(i in 1:length(trial_ids)){
			trial_df = filter(data, trial == trial_ids[i])
			
			delta_coefs[, i] = (C %*% t(Dstar)) %*% (trial_df$value - trial_df$y0_star - 
																							 	xstar_ls[[i]] %*% as.vector(beta_coefs))/sigma_cur
		}
		
		deltaStar = as.vector(Dstar %*% delta_coefs)
		data = data %>% mutate(delta = as.vector(spline_basis %*% delta_coefs),
													 deltaStar = deltaStar)

		
		# sse
		sse_mat[iter, 5] = alpha_loss_re(data, x_data, spline_basis, basis_time, beta_cur, alpha_cur, delta_coefs, C)
		
		##################################################
		# update initial position
		if(update_y0){
			data = data %>% 
				mutate(yhat = y0_star + as.numeric(xstar %*% as.vector(beta_coefs))) %>%
				as_tibble() %>% nest(-trial) %>%
				mutate(data = map2(data, xstar_ls, make_y0, 
													 alpha = alpha_cur, 
													 beta_coefs = as.vector(beta_coefs),
													 slice_num = 1,
													 mle = mle,
													 use_re = TRUE),
							 ) %>%  
				unnest() %>%
				# update y0_star
				mutate(y0_star = Yo * exp(- alpha_cur * time))
			
			y0_mat[, iter + 1] = unique(data$Yo)
			
		}
		
		
		
		##################################################
		sse = sse_vec[iter] = sse_mat[iter, 6] = alpha_loss_re(data, x_data, spline_basis, basis_time, beta_cur, alpha_cur, delta_coefs, C)
		
		message(paste0("current iteration: ", iter))
		message(paste0("current alpha: ", alpha_cur))
		message(paste0("current error: ", error_cur))
		message(paste0("current sse: ", sse))
		message(paste0("current sigma: ", sigma_cur))
		message(paste0("current lambda_d: ", lambda_d_cur))
		message(paste0("current lambda_b: ", lambda_b_cur))
		
		############################################################
		#
		# change error calculation to involve more parameters?
		#error_cur = sse
		iter = iter + 1
	} # end while loop
	##################################################

	data = data %>% mutate(yhat = y0_star + as.numeric(xstar %*% as.vector(beta_coefs)) + deltaStar)

	############################################################
	# get surface for each forcing function
	surfaces = as.list(rep(NA, P))
	for(p in 1:P){
		beta_p = beta_cur[, p]
		surfaces[[p]] = make_surface(basis_time, alpha_cur, beta = beta_p) 
	}

	colnames(beta_cur) = paste0("beta", 0:(P-1))
	beta_cur = as_tibble(beta_cur) %>% mutate(time = basis_time)
	
	list(data = data, 
			 beta = beta_cur, 
			 alpha = as.numeric(alpha_cur),
			 sigma = sigma_cur,
			 lambda_d = lambda_d_cur,
			 lambda_b = lambda_b_cur,
			 delta_coefs = delta_coefs,
			 loss = sse, 
			 surface = surfaces, 
			 maxiter = iter,
			 y0_mat = y0_mat,
			 sse_vec = sse_vec,
			 sse_mat = sse_mat, 
			 alpha_vec = alpha_vec,
			 lambda_d_vec = lambda_d_vec)
}
