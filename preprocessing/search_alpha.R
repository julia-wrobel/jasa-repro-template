
# forcing functions are a design matrix, if you want to include an intercept it needs to include one
# forcing functions are a matrix that is part of the input data frame, with 1 for intercept included
search_alpha = function(data, Kt = 20,  forcing_functions = c("int", "x"),
												alpha_min = 2, alpha_max = 14){
	
	# process data to extract forcing functions
	x_data = data %>% select(trial, all_of(forcing_functions))  %>% nest(data = c(-trial)) %>% pull(data)
	
	basis_time = filter(data, trial == 1)$time # assumes grid is even across subjects
	spline_basis = bs(basis_time, df = Kt, intercept = TRUE)
	
	N = length(unique(data$trial))
	D = length(basis_time)
	P = length(forcing_functions)
	
	get_alpha_loss = function(data, alpha0){
		# make yiStar term
		data = data %>% mutate(Yo = y0, 
													 y0_star = Yo * exp(- alpha0 * time))
		
		##################################################
		# find beta
		thetaStar = integrate_alpha(basis_time, spline_basis, alpha0, covar = NULL)
		xstar_ls = lapply(x_data, make_Xstar, alpha = alpha0, spline_basis = spline_basis, basis_time = basis_time,
											thetaStar = thetaStar)
		xstar = do.call(rbind, xstar_ls)
		beta_coefs = solve(crossprod(xstar)) %*% (t(xstar) %*% (data$value - data$y0_star))
		beta_coefs = matrix(beta_coefs, nrow = Kt, ncol = P)
		
		##################################################
		# calculate loss based on current alpha and beta
		sum((data$value - data$y0_star - as.numeric(xstar %*% as.vector(beta_coefs)))^2)
	}
	# apply across grid of values and return results
	
	alpha_grid = seq(alpha_min, alpha_max, by = 2)
	loss = sapply(alpha_grid, get_alpha_loss, data = data)
	best_alpha = alpha_grid[which.min(loss)]
	alpha_grid2 = seq(best_alpha - 1, best_alpha + 1, by = 0.5)
	loss2 = sapply(alpha_grid2, get_alpha_loss, data = data)
	best_alpha = alpha_grid2[which.min(loss2)]
	alpha_grid3 = seq(best_alpha - .5, best_alpha + .5, by = 0.1)
	loss3 = sapply(alpha_grid3, get_alpha_loss, data = data)
	best_alpha = alpha_grid3[which.min(loss3)]
	
	loss_df = tibble(alpha = c(alpha_grid, alpha_grid2, alpha_grid3), loss = c(loss, loss2, loss3)) %>%
		arrange(loss)
	
	list(best_alpha = best_alpha, loss = loss_df)
}
