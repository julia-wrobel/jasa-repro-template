
# right now only step functions are incorporated, no intercept
simulate_flode = function(N = 50, D = 100, sigma_y0 = 50, sigma = 0.5, seed = 546474, alpha = 2,
													t_max = 1, rand_int = FALSE, lambda0 = 40, 
													lambda_fixed = TRUE, constant_beta = FALSE,
													x = c("step", "periodic")){
	
	set.seed(seed)
	Kt = 10
	
	# define grid and basis
	# intercept only included if x = "periodic"
	grid = seq(0, t_max, length.out = D)
	
	
	if(constant_beta){
		beta_coefs_x = rep(50, length.out = Kt)
		beta_coefs_0 = c(0, 0, seq(0.5, 42, length.out = 4), 40,  40, 40, 40)
	}else{
		beta_coefs_x = 1.5 * c(2, 2, 30, 28, 25, 18, 15, 10, 7, 7)
		beta_coefs_0 = c(-5, 5, seq(10, 45, length.out = 5),  40, 40, 40)
	}
	
	
	beta_coefs_x2 = seq(20, 10, length.out = Kt)
	beta_basis = bs(grid, df = Kt, intercept = TRUE)
	beta_0 = beta_basis %*% beta_coefs_0
	beta_x = beta_basis %*% beta_coefs_x
	beta_x2 = beta_basis %*% beta_coefs_x2
	
	# get P and forcing functions
	P = length(x) + 1
	xmat = xmat2 = matrix(NA, nrow = N, ncol = D)
	jump_points = runif(N, quantile(grid, .1), quantile(grid, .9))
	#sine_shift = runif(N, -0.3, 0.2)
	sine_shift = runif(N, -0.8, 1.4)
	sine_scale = rnorm(N, 1, 0.5)
	sqrt_scale = runif(N, 2, 5)

	# define other elements
	Y0 = rnorm(N, 0, sqrt(sigma_y0))
	epsilon = matrix(rnorm(N*D, 0, sqrt(sigma)), nrow = N, ncol = D)
	Y = matrix(NA, N, D)
	colnames(Y) = paste0("time_", 1:D)
	intercept_surface = make_sim_surface(alpha = alpha, beta = beta_0, D = D, t_max = t_max)
	surface = make_sim_surface(alpha = alpha, beta = beta_x, D = D, t_max = t_max)
	surface2 = make_sim_surface(alpha = alpha, beta = beta_x2, D = D, t_max = t_max)
	#ds = c(diff(grid)[1], diff(grid))
	ds = c(0, diff(grid))
	
	#lambda = rep(c(seq(10, 2, length.out = Kt/2), seq(5, 15, length.out = Kt/2)), each = N)
	if(lambda_fixed){
		lambda = rep(lambda0, Kt * N)
	}else{
		lambda = rep(c(rep(0.5, Kt/2), seq(100, 50, length.out = Kt/2)), each = N)
	}
	
	re_coef = matrix(rnorm(N * Kt, 0, sqrt(lambda)), nrow = N, ncol = Kt)
	re_mat = re_coef %*% t(beta_basis)
	intercept_surface = t(t(intercept_surface) * ds)
	re_mat_data = matrix(NA, N, D)	
	xstar_beta = matrix(NA, N, D)
	
	for(i in 1:N){
		if(x == "step"){
			xmat[i,] = x_f(grid, jump_points[i])
		}else if(x == "periodic"){
			xmat[i,] = sine_scale[i] * sin(pi*seq(0, 1, length.out = D) + 1 + sine_shift[i])
		}else if(x == "multiple"){
			xmat[i,] = sine_scale[i] * sin(pi*seq(0, 1, length.out = D) + 1 + sine_shift[i])
			xmat2[i,] = sin(1*pi*seq(0, 1, length.out = D) + sine_shift[i])
			# http://www.math.ttu.edu/~gilliam/ttu/s10/m3351_s10/sincos_ints_fs.pdf
		}
		
		
		forcing_surface = t(t(surface) * xmat[i,] * ds)
		forcing_surface2 = t(t(surface2) * xmat2[i,] * ds)
		
		if(x == "multiple"){
			Y[i,] = Y0[i] * exp(-alpha * grid) + rowSums(forcing_surface) + 
				rowSums(forcing_surface2) + epsilon[i,]	
			
			if(rand_int){
				surface_re_i = make_sim_surface(alpha = alpha, beta = re_mat[i,], D = D, t_max = t_max)
				surface_re_i = t(t(surface_re_i) * ds)
				
				Y[i,] = Y0[i] * exp(-alpha * grid) + rowSums(forcing_surface) + rowSums(forcing_surface2) + 
					rowSums(surface_re_i) + epsilon[i,]	
			}
			
		}else{
			xstar_beta[i, ] = rowSums(forcing_surface)
			Y[i,] = Y0[i] * exp(-alpha * grid) + rowSums(intercept_surface) + 
				rowSums(forcing_surface) + epsilon[i,]	
			
			if(rand_int){
				surface_re_i = make_sim_surface(alpha = alpha, beta = re_mat[i,], D = D, t_max = t_max)
				surface_re_i = t(t(surface_re_i) * ds)
				re_mat_data[i, ] = rowSums(surface_re_i)
				
				Y[i,] = Y0[i] * exp(-alpha * grid) + rowSums(intercept_surface) + rowSums(forcing_surface) + 
					rowSums(surface_re_i) + epsilon[i,]	
			}
		}
		
	} # end for loop
	
	
	simulated_data = data.frame(
		trial = rep(1:N, each = D),
		y0 = rep(Y0, each = D),
		time = rep(grid, N),
		value = as.vector(t(Y)),
		x = as.vector(t(xmat)),
		x2 = as.vector(t(xmat2)),
		re = as.vector(t(re_mat)),
		reStar = as.vector(t(re_mat_data)),
		xStarbeta = as.vector(t(xstar_beta))
	)
	
	coef_df = data.frame(
		time = grid,
		beta0 = beta_0,
		beta0_dataScale = rowSums(intercept_surface),
		beta1 = beta_x,
		beta2 = beta_x2,
		trial = 0
	)
	
	list(data = simulated_data, alpha = alpha, epsilon = epsilon,
			 y0 = Y0, coefficient_fns = coef_df, surface = surface, surface2 = surface2,
			 re_mat_data = re_mat_data)
	
} # end function


# returns a matrix which is a surface where the columns are S and the rows are T
# equals zero when s > t
make_sim_surface = function(alpha, beta, D, t_max){
	t = s = seq(0, t_max, length.out = D)
	surface = matrix(NA, D, D)
	colnames(surface) = paste0("s", s)
	rownames(surface) = paste0("t", t)
	
	for(time_index in 1:D){
		surface[,time_index] = exp(-alpha * (t-s[time_index]) ) * beta[time_index] * (s[time_index] <= t)
	}
	
	surface
}

x_f = function(grid, jump_point){
	ifelse(grid < jump_point, 0, 1)
}


