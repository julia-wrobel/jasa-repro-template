####################################################################
# Julia Wrobel
# December 30, 2021
# 
# This file produces simulations for flode paper.
# Compares flode to functional historical regression
# 
####################################################################

library(splines)
library(tibble)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(mvtnorm)
library(pracma)
library(tictoc)
library(refund)




# get parallelization set up
library(foreach)
library(doParallel)
n_cores = detectCores() - 1
registerDoParallel(n_cores)

###############################################################
## define or source functions used in code below
###############################################################


source(here::here("source", "make_pffr.R"))
source(here::here("source", "simulated_paw.R"))
source(here::here("source", "estimate_flode.R"))
source(here::here("source", "utils.R"))
source(here::here("source", "search_alpha.R"))

###############################################################
## set simulation design elements 
###############################################################
N_vec = c(100)
D_vec = c(50)
lambda_d_vec = c(50)
alpha_vec = c(2, 4, 6, 8, 10, 12, 0.1, 0.5, 1)

seed_start = 1000 
Kt = 20
lambda_b = 100
N_iter = 30
miter = 200

params = expand.grid(seed_start = seed_start, 
										 N_vec = N_vec,
										 D_vec = D_vec, 
										 lambda_d_vec = lambda_d_vec,
										 alpha_vec = alpha_vec)


## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())


###############################################################
## start simulation code
###############################################################


## need to parallelize this part
## actually, parallelize the iterations?
for(scenario in 7:9){
	
	results_mat = matrix(NA, nrow = 1, ncol = 19)
	colnames(results_mat) = c("iter","scenario", "seed", "alpha_est", "lambda_d_est", "sse_flode", "lambda_b_est", 
													"maxiter", "N", "D", "alpha_true", "lambda_d_true", "surfaceErr_flode", "time_flode", 
													"sse_fhist", "surfaceErr_fhist", "time_fhist", "sigma_flode", "sigma_fhist")

	results_mat[, 2] = scenario

	
	###############################################################
	## set simulation design elements 
	###############################################################
	N = params$N_vec[scenario]
	D = params$D_vec[scenario]
	SEED.START = params$seed_start[scenario]
	alpha0 = params$alpha_vec[scenario]
	lambda_d0 = params$lambda_d_vec[scenario]
	
	results_mat[, 9] = N
	results_mat[, 10] = D
	results_mat[, 11] = alpha0
	results_mat[, 12] = lambda_d0
	
	results = foreach(iter = 1:N_iter, .combine = 'rbind') %dopar% {
		
		seed.iter = (SEED.START - 1)*N_iter + iter
		results_mat[, 3] = seed.iter
		
		## generate data
		simulated_data = simulate_flode(I = N, 
																		D = D, 
																		sigma_y0 = 5, 
																		sigma = 0.1, 
																		alpha = alpha0, 
																		rand_int = TRUE, 
																		lambda0 = lambda_d0,
																		seed = seed.iter)
		
		
		################################################################################
		## flode 
		dat = simulated_data$data %>% mutate(int = 1)
		
		set.seed(seed.iter)
		
		tic(quiet = TRUE)
		
		# grid search to initialize alpha
		alpha_start = search_alpha(data = dat, Kt = Kt)$best_alpha
		
		# run flode
		flode_results = estimate_flode(dat, 
																	 Kt = Kt, 
																	 alpha0 = alpha_start,  
																	 max_iter = miter, 
																	 forcing_functions = c("int", "x"),
																	 tol =  0.00001, 
																	 lambda_b = lambda_b,
																	 lambda_d = lambda_d0,
																	 random_effects = TRUE,
																	 initial_position = TRUE,
																	 y0_vec = simulated_data$y0)
		time_flode = toc()
		
		results_mat[, 1] = iter
		
		results_mat[, 4] = flode_results$alpha
		results_mat[, 5] = flode_results$lambda_d
		results_mat[, 6] = flode_results$sse
		results_mat[, 7] = flode_results$lambda_b
		results_mat[, 8] = flode_results$maxiter
		results_mat[, 18] = flode_results$sigma
		
		################################################################################
		## get surface errors and computation time for flode
		
		# flode surface error
		results_mat[, 13] = sum((simulated_data$surface - flode_results$surface[[2]])^2)
		
		# get computation time
		results_mat[, 14] = time_flode$toc -time_flode$tic
		
		################################################################################
		## fhist using random intercept
		tic(quiet = TRUE)
		fhist_results_re = make_pffr(dat, random_int = TRUE)
		time_fhist = toc()
		
		fhist_surface_re = coef(fhist_results_re, n1 = D, n2 = D)$smterms$"ff(x,s)"$coef %>%
			select(s = x.smat, t = x.tmat, value = value) %>% as_tibble() 
		
		## interpolate surface to get on same grid as flode
		interp_obj = list(x = unique(fhist_surface_re$s), y = unique(fhist_surface_re$t),
											z = matrix(fhist_surface_re$value, D, D, byrow = TRUE))
		
		t = s = seq(0, 1, length.out = D)
		loc_mat = list(x = s, y = t)
		surface_re = fields::interp.surface.grid(interp_obj, loc_mat)$z
		surface_re[D,] = surface_re[D-1,] # carry last value forward so not NA
		surface_re = surface_re * lower.tri(surface_re, diag = TRUE)# zero out upper triangle of the matrix
		
		fhist_intercept_re = coef(fhist_results_re, n1 = D, n2 = D)$smterms$"Intercept(t)"$coef %>%
			select(t = t.vec, value = value) %>% as_tibble() 
		
		
		results_mat[, 15] = sum((fhist_results_re$residuals)^2)
		
		################################################################################
		## get surface and computation time for fhist
		results_mat[, 16] = sum((simulated_data$surface - surface_re)^2)
		
		# get computation time
		results_mat[, 17] = time_fhist$toc -time_fhist$tic
		
		# get sigma
		results_mat[, 19] = fhist_results_re$scale
		################################################################################

		
		results_mat
	} # end N_iter foreach loop

	filename = paste0("/Users/juliawrobel/Documents/projects/2018/201803_Hantman_Trajectories/results/", "20211230" ,"/", "flodeVSfhist_", scenario, ".RDA")
	save(results,
			 file = filename)
} # end scenario loop 

###############################################################
## end sim
###############################################################


