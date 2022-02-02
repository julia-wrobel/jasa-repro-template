make_pffr = function(data, covar = "x", random_int = TRUE){
	N = length(unique(data$trial))
	D = length(filter(data, trial == 1)$time)
	
	y_fhist = matrix(data$value, nrow = N, ncol = D, byrow = TRUE)
	x_fhist = matrix(data[[covar]], nrow = N, ncol = D, byrow = TRUE)
	
	y0 = data %>%
		group_by(trial) %>%
		summarize(y0 = first(y0)) %>% ungroup() %>% pull(y0)
	
	fhist_df = data.frame(trial = factor(1:N), y0 =  y0)
	fhist_df$y = y_fhist
	fhist_df$x = x_fhist
	t = s = seq(0, 1, length.out = D)
	
	if(random_int){
		mod = pffr(y ~ s(trial, bs = "re") + ff(x, xind = s, limits = "s<t"), 
							 yind = t, fhist_df)
	}else{
		mod = pffr(y ~ y0 + ff(x, xind = s, limits = "s<t"), yind = t, fhist_df)
	}
	mod
}