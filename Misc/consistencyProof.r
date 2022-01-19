n <- function(t) {
	y = mean(t) * norm(floor(t), type = "2") - mean(floor(t)) * (t(t) %*% floor(t))
	return(y)
}

d <- function(t) {
	temp = n(t) + (1 - 1/nrow(t)) * mean(floor(t)) * (t(t) %*% floor(t))
	y = (1/mean(t)) * temp

	return(y)
}
