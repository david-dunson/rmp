print.rmpobject <- function(x, ...) 
{
cat(x$name,"\n")
cat(x$mcmc$nrep, "iterations completed in ", as.double(x$mcmc$time), attr(x$mcmc$time, "units"), "\n")
}
#----------------------------------------
summary.rmpobject <- function(object, ...) 
{
cat("------------------------------\n",
	object$name, "\n------------------------------\n",
   "Average number of clusters in the mixture: ", object$clustering$groups,
	"\nPosterior pmf and 95% pointwise credible interval: \n", sep="")
print(data.frame(lower = round(object$pmf$lower.95,4), pmf=round(object$pmf$post.pmf,4), upper = round(object$pmf$upper.95,4), row.names=object$pmf$domain))
cat("-----------------------\n")
}

