dpmpoiss <-
function(y, k, nrep, nb, alpha=1, a, b, lb=NULL, ub=NULL, print = 1, 
	ndisplay = nrep/4, plot = FALSE, pdfwrite = FALSE, ... )
{
	if(is.null(lb)) lb <- max(0,min(y)-10)
	if(is.null(ub)) ub <- max(y)+10
	n <- length(y)
	v.alpha <- rep(alpha/k,k)
	print <-  as.integer(table(factor(print, levels=1:5)))
	start.loop <- Sys.time()	
	res <- .C("dpmpois",
				as.integer(c(n,k,ub+1,nrep)),
				as.double(v.alpha), 
				as.double(c(a,b)),
				as.double(y),
				as.integer(c(print,ndisplay)), 
				lambda=as.double(rep(1,k*nrep)), 
				pi=as.double(rep(1/k,k*nrep)), 
				probability=as.double(rep(0,(ub+1)*nrep)),
				enne=as.double(rep(1,k*nrep)),
				PACKAGE="rmp"
				)
	end.loop <- Sys.time()
	lambda<-matrix(res$lambda,byrow=T,nrow=nrep,ncol=k)
	pi<-matrix(res$pi,byrow=T,nrow=nrep,ncol=k)
	enne<-matrix(res$enne,byrow=T,nrow=nrep,ncol=k)	
	pmf<-matrix(res$probability,byrow=T,nrow=nrep,ncol=ub+1)
	
	post.pmf<-apply(pmf[(nb+1):nrep,],2,mean)
	sorted.p<-apply(pmf[(nb+1):nrep,],2,sort)
	pmf2.5 <-sorted.p[(0.025*(nrep-nb)),]
	pmf97.5 <-sorted.p[(0.975*(nrep-nb)),]
	
	sorted.n<-t(apply(enne[(nb+1):nrep,],1,sort,decreasing=T))
	n.mean <- apply(sorted.n,2,mean)
	arezero = function(x) sum(x==0)
	groups = k - mean( apply(sorted.n,1,arezero))

	post.lambda<-apply(lambda[(nb+1):nrep,],2,mean)
	post.pi<-apply(pi[(nb+1):nrep,],2,mean)

	empirical.pmf <- table(factor(y, levels=lb:ub))/length(y)
  
	if(plot){
	plot(lb:ub,empirical.pmf,type='h', main="", ylab="pmf", xlab="estimated pmf (red) and empirical pmf (black)") 
	points(lb:ub+.4,post.pmf[lb:ub+1],col=2,type='h')
	}
	
	if(pdfwrite){
		cat("\nWriting pdf files with traceplots ...\n")
		filename<-paste("traceplots ",date(),".pdf",sep="")
		pdf(filename)
		par(mfrow=c(3,3))
		plot(0,main ="TRACEPLOTS FOR LAMBDA -->")  
		for(i in 1:k) 	plot(res$lambda,type='l') 
		plot(0,main ="TRACEPLOTS FOR PI -->")   
		for(i in 1:k)  plot(res$pi,type='l') 
		par(mfrow=c(1,1))
		dev.off()

		filename<-paste("Estim pmf ",date(),".pdf",sep="")
		pdf(filename)
		plot(lb:ub,empirical.pmf,type='h',main="", ylab="pmf", xlab="estimated pmf (red) and empirical pmf (black)") 
		points(lb:ub+.4,post.pmf[lb:ub+1],col=2, type='h')
		dev.off()
		cat(paste("Pdf files with traceplots and posterior summaries written in ",getwd(),"\n\n",sep=""))
	}

    out<-structure(
		list(name = "DP mixture of Poisson kernels", mcmc = list(nrep = nrep, nb = nb, time = (end.loop - start.loop)),  
		  pmf = list(post.pmf=post.pmf[lb:ub+1], empirical=empirical.pmf, lower.95 = pmf2.5[lb:ub+1], upper.95 = pmf97.5[lb:ub+1], domain = lb:ub),
		  parameters = list(post.lambda=post.lambda, post.pi=post.pi),
  	        clustering = list(nmean=n.mean,groups=groups)), class  = "rmpobject")
	cat("\n")
    invisible(out)
    	}

