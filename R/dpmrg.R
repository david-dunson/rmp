dpmrg <-
function(ydis, k, nrep, nb, alpha, alpha_r= FALSE, mu0=mean(ydis), kap=var(ydis), 
	atau, btau, a_a=1, b_a=1, lb=NULL, ub=NULL, print = 1, ndisplay = nrep/4,
	plot = FALSE, pdfwrite = FALSE, ...)
{
	if(is.null(lb)) lb <- max(0,min(ydis)-10)
	if(is.null(ub)) ub <- max(ydis)+10
	n <- length(ydis)
	v.alpha <- rep(alpha/k,k)
	print <-  as.integer(table(factor(print, levels=1:6)))
	start.loop <- Sys.time()			
	res <- .C("rmg"	,
			as.integer(c(n,k,ub+1,nrep,alpha_r)),
			as.double(c(mu0, kap, atau, btau,a_a,b_a)),
			as.double(ydis),
			as.integer(c(print,ndisplay)), 
			alpha=as.double(rep(alpha, nrep)),
			mu=as.double(rep(mean(ydis),k*nrep)), 
			tau=as.double(rep(1/var(ydis),k*nrep)), 
			pi=as.double(rep(1/k,k*nrep)), 
			probability=as.double(rep(0,(ub+1)*nrep)),
			enne=as.double(rep(1,k*nrep)),
			PACKAGE="rmp"
			)
	end.loop <- Sys.time()
	mu<-matrix(res$mu,byrow=T,nrow=nrep,ncol=k)
	tau<-matrix(res$tau,byrow=T,nrow=nrep,ncol=k)
	pi<-matrix(res$pi,byrow=T,nrow=nrep,ncol=k)
	enne<-matrix(res$enne,byrow=T,nrow=nrep,ncol=k)	
	cmf<-matrix(res$probability,byrow=T,nrow=nrep,ncol=ub+1)
	pmf<-cmf-cbind(rep(0,nrep),cmf[,1:(dim(cmf)[2]-1)])

	post.pmf<-apply(pmf[(nb+1):nrep,],2,mean)
	sorted.p<-apply(pmf[(nb+1):nrep,],2,sort)
	pmf2.5 <-sorted.p[(0.025*(nrep-nb)),]
	pmf97.5 <-sorted.p[(0.975*(nrep-nb)),]
	
	sorted.n<-t(apply(enne[(nb+1):nrep,],1,sort,decreasing=T))
	n.mean <- apply(sorted.n,2,mean)
	
	arezero = function(x) sum(x==0)
	groups = k - mean( apply(sorted.n,1,arezero))

	post.mu<-apply(mu[(nb+1):nrep,],2,mean)
	post.tau<-apply(tau[(nb+1):nrep,],2,mean)
	post.pi<-apply(pi[(nb+1):nrep,],2,mean)
	post.alpha<-mean(res$alpha[(nb+1):nrep])

	empirical.pmf <- table(factor(ydis, levels=lb:ub))/length(ydis)

	if(plot){
	plot(lb:ub,empirical.pmf,type='h', main="", ylab="pmf", xlab="estimated pmf (blue) and empirical pmf (black)") 
	points(lb:ub+.2,post.pmf[lb:ub+1], col=4, type='h')
	}
	
	if(pdfwrite){
		cat("\nWriting pdf files with traceplots ...\n")
		filename<-paste("traceplots ",date(),".pdf",sep="")
		pdf(filename)
		par(mfrow=c(3,3))
		plot(0,main ="TRACEPLOTS FOR MU -->")  
		for(i in 1:k) 	plot(res$mu,type='l')
		plot(0,main ="TRACEPLOTS FOR TAU -->")  
		for(i in 1:k)  plot(res$tau,type='l')
		plot(0,main ="TRACEPLOTS FOR PI -->") 
		for(i in 1:k)  plot(res$pi,type='l') 
		par(mfrow=c(1,1))
		dev.off()

		if(plot){
		filename<-paste("Estim pmf ",date(),".pdf",sep="")
		pdf(filename)
		plot(lb:ub,empirical.pmf,type='h', main="", ylab="pmf", xlab="estimated pmf (blue) and empirical pmf (black)") 
		points(lb:ub+.2,post.pmf[lb:ub+1],col=4,type='h')
		dev.off()
		}
		cat(paste("Pdf files with traceplots and posterior summaries written in ",getwd(),"\n\n",sep=""))
	}
    out<-structure(
		list(name = "DP mixture of rounded Gaussians", mcmc = list(nrep = nrep, nb = nb, time = (end.loop - start.loop)),  
		  pmf = list(post.pmf=post.pmf[lb:ub+1], empirical=empirical.pmf, lower.95 = pmf2.5[lb:ub+1], upper.95 = pmf97.5[lb:ub+1], domain = lb:ub),
		   parameters = list(post.mu=post.mu, post.tau=post.tau, post.pi=post.pi, post.alpha=post.alpha), 
  	        clustering = list(nmean=n.mean,groups=groups)), class  = "rmpobject")
	cat("\n")
    invisible(out)
	}
