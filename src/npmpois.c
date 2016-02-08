/*============================================================================*/
// Function performing univariate probability mass function estimation
// with nonparametric (DP or 2PD) mixture of Poissons 
// C code by Antonio Canale <antonio.canale@unito.it> 
/*============================================================================*/
#include<R.h>
#include<Rmath.h>
void npmpois(int *ipar, double *alpha, double *pdp_par, double *hyperpar, double *y,
		int *print, double *lambda, double *pi, double *ach, double *bch,
		double *minslice,  double *probability, double *enne)
{
	int n = ipar[0] , 
	    k = ipar[1] , 
	    grid = ipar[2] , 
	    nrep = ipar[3] ,
	    mixing_hyperprior = ipar[4],  
	    HYPRIOR_BASEM  = ipar[5],  
    	    SLICE = ipar[6] ,  
	    print_ite = print[0] , 
	    print_multinomial = print[1] , 
	    print_posterior = print[2] , 
	    print_dirichlet = print[3] , 
	    print_post_prob = print[4] , 
	    every = print[5] , 
	    flag = 0,
	    i , j ,  h , l, 
	    take,
   	    class = 1,
   	    occcl,
	    *S,
		*nh, *nh_complement;

	S = (int*) R_alloc(n*nrep, sizeof(int));
	nh = (int*) R_alloc(k, sizeof(int));
	nh_complement = (int*) R_alloc(k, sizeof(int));


	double 
		a = hyperpar[0],
		b = hyperpar[1],
		hyperpriortype = hyperpar[2], // 1 for gamma(a,b) [mean a/b], 2 for normal(a,b) with b sd 
	       ZERO = 0 , 
		eta = 0.5,
		a_a = hyperpar[3],
		b_a = hyperpar[4],
		sigma = pdp_par[0],
		theta = pdp_par[1],	
   	       pieta = 1.0,
		rand,
		mass,
		cum_p,
		*slices,
		*pu, //for polya urn weights
		*V, *ONE_min_V,
		*probs,
		*ysum,
		*apost,
		*bpost,
		mhprob,
		phistar,
		newlambda; 
	probs = (double *) R_alloc(k, sizeof(double));
	ysum = (double *) R_alloc(k, sizeof(double));
	V = ( double* ) R_alloc(k+1, sizeof(double));
	ONE_min_V = ( double* ) R_alloc(k+1, sizeof(double));
	apost = (double *) R_alloc(k, sizeof(double));
	bpost = (double *) R_alloc(k, sizeof(double));
	slices = ( double* ) R_alloc(n, sizeof(double));
	pu = (double*) R_alloc(k, sizeof(double));
/*
	double ysum[k];
	double nh[k];
	double apost[k];
	double bpost[k];
*/

/* Initialization of the state */
for ( i = 0 ; i < (n*nrep) ; i++ ) S[i] = 1 ; 
enne[0] = n-k+1;
for ( i = 1 ; i < k ; i++ ) enne[i] = 1;
occcl = k;
for(h = 0; h < k; h++)
{
	nh[h] = 0;
	for(j = 0; j < n; j++)
	{
		if(S[j] == (h+1))
		{
			nh[h] += 1;
		}
	}
}
GetRNGstate();
if(print_ite & !(print[1] | print[2] | print[3] | print[4]))  Rprintf("\nIterations: ");

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
/*
hyper-hyper-prior sugli iperparametri del Kernel
*/
double meanlambda, sumlamminach2;
ach[0] = a;
bch[0] = b;
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


/* Main Loop*/
for ( i = 1 ; i < nrep ; i++ ) 
{
R_CheckUserInterrupt(); 
if((i % every) == 0) flag = 1;
else flag = 0;
if(flag & print_ite & !(print[1] | print[2] | print[3] | print[4]))  Rprintf("%d/%d - ",i, nrep);
if(flag & print_ite & (print[1] | print[2] | print[3] | print[4]))  Rprintf("\nIteration %d over %d \n",i,nrep);
if(SLICE)
{
if(flag & print_multinomial){ 
Rprintf("################################################\n");
Rprintf(" Step 2.0: First do the slice sampler \n");
Rprintf("################################################\n");
}

for(j = 0; j < n; j++)
{
	slices[j] = runif(0,1)*pi[(i-1)*k + S[(i-1)*n + j]-1];
}
minslice[i] = slices[0]; 
for(j = 1; j < n; j++)
{
	if(slices[j]<minslice[i]) minslice[i] = slices[j];
}
if(flag & print_multinomial){ 
Rprintf("################################################\n");
Rprintf(" Step 2: Update S from multinomial distribution \n");
Rprintf("################################################\n");
}
for(j = 0; j < n; j++)
{
		R_CheckUserInterrupt(); 
	mass = 0;
	for(h = 0; h < k; h++)
	{
		if(pi[(i-1)*k + h]>slices[j])
		{
			probs[h] =  dpois(y[j], lambda[(i-1)*k + h], 0);
		}
		else probs[h] = 0;
		mass += probs[h];
	}
	if(flag & print_multinomial){
	Rprintf("The total mass is = %e \n",mass);
	Rprintf("Class from  %d ",S[((i-1)*n)+j]);
	}
	rand = runif(0,1);
	cum_p = 0;
	for (h = 0; h < k; h++){
		probs[h] =  probs[h]/mass;	
		cum_p += probs[h];
		take = (rand < cum_p);
	  	if(take) { 
			class = h+1;
			h = k;
			}
	}
	if(flag & print_multinomial)	Rprintf(" to %d \n",class);
	S[i*n + j] = class;

}	
}
else
{
if(flag & print_multinomial){ 
Rprintf("################################################\n");
Rprintf(" Step 2.0: Polya-Urn sampler \n");
Rprintf("################################################\n");
}
for(j = 0; j < n; j++)
{	
	if(flag & print_multinomial)
	{
	Rprintf("Subj %i:",j+1);
	}
	//take the index of the j-th subject
	h = S[(i-1)*n + j] - 1;
	//decrease the related cluster size
	nh[h] = nh[h] -1;	
	//recalculate the number of occupied clusters
	occcl = k;
	for(h =0; h<k; h++)
	{
		if(nh[h] == 0) occcl = occcl - 1;
	}	
	//compute the weihgs. everywhere proportional to the cluster size
	for(h = 0; h < k; h++)
	{
		pu[h] = nh[h] - sigma;
		if(nh[h] == 0) pu[h] = 0;
	}
	//if cluster size of subjet j was one, place the "new atom" in the same location, ...
	if(nh[ S[(i-1)*n + j] - 1 ] == 0)
	{
		if(flag & print_multinomial) Rprintf("it left a cluster empty:");
		pu[S[(i-1)*n + j] - 1] = theta + occcl*sigma;	
	}
	//... otherwise for the first empty cluster location, allocate a "new atom" from the base measure
	//and assign the new weight proportional to theta + K*sigma (total mass alpha in the DP)
	else
	{
	for(h = 0; h < k; h++)
	{
		if(nh[h] == 0)
		{
			if(flag & print_multinomial) Rprintf("1st empty cluster is %i:",h+1);
			pu[h] = theta + occcl*sigma;
			if(hyperpriortype==1) 
			{
			lambda[((i-1)*k) + h] = rgamma(a,1/b);
			}
			if(hyperpriortype==2)
			{	
			lambda[((i-1)*k) + h] = exp(rnorm(ach[i-1], sqrt(bch[i-1])));
			}
			h = k;
		}
	}
	}
	mass = 0;
	for(h = 0; h < k; h++)
	{
		probs[h] = pu[h]*dpois(y[j], lambda[(i-1)*k + h], 0);
		mass += probs[h];
	}
	if(flag & print_multinomial)
	{
	Rprintf("   -- Class from  %d ",S[((i-1)*n)+j]);
	}
	rand = runif(0,1);
	if(flag & print_multinomial) Rprintf(" - %f : p, ",rand);
	cum_p = 0;
	for (h = 0; h < k; h++)
	{
		probs[h] =  probs[h]/mass;	
		cum_p += probs[h];
		take = (rand < cum_p);
	if(flag & print_multinomial & nh[h] > 0) Rprintf("%f,",probs[h]);
	  	if(take)
		{ 
			class = h+1;
			h = k;
		}

	}
	if(flag & print_multinomial)	Rprintf("  - to %d (w= %f = n_h(%i)*f(y)",class, probs[class-1], nh[class-1]);
	S[i*n + j] = class;
	nh[class-1] = nh[class-1] + 1;
	if(flag & print_multinomial) Rprintf(" (now class %i has %i subjects)",class, nh[class-1]);
	if(flag & print_multinomial)	Rprintf("\n");
}
}
if(flag & print_posterior){
Rprintf("#######################################################################################\n");
Rprintf("Step 3 and 4: Update lambda and pi from the conditional posterior \n");
Rprintf("#######################################################################################\n");
}
for(h = 0; h < k; h++)
{
	ysum[h] = ZERO;
	nh[h] = ZERO;
}
nh_complement[0] = n;
for(h = 0; h < k; h++)
{
	for(j = 0; j < n; j++)
	{
		if(S[i*n + j] == (h+1))
		{
			nh[h] += 1;
			ysum[h] += y[j];
		}
	}
	apost[h] = a + ysum[h];
	bpost[h] = b + nh[h];
	enne[(i*k)+h] = nh[h];
	if(hyperpriortype==1) 
	{
		lambda[(i*k)+h] = rgamma(apost[h],1/bpost[h]);
	}
	if(hyperpriortype==2)
	{	
		phistar = rnorm(log(lambda[((i-1)*k)+h]), bch[i-1]);
		mhprob = exp( (phistar - log(lambda[((i-1)*k)+h])) * ysum[h] - 
			(exp(phistar) - lambda[((i-1)*k)+h]) * nh[h] - 
			0.5*(phistar                 -ach[i-1])*(phistar                 -ach[i-1])/(bch[i-1]*bch[i-1]) +
			0.5*(log(lambda[((i-1)*k)+h])-ach[i-1])*(log(lambda[((i-1)*k)+h])-ach[i-1])/(bch[i-1]*bch[i-1]) );
		if(mhprob>=1) lambda[(i*k)+h] = exp(phistar);
		else
		{
			rand = runif(0,1);
			if(rand<mhprob) lambda[(i*k)+h] = exp(phistar);
			else{ lambda[(i*k)+h] = lambda[((i-1)*k)+h];}
		}
	}	
	if(h==0) nh_complement[h] = n - nh[h];	
	else nh_complement[h] = nh_complement[h-1] - nh[h];		
	if(flag & print_posterior & (nh[h]>0))	Rprintf("h=%i, nh=%i, m=%f, lambda=%f\n", h+1, nh[h],ysum[h],lambda[(i*k)+h] );
}
if(flag & print_posterior)	Rprintf("############################################### \n");


if(mixing_hyperprior){
	occcl = k;
	// sample alpha precision of the DP
	for(h =0; h<k; h++){
		if(flag & print_dirichlet) Rprintf("%i, ", nh[h]);
		if(nh[h] == 0) occcl = occcl -1;
		}
	if(flag & print_dirichlet) Rprintf("Occ cluster = %i", occcl);		
	eta = rbeta(alpha[i-1]+1,n);
	if(flag & print_dirichlet) Rprintf("\n eta = %f", eta);
	pieta = (a_a + occcl - 1) / (n*b_a - n*log(eta) + a_a + occcl - 1);
	if(flag & print_dirichlet) Rprintf("\n pi_eta = %f", pieta);
	
	rand = runif(0,1);
	
	if(rand<pieta)
	{
	alpha[i] = rgamma(a_a + occcl, 1/(b_a - log(eta)));
	}
	else
	{
	alpha[i] = rgamma(a_a + occcl - 1, 1/(b_a - log(eta)));	
	}
	if(flag & print_dirichlet) Rprintf("alpha = %f", alpha[i]);
	theta = alpha[i];
}
if(SLICE)
{
ONE_min_V[0] = 1;
V[k+1] = 1;
for(h = 0; h < k ; h++){
		if(flag & print_dirichlet) Rprintf("nh = %i, nh_compl=%i, theta = %f, sigma= %f",nh[h],nh_complement[h],theta,sigma);
		V[h] = rbeta(1-sigma+nh[h],  theta + (h+1)*sigma + nh_complement[h]);
		if(flag & print_dirichlet) Rprintf("V_%i = %f", h+1, V[h]);
		ONE_min_V[h+1] = ONE_min_V[h] * (1 - V[h]);
		if(flag & print_dirichlet) Rprintf(" - (1 - V_%i) = %f", h+1, ONE_min_V[h]);
		pi[(i*k)+h] = V[h] * ONE_min_V[h];
		if(flag & print_dirichlet) Rprintf(" - pi_%i = %f\n", h+1, pi[(i*k)+h]);
		}
}
// Se invece usiamo il marginal Sampler
else
{
for(h = 0; h<k ; h++)
{
	if(nh[h]>0)
	{
		if(flag & print_dirichlet) Rprintf("nh = %i, theta = %f, sigma= %f",nh[h],theta,sigma);
		pi[(i*k)+h] = (nh[h]-sigma) /(n+theta);
		if(flag & print_dirichlet) Rprintf(" - pi_%i = %f\n", h+1, pi[(i*k)+h]);
	}
	else
	{
		pi[(i*k)+h] = 0;
	}
}
(ONE_min_V[k+1]) = (theta + occcl*sigma)/(n + theta);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
/*
Hyperprior on base measure's parameters {only if exp(lambda) ~ N(.,.) }
*/
if(HYPRIOR_BASEM & (hyperpriortype==2))
{
occcl = k;
meanlambda =0;
sumlamminach2 = 0;
for(h =0; h<k; h++)
{
	if(nh[h] == 0)
	{
		occcl = occcl -1;
	}
	else
	{
		meanlambda = meanlambda + log(lambda[(i*k)+h]);
		sumlamminach2 = sumlamminach2 + pow(log(lambda[(i*k)+h]) - ach[i-1],2);
	}
}
meanlambda = meanlambda/occcl;
eta = rnorm(log(bch[i-1]),0.1);
mhprob = exp( eta*(1-occcl)           - 0.5*(exp(-2*eta)*sumlamminach2) 
	   -log(bch[i-1])*(1-occcl) + 0.5*(exp(-2*log(bch[i-1]))*sumlamminach2)  )	;
			if(mhprob>=1) bch[i] = exp(eta);
			else
			{
				rand = runif(0,1);
				if(rand<mhprob) bch[i] = exp(eta);
				else{ bch[i] = bch[i-1];}
			}
ach[i] = rnorm(meanlambda, sqrt(bch[i]/occcl)); //3.58;//
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

if(flag & print_post_prob) Rprintf("\n The probabilities are ", probability[(i*grid)]);
for(l = 0; l < grid; l++)
{
	probability[(i*grid)+l] = 0;
	for(h =0; h < k; h++)
	{
		probability[(i*grid)+l] += pi[(i*k)+h] * dpois(l, lambda[(i*k)+h], 0);
	}
	newlambda = exp(rnorm(ach[i], sqrt(bch[i])));
	probability[(i*grid)+l] +=  (ONE_min_V[k+1]) * dpois(l, newlambda, 0);
	if(flag & print_post_prob) Rprintf("%f, ", probability[(i*grid)+l]);
}

}	
// end of main loop 
PutRNGstate();
} 
// end of all!
