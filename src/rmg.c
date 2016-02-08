/*============================================================================*/
// Function performing univariate probability mass function estimation
// with rounded mixture of Gaussian kernels as in 
// "Canale, A. and Dunson, D. B. (2011), _Bayesian Kernel Mixtures for Counts_, 
// Journal of American Statistical Association, 106, 1528-1539."
// C code by Antonio Canale <antonio.canale@unito.it> 
// part of the code is taken from Nicola Lunardon <lunardon@stat.unipd.it> 
//
// Version 2.0 - 2015 Gen-Feb (by Antonio Canale)
// Two-parameter Poisson Dirichlet prior added
// Slice sampler added 
// Hyper-parameters for the base measure of the 2PD or DP prior
/*============================================================================*/
#include<R.h>
#include<Rmath.h>
void rmg(int *ipar, double *hyperpar, double *ydis, int *print, double *alpha, 
	double *pdp_par, double *mu0, double *kappa, double *atau, double *btau,
	double *mu, double *tau, double *pi, double *minslice, double *probability, double *enne)

{
	int n = ipar[0] , 
	    k = ipar[1] , 
	    grid = ipar[2] , 
	    nrep = ipar[3] ,
    	    HYPRIOR_STICK = ipar[4] , 
    	    HYPRIOR_BASEM = ipar[5] ,  
    	    SLICE = ipar[6] ,  
	    print_ite = print[0] , 
	    print_dataug = print[1], 
	    print_multinomial = print[2] , 
	    print_posterior = print[3] , 
	    print_dirichlet = print[4] , 
	    print_post_prob = print[5] , 
	    every = print[6] , 
	    flag = 1, 
	    i , j ,  h , l, 
	    take,
	    class = 1,
   	    occcl,
//	    S[n*nrep] , 
	    uptail1,
	    uptail2,
		 *nh,
		 *nh_complement,
		 *S;

	S = (int*) R_alloc(n*nrep, sizeof(int));
	nh = (int*) R_alloc(k, sizeof(int));
	nh_complement = (int*) R_alloc(k, sizeof(int));

	mu0[0] = hyperpar[0]; 
	kappa[0] = hyperpar[1];
	atau[0] = hyperpar[2]; 
	btau[0] = hyperpar[3];

	double 
	//mu0 = hyperpar[0], in Version 1.0 fixed hyperpar
	//kappa = hyperpar[1],
	//atau = hyperpar[2], 
	//btau = hyperpar[3],
	a_a = hyperpar[4],
	b_a = hyperpar[5],
	sigma = pdp_par[0],
	theta = pdp_par[1],	
	eta,
	pieta = 1.0,
       	ZERO = 0 , 
	rand ,
	mass,
	cum_p,
	sumtau_thetaminmu2 = 0,
	sumtau = 0,
	summutau = 0,
	sumlogtau = 0,
	*y,
	*u,
	*slices,
	*pu, //for polya urn
	*tresh,
	*probs,
	*ybar,
	*V, *ONE_min_V,
	*ydev,
	*ahattau,
	*bhattau,
	*kappahat,
	*muhat,
	newtau, newmu;

	y = ( double* ) R_alloc(n, sizeof(double));
	u = ( double* ) R_alloc(n, sizeof(double));
	slices = ( double* ) R_alloc(n, sizeof(double));
	pu = (double*) R_alloc(k, sizeof(double));
	tresh = ( double* ) R_alloc(2*n, sizeof(double));
	probs = ( double* ) R_alloc(k, sizeof(double));
	ybar = ( double* ) R_alloc(k, sizeof(double));
	V = ( double* ) R_alloc(k+1, sizeof(double));
	ONE_min_V = ( double* ) R_alloc(k+1, sizeof(double));
	ydev = ( double* ) R_alloc(k, sizeof(double));
	ahattau = ( double* ) R_alloc(k, sizeof(double));
	bhattau = ( double* ) R_alloc(k, sizeof(double));
	kappahat = ( double* ) R_alloc(k, sizeof(double));
	muhat = ( double* ) R_alloc(k, sizeof(double));

/*

	double ybar[k];
	double V[k+1], ONE_min_V[k];
	double ydev[k];
	int nh[k];
	int nh_complement[k];
	double ahattau[k];
	double bhattau[k];
	double kappahat[k];
	double muhat[k];
*/

/* Initialization of the state */
for ( i = 0 ; i < (n*(nrep)) ; i++ ) S[i]  = 1; 
for ( i = 0 ; i < k ; i++ ) S[i]  = i+1; 
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
/* Initialization of the atoms */
for(h = 0; h < k; h++)
{	
	ybar[h] = ydis[h];
	ydev[h] = 0;	
	ahattau[h] = atau[0] + 1/2;
	bhattau[h] = btau[0] + 0.5 * ( ydev[h] + ( 1/(1+kappa[0]) )* ( (ybar[h] - mu0[0])*(ybar[h] - mu0[0]) ) );
	kappahat[h] = 1/(1/kappa[0] +1);
	muhat[h] = kappahat[h] * ( (1/kappa[0]) * mu0[0] +ybar[h] );
	tau[h] = rgamma(ahattau[h], 1/bhattau[h]);
	mu[h]  = rnorm(muhat[h], sqrt(kappahat[h]/tau[h]) );
}
/* Main Loop*/
for ( i = 1 ; i < nrep ; i++ )
{ 
R_CheckUserInterrupt(); 
if((i % every) == 0) flag = 1;
else flag = 0;
if(flag & print_ite & !(print[1] | print[2] | print[3] | print[4]))  Rprintf("%d/%d - ",i, nrep);
if(flag & print_ite & (print[1] | print[2] | print[3] | print[4]))  Rprintf("\nIteration %d over %d \n",i,nrep);

if(flag & print_dataug)
{
	for(j = 0; j < n; j++) 
	{
		Rprintf("ydis %f Class S = %d mu %f tau %f \n", ydis[j], S[((i-1)*n) +j],  
						 mu[((i-1)*k) + S[((i-1)*n)+j] -1], 
						tau[((i-1)*k) + S[((i-1)*n)+j] -1]) ;
	}
}
/* generate y^*  */
if(flag & print_dataug) Rprintf("###########################\nStep 1: Data augmentation  \n###########################\n");
for ( j = 0 ; j < n ; j++ )		
{
	uptail1 = ( ydis[j]+1 - mu[((i-1)*k) + S[((i-1)*n)+j] -1]) * sqrt( tau[((i-1)*k) + S[((i-1)*n)+j] -1]) > 5;
	uptail2 = ( ydis[j]   - mu[((i-1)*k) + S[((i-1)*n)+j] -1]) * sqrt( tau[((i-1)*k) + S[((i-1)*n)+j] -1]) > 5;

	if( (uptail1 + uptail2)==0 )
	{
		tresh[n+j]   =  pnorm(ydis[j]+1, mu[((i-1)*k) + S[((i-1)*n)+j] -1], 1/sqrt(tau[((i-1)*k) + S[((i-1)*n)+j] -1]), 1, 1 );
		tresh[j]   =    pnorm(ydis[j],   mu[((i-1)*k) + S[((i-1)*n)+j] -1], 1/sqrt(tau[((i-1)*k) + S[((i-1)*n)+j] -1]), 1, 1 );	
		u[j] = runif(tresh[j] , tresh[n+j] );
		y[j] = qnorm(	u[j],              mu[((i-1)*k) + S[((i-1)*n)+j] -1], 1/sqrt(tau[((i-1)*k) + S[((i-1)*n)+j] -1]), 1, 1 );	
	}
	else
	{
		tresh[n+j]   =  pnorm(ydis[j]+1, mu[((i-1)*k) + S[((i-1)*n)+j] -1], 1/sqrt( tau[((i-1)*k) + S[((i-1)*n)+j] -1]), 0, 1 );
		tresh[j]   =    pnorm(ydis[j],   mu[((i-1)*k) + S[((i-1)*n)+j] -1], 1/sqrt( tau[((i-1)*k) + S[((i-1)*n)+j] -1]), 0, 1 );	
		u[j] = runif(tresh[n+j] , tresh[j] );
		y[j] = qnorm(u[j],               mu[((i-1)*k) + S[((i-1)*n)+j] -1], 1/sqrt( tau[((i-1)*k) + S[((i-1)*n)+j] -1]), 0, 1);
	}
if(flag & print_dataug) Rprintf("Up1 %i t1 %e Up2 %i t2 %e u = %e ydis = %f y = %f \n", uptail1,tresh[j],uptail2,tresh[n+j],u[j],ydis[j],y[j]);
}
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
	//for each observation
	mass =0;
	for(h = 0; h < k; h++)
	{
		// For each h in 1...k compute the p(S = h)
		if(pi[(i-1)*k + h]>slices[j])
		{
			if(ydis[j] !=0) 
			{
				probs[h] =  (pnorm(ydis[j]+1, mu[(i-1)*k + h], sqrt(1/tau[(i-1)*k + h]), 1, 0) - 
				             pnorm(ydis[j],   mu[(i-1)*k + h], sqrt(1/tau[(i-1)*k + h]), 1, 0));
			}
			else
			{
				probs[h] =  pnorm(1, mu[(i-1)*k + h], sqrt(1/tau[(i-1)*k + h]), 1, 0);
			
			}
		}
		else probs[h] =0;
		mass += probs[h];
	}

	if(flag & print_multinomial)
	{
	Rprintf("Class from  %d ",S[((i-1)*n)+j]);
	}
	rand = runif(0,1);
	cum_p = 0;
	for (h = 0; h < k; h++)
	{
		probs[h] =  probs[h]/mass;	
		cum_p += probs[h];
		take = (rand < cum_p);
	  	if(take)
		{ 
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
	//if cluster size of subjet j was one, the "new atom" is the same, otherwise
	//for the first empty cluster, allocate a new atom from the base measure
	//and assign the new weight proportional to theta (total mass alpha in the DP)
	if(nh[ S[(i-1)*n + j] - 1 ] == 0)
	{
		if(flag & print_multinomial) Rprintf("it left a cluster empty:");
		pu[S[(i-1)*n + j] - 1] = theta + occcl*sigma;	
	}
	else
	{
	for(h = 0; h < k; h++)
	{
		if(nh[h] == 0)
		{
			if(flag & print_multinomial) Rprintf("1st empty cluster is %i:",h+1);
			pu[h] = theta + occcl*sigma;
			tau[((i-1)*k) + h]  = rgamma(atau[i-1], 1/btau[i-1]);
			 mu[((i-1)*k) + h]  = rnorm(mu0[i-1], sqrt(kappa[i]/tau[((i-1)*k)+h]));
			h = k;
		}
	}
	}
	mass = 0;
	for(h = 0; h < k; h++)
	{
		if(ydis[j] !=0) 
		{
			probs[h] = pu[h]* (pnorm(ydis[j]+1, mu[(i-1)*k + h], sqrt(1/tau[(i-1)*k + h]), 1, 0) - 
			                   pnorm(ydis[j],   mu[(i-1)*k + h], sqrt(1/tau[(i-1)*k + h]), 1, 0));
		}
		else
		{
			probs[h] =  pu[h]* pnorm(1, mu[(i-1)*k + h], sqrt(1/tau[(i-1)*k + h]), 1, 0);
		
		}
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
if(flag & print_dirichlet){
Rprintf("#######################################################\n");
Rprintf("Step 3: Sample V_h from beta and udpdate the weights	\n");
Rprintf("#######################################################\n");
}

for(h = 0; h < k; h++)
{
	ybar[h] = ZERO;
	ydev[h] = ZERO;
	nh[h] = ZERO;
	ahattau[h] = ZERO;
	bhattau[h] = ZERO;
	kappahat[h] = ZERO;
	muhat[h] = ZERO;
}
nh_complement[0] = n;
//deviance, sample cluster means, cluster size and stuff like that
for(h = 0; h < k; h++)
{
	ybar[h] = 0;
	ydev[h] = 0;
	nh[h] = 0;
	for(j = 0; j < n; j++)
	{
		if(S[i*n + j] == (h+1))
		{
			nh[h] += 1;
			ybar[h] += y[j];
			ydev[h] += (y[j]*y[j]) ;
		}
	}
	nh_complement[h] = nh_complement[h-1] - nh[h];		
	if(h==0) nh_complement[h] = n - nh[h];	
}
//number of occupied clusters (used if alpha \sim Ga() in DP mixture and generally if
//hyperpriors for the base measure are assumed
occcl = k;
for(h =0; h<k; h++)
{
	if(nh[h] == 0) occcl = occcl - 1;
}
if(HYPRIOR_STICK)
{
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
if(flag & print_dirichlet) Rprintf("\n");
//-----
if(SLICE)
{
ONE_min_V[0] = 1;
V[k+1] = 1;
for(h = 0; h<k ; h++)
{
	if(flag & print_dirichlet) Rprintf("nh = %i, nh_compl=%i, theta = %f, sigma= %f",nh[h],nh_complement[h],theta,sigma);
	V[h] = rbeta(1-sigma+nh[h],  theta + (h+1)*sigma + nh_complement[h]);
	if(flag & print_dirichlet) Rprintf("V_%i = %f", h+1, V[h]);
	ONE_min_V[h+1] = ONE_min_V[h] * (1 - V[h]);
	if(flag & print_dirichlet) Rprintf(" - (1 - V_%i) = %f", h+1, ONE_min_V[h]);
	pi[(i*k)+h] = V[h] * ONE_min_V[h];
	if(flag & print_dirichlet) Rprintf(" - pi_%i = %f\n", h+1, pi[(i*k)+h]);
}	
}
// Se Marginal Sampler
else
{
for(h = 0; h<k ; h++)
{
	if(nh[h]>0)
	{
		if(flag & print_dirichlet) Rprintf("nh = %i, nh_compl=%i, theta = %f, sigma= %f",nh[h],nh_complement[h],theta,sigma);
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
if(flag & print_posterior){
Rprintf("##########################################################\n");
Rprintf("Step 4: Update mu and tau from the conditional posterior  \n"); //also known as reshuffling!
Rprintf("##########################################################\n");
}
for(h = 0; h < k; h++)
{	
	if(nh[h] > 1)
	{
		ybar[h] = ybar[h] / nh[h];
		ydev[h] = ydev[h] - nh[h] * (ybar[h]*ybar[h]);
	}
	else
	{
		if(nh[h]==1) ydev[h] = 0;	
//		ydev[h] = 0;
//		ybar[h] = 0;	
	}
	ahattau[h] = atau[i-1] + nh[h]/2;
	bhattau[h] = btau[i-1] + 0.5 * ( ydev[h] + ( nh[h]/(1+kappa[i-1]*nh[h]) )* ( (ybar[h] - mu0[i-1])*(ybar[h] - mu0[i-1]) ) );
	kappahat[h] = 1/(1/kappa[i-1] +nh[h])	;
	muhat[h] = kappahat[h] * ( (1/kappa[i-1]) * mu0[i-1] +nh[h]*ybar[h] );
	tau[(i*k)+h] = rgamma(ahattau[h], 1/bhattau[h]);
	mu[(i*k)+h]  = rnorm(muhat[h], sqrt(kappahat[h]/tau[(i*k)+h]) )	;
	enne[(i*k)+h] = nh[h];
if(flag & print_posterior)	Rprintf("%f %f %f %f %f  \n", nh[h],ybar[h],ydev[h],mu[(i*k)+h],tau[(i*k)+h] );
}


if(HYPRIOR_BASEM)
{
if(flag & print_posterior & HYPRIOR_BASEM){
Rprintf("##########################################################\n");
Rprintf("Step 4b: Update (mu, tau)'s hyperprior hyperparameters    \n");
Rprintf("##########################################################\n");
}
sumtau_thetaminmu2 = 0;
sumtau = 0;
summutau = 0;
sumlogtau = 0;
for(h = 0; h < k; h++)
{	
	if(nh[h] > 1)
	{
		sumtau_thetaminmu2 = sumtau_thetaminmu2 + tau[(i*k)+h]*pow(mu[(i*k)+h]-mu0[i-1],2);
		sumtau = sumtau + tau[(i*k)+h];
		summutau = summutau + mu[(i*k)+h]*tau[(i*k)+h];
		sumlogtau = sumlogtau + log(tau[(i*k)+h]);

	}
}
//kappa[i] = 1/rgamma(1+0.5*occcl,1/(0.5*sumtau_thetaminmu2) );
//mu0[i] = rnorm(summutau/sumtau, sqrt(kappa[i]/sumtau));
btau[i] = rgamma(0.5+occcl*atau[i-1],1/(0.5 + sumtau));
}
if(flag & print_post_prob){
Rprintf("##########################################################\n");
Rprintf("Step 5: Compute posterior PMF    \n");
Rprintf("##########################################################\n");
}

if(flag & print_post_prob) Rprintf("\n The probabilities are ", probability[(i*grid)]);
newtau= rgamma(atau[i], 1/btau[i]);
newmu =rnorm(mu0[i], sqrt(kappa[i]/newtau));
for(l = 0; l < grid; l++){
	probability[(i*grid)+l] = 0;
//NEW}	
//NEW	rand = runif(0,1);
//NEW	cum_p = 0;
//NEW	for (h = 0; h < k; h++)
//NEW	{
//NEW		cum_p += pi[(i*k)+h];
//NEW	  	if((rand < cum_p)) 
//NEW		{ 
//NEW			for(l = 0; l < grid; l++)
//NEW			{
//NEW			probability[(i*grid)+l] = pnorm(l+1, mu[(i*k)+h], sqrt(1/tau[(i*k)+h]), 1, 0 );
//NEW			}
//NEW		h = k;
//NEW		}
//NEW	}
//NEW  	if((rand > cum_p))
//NEW	{
//NEW		newtau= rgamma(atau[i], 1/btau[i]);
//NEW		newmu =rnorm(mu0[i], sqrt(kappa[i]/newtau))	;
//NEW		for(l = 0; l < grid; l++)
//NEW		{
//NEW			probability[(i*grid)+l] = pnorm(l+1, newmu, sqrt(1/newtau), 1, 0 );
//NEW		}
//NEW	}
	for(h =0; h < k; h++){
		probability[(i*grid)+l] += pi[(i*k)+h] * pnorm(l+1, mu[(i*k)+h], sqrt(1/tau[(i*k)+h]), 1, 0 );
		}
	probability[(i*grid)+l] +=  (ONE_min_V[k+1]) * pnorm(l+1, newmu, sqrt(1/newtau), 1, 0 );
if(flag & print_post_prob) Rprintf("%f, ", probability[(i*grid)+l]);
}
}	
// end of main loop
PutRNGstate();
} 
// end of all!



