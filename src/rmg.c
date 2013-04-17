/*============================================================================*/
// Function performing univariate probability mass function estimation
// with rounded mixture of Gaussian kernels as in 
// "Canale, A. and Dunson, D. B. (2011), _Bayesian Kernel Mixtures for Counts_, 
// Journal of American Statistical Association, 106, 1528-1539."
// C code by Antonio Canale <antonio.canale@unito.it> 
// part of the code is taken from Nicola Lunardon <nicola.lunardon@econ.units.it> 
/*============================================================================*/
#include<R.h>
#include<Rmath.h>
void rmg(int *ipar, double *hyperpar, double *ydis,
		int *print, double *alpha, double *mu, double *tau, double *pi, double *probability, double *enne)

{
	int n = ipar[0] , 
	    k = ipar[1] , 
	    grid = ipar[2] , 
	    nrep = ipar[3] ,
    	    alpha_r = ipar[4] ,  
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

	double 
		mu0 = hyperpar[0],
		kappa = hyperpar[1],
		atau = hyperpar[2], 
		btau = hyperpar[3],
		a_a = hyperpar[4],
		b_a = hyperpar[5],
		eta,
//   	       newalpha[k],
   	       pieta = 1.0,
/*
	       y[n]  ,
	       u[n]  ,
	       tresh[2*n] ,
*/
//		probs[k],
	       ZERO = 0 , 
		rand ,
		mass,
		cum_p,

		*y,
		*u,
		*tresh,
		*newalpha,
		*probs,
		*ybar,
		*V, *ONE_min_V,
		*ydev,
		*ahattau,
		*bhattau,
		*kappahat,
		*muhat;


	newalpha = ( double* ) R_alloc(k, sizeof(double));
	y = ( double* ) R_alloc(n, sizeof(double));
	u = ( double* ) R_alloc(n, sizeof(double));
	tresh = ( double* ) R_alloc(2*n, sizeof(double));
	probs = ( double* ) R_alloc(k, sizeof(double));
	ybar = ( double* ) R_alloc(k, sizeof(double));
	V = ( double* ) R_alloc(k+1, sizeof(double));
	ONE_min_V = ( double* ) R_alloc(k, sizeof(double));
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
for ( i = 0 ; i < (n*nrep) ; i++ ) S[i]  = 1 ; 
GetRNGstate();
if(print_ite & !(print[1] | print[2] | print[3] | print[4]))  Rprintf("\nIterations: ");

/* Main Loop*/
for ( i = 1 ; i < nrep ; i++ ) { 
R_CheckUserInterrupt(); 
if((i % every) == 0) flag = 1;
else flag = 0;
if(flag & print_ite & !(print[1] | print[2] | print[3] | print[4]))  Rprintf("%d/%d - ",i, nrep);
if(flag & print_ite & (print[1] | print[2] | print[3] | print[4]))  Rprintf("\nIteration %d over %d \n",i,nrep);

if(flag & print_dataug)
{
for(j = 0; j < n; j++) {
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

	if( (uptail1 + uptail2)==0 ){
		tresh[n+j]   =  pnorm(ydis[j]+1, mu[((i-1)*k) + S[((i-1)*n)+j] -1], 1/sqrt(tau[((i-1)*k) + S[((i-1)*n)+j] -1]), 1, 1 );
		tresh[j]   =    pnorm(ydis[j],   mu[((i-1)*k) + S[((i-1)*n)+j] -1], 1/sqrt(tau[((i-1)*k) + S[((i-1)*n)+j] -1]), 1, 1 );	
		u[j] = runif(tresh[j] , tresh[n+j] );
		y[j] = qnorm(	u[j],              mu[((i-1)*k) + S[((i-1)*n)+j] -1], 1/sqrt(tau[((i-1)*k) + S[((i-1)*n)+j] -1]), 1, 1 );	
	}
	else{
		tresh[n+j]   =  pnorm(ydis[j]+1, mu[((i-1)*k) + S[((i-1)*n)+j] -1], 1/sqrt( tau[((i-1)*k) + S[((i-1)*n)+j] -1]), 0, 1 );
		tresh[j]   =    pnorm(ydis[j],   mu[((i-1)*k) + S[((i-1)*n)+j] -1], 1/sqrt( tau[((i-1)*k) + S[((i-1)*n)+j] -1]), 0, 1 );	
		u[j] = runif(tresh[n+j] , tresh[j] );
		y[j] = qnorm(u[j],               mu[((i-1)*k) + S[((i-1)*n)+j] -1], 1/sqrt( tau[((i-1)*k) + S[((i-1)*n)+j] -1]), 0, 1);
	}
if(flag & print_dataug) Rprintf("Up1 %i t1 %e Up2 %i t2 %e u = %e ydis = %f y = %f \n", uptail1,tresh[j],uptail2,tresh[n+j],u[j],ydis[j],y[j]);
}

if(flag & print_multinomial){ 
Rprintf("################################################\n");
Rprintf(" Step 2: Update S from multinomial distribution \n");
Rprintf("################################################\n");
}
for(j = 0; j < n; j++){	
	//for each observation
	mass =0;
	for(h = 0; h < k; h++){
	// For each h in 1...k compute the p(S = h)
		if(ydis[j] !=0) {
			probs[h] =  pi[(i-1)*k + h] * (pnorm(ydis[j]+1, mu[(i-1)*k + h], sqrt(1/tau[(i-1)*k + h]), 1, 0) - 
				                          pnorm(ydis[j],   mu[(i-1)*k + h], sqrt(1/tau[(i-1)*k + h]), 1, 0));
			}
		else{
			probs[h] =  pi[(i-1)*k + h] * pnorm(1, mu[(i-1)*k + h], sqrt(1/tau[(i-1)*k + h]), 1, 0);
			}
		mass += probs[h];
		}
	if(flag & print_multinomial){

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

if(flag & print_dirichlet){
Rprintf("#######################################################\n");
Rprintf("Step 3: Sample V_h from beta and udpdate the weights	\n");
Rprintf("#######################################################\n");
}

for(h = 0; h < k; h++){
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
for(h = 0; h < k; h++){
	ybar[h] = 0;
	ydev[h] = 0;
	nh[h] = 0;
	for(j = 0; j < n; j++){
		if(S[i*n + j] == (h+1)){
			nh[h] += 1;
			ybar[h] += y[j];
			ydev[h] += (y[j]*y[j]) ;
			}
		}
	nh_complement[h] = nh_complement[h-1] - nh[h];		
	if(h==0) nh_complement[h] = n - nh[h];	
	}

if(alpha_r){
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
}

if(flag & print_dirichlet) Rprintf("\n");
ONE_min_V[0] = 1;
V[k+1] = 1;
for(h = 0; h<k ; h++){
		if(flag & print_dirichlet) Rprintf("nh = %i, %f+%i = %f   ",nh[h], alpha[i], nh_complement[h], alpha[i] + nh_complement[h]);
		V[h] = rbeta(1+nh[h],  alpha[i] + nh_complement[h]);
		if(flag & print_dirichlet) Rprintf("V_%i = %f", h+1, V[h]);
		ONE_min_V[h+1] = ONE_min_V[h] * (1 - V[h]);
		if(flag & print_dirichlet) Rprintf(" - (1 - V_%i) = %f", h+1, ONE_min_V[h]);
		pi[(i*k)+h] = V[h] * ONE_min_V[h];
		if(flag & print_dirichlet) Rprintf(" - pi_%i = %f\n", h+1, pi[(i*k)+h]);
		}	

if(flag & print_posterior){
Rprintf("##########################################################\n");
Rprintf("Step 4: Update mu and tau from the conditional posterior  \n");
Rprintf("##########################################################\n");
}

for(h = 0; h < k; h++){	
	if(nh[h] > 1){
	ybar[h] = ybar[h] / nh[h];
	ydev[h] = ydev[h] - nh[h] * (ybar[h]*ybar[h]);
	}
	else{
		if(nh[h]==1) ydev[h] = 0;	
	}
	ahattau[h] = atau + nh[h]/2;
	bhattau[h] = btau + 0.5 * ( ydev[h] + ( nh[h]/(1+kappa*nh[h]) )* ( (ybar[h] - mu0)*(ybar[h] - mu0) ) );
	kappahat[h] = 1/(1/kappa +nh[h])	;
	muhat[h] = kappahat[h] * ( (1/kappa) * mu0 +nh[h]*ybar[h] );
	tau[(i*k)+h] = rgamma(ahattau[h], 1/bhattau[h]);
	mu[(i*k)+h]  = rnorm(muhat[h], sqrt(kappahat[h]/tau[(i*k)+h]) )	;
	newalpha[h] = alpha[i] + nh[h];
	enne[(i*k)+h] = nh[h];
if(flag & print_posterior)	Rprintf("%f %f %f %f %f  \n", nh[h],ybar[h],ydev[h],mu[(i*k)+h],tau[(i*k)+h] );
}
if(flag & print_posterior)	Rprintf("############################################### \n");


if(flag & print_post_prob) Rprintf("\n The probabilities are ", probability[(i*grid)]);
for(l = 0; l < grid; l++){
	probability[(i*grid)+l] = 0;
	for(h =0; h < k; h++){
		probability[(i*grid)+l] += pi[(i*k)+h] * pnorm(l+1, mu[(i*k)+h], sqrt(1/tau[(i*k)+h]), 1, 0 );
		}
if(flag & print_post_prob) Rprintf("%f, ", probability[(i*grid)+l]);
}
}	
// end of main loop
PutRNGstate();
} 
// end of all!



