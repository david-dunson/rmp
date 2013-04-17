/*============================================================================*/
// Function performing univariate probability mass function estimation
// with Dirichlet process mixture of Poissons 
// C code by Antonio Canale <antonio.canale@unito.it> 
/*============================================================================*/
#include<R.h>
#include<Rmath.h>
void dpmpois(int *ipar, double *alpha, double *hyperpar, double *y,
		int *print, double *lambda, double *pi, double *probability, double *enne)

{
	int n = ipar[0] , 
	    k = ipar[1] , 
	    grid = ipar[2] , 
	    nrep = ipar[3] , 
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
	    S[n*nrep];

	double 
		a = hyperpar[0],
		b = hyperpar[1],
//		newalpha[k],
//		probs[k],
	       ZERO = 0 , 
		rand,
		mass,
		cum_p,
		*newalpha,
		*probs,
		*ysum,
		*nh,
		*apost,
		*bpost; 

	probs = (double *) R_alloc(k, sizeof(double));
	newalpha = (double *) R_alloc(k, sizeof(double));
	ysum = (double *) R_alloc(k, sizeof(double));
	nh = (double *) R_alloc(k, sizeof(double));
	apost = (double *) R_alloc(k, sizeof(double));
	bpost = (double *) R_alloc(k, sizeof(double));

/*
	double ysum[k];
	double nh[k];
	double apost[k];
	double bpost[k];
*/

/* Initialization of the state */
for ( i = 0 ; i < (n*nrep) ; i++ ) S[i] = 1 ; 
GetRNGstate();
if(print_ite & !(print[1] | print[2] | print[3] | print[4]))  Rprintf("\nIterations: ");

/* Main Loop*/
for ( i = 1 ; i < nrep ; i++ ) 
{
R_CheckUserInterrupt(); 
if((i % every) == 0) flag = 1;
else flag = 0;
if(flag & print_ite & !(print[1] | print[2] | print[3] | print[4]))  Rprintf("%d/%d - ",i, nrep);
if(flag & print_ite & (print[1] | print[2] | print[3] | print[4]))  Rprintf("\nIteration %d over %d \n",i,nrep);

if(flag & print_multinomial){ 
Rprintf("################################################\n");
Rprintf(" Step 2: Update S from multinomial distribution \n");
Rprintf("################################################\n");
}

for(j = 0; j < n; j++){
		R_CheckUserInterrupt(); 
	mass = 0;
	for(h = 0; h < k; h++){
		probs[h] =  pi[(i-1)*k + h] * dpois(y[j], 1/lambda[(i-1)*k + h],0);
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
if(flag & print_posterior){
Rprintf("#######################################################################################\n");
Rprintf("Step 3 and 4: Update lambda and pi from the conditional posterior \n");
Rprintf("#######################################################################################\n");
}
for(h = 0; h < k; h++){
	ysum[h] = ZERO;
	nh[h] = ZERO;
	}

for(h = 0; h < k; h++){
	for(j = 0; j < n; j++){
		if(S[i*n + j] == (h+1)){
			nh[h] += 1;
			ysum[h] += y[j];
			}
	apost[h] = a + ysum[h];
	bpost[h] = (1/b) / (nh[h]/b + 1);
	newalpha[h] = alpha[h] + nh[h];
	enne[(i*k)+h] = nh[h];
	lambda[(i*k)+h] = 1/rgamma(apost[h],bpost[h]);
	}
		
	
	
if(flag & print_posterior)	Rprintf("%f %f %f  \n", nh[h],ysum[h],lambda[(i*k)+h] );
}
if(flag & print_posterior)	Rprintf("############################################### \n");

double res_parziale[k] , somma=0;
	for(h=0 ; h< k ; h++){	
		res_parziale[h]=rgamma(newalpha[h],1);
		somma +=res_parziale[h];
		}
	for(h=0 ; h< k ; h++){	
		pi[(i*k)+h]=res_parziale[h]/somma;
if(flag & print_dirichlet)	Rprintf("New pi_%i=%f \n", (h+1), pi[(i*k)+h]);
		}


if(flag & print_post_prob) Rprintf("\n The probabilities are ", probability[(i*grid)]);
for(l = 0; l < grid; l++)
{
	probability[(i*grid)+l] = 0;
	for(h =0; h < k; h++){
		probability[(i*grid)+l] += pi[(i*k)+h] * dpois(l+1, 1/lambda[(i*k)+h], 0 );
		}
	if(flag & print_post_prob) Rprintf("%f, ", probability[(i*grid)+l]);
}
}	
// end of main loop 
PutRNGstate();
} 
// end of all!
