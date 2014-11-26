/*--------------------------------------------------------------------------------------*/
/* linesearch.c																			*/
/*	Implements various line search techniques.  These are specific to the functions		*/
/* that are to be minimized.   This is because each function we are minimizing takes	*/
/* different parameters.   It would be a pain to make a general function that took a	*/
/* variable length parameter list.   I didn't really feel like doing this....			*/
/*--------------------------------------------------------------------------------------*/

#ifndef LINESEARCH_C_
#define LINESEARCH_C_

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "linesearch.h"


/*--------------------------------------------------------------*/
/*	Purpose: bisection method for calculating theta				*/
/*  Returns: Optimal value of THETA								*/
/*  Args   :												    */
/*		(double) A_in,B_in : initial interval					*/
/*		(double) epsilon : quits when range < epsilon			*/
/*		(double) SIG_in	 : SIGMA								*/
/*		(double) m_in	 : number of edges						*/
/*		(double*) f_in	 : scaled edge congestion				*/
/*		(int)	max_iters: maximum number of iterations for		*/
/*			bisection method.   This is a failsafe				*/
/*																*/
/*	Other Notes:												*/
/*		This function will run to max_iters when epsilon is too */
/*	small.   I believe this is because the solution lies between*/
/*	a two representable numbers.   Thus, we keep iterating		*/
/*	forever.													*/
/*--------------------------------------------------------------*/

double solve_BISECT_THETA(double A_in, double B_in, double epsilon, double SIG_in, int m_in, double *f_in, int max_iters) {
	int  i;
	
	double objective(double x, double SIG, int m, double *f) {
		int i;
		double partial;
	
		partial = 0.0;
		for (i = 1; i <= m+1; i++)
			partial += x / (x - f[i]);
		partial *= SIG/(6.0*(m+1));
		partial -= 1.0;
	
		return(partial);
	}
	
	double a,b,f_a,f_b,ret;
	
	
	a = A_in;
	b = B_in;
	
	f_a = objective(a,SIG_in,m_in,f_in);
	f_b = objective(b,SIG_in,m_in,f_in);

	i = 0;
		
	if ( f_a*f_b >= 0 ) {
		printf("\t\tInvalid starting points to bisection method\n");
		
		printf("\t\ta = %e , b = %e\n",a,b);
		printf("\t\tf_a = %e , f_b = %e\n",f_a , f_b);
		
		while ( f_a*f_b >= 0 ) {
			printf("\t\t\tEnter new values for a and b\n");
			scanf("%lf %lf",&a,&b);
			printf("\t\t\tECHO a = %f, b = %f\n",a,b);
			
			f_a = objective(a,SIG_in,m_in,f_in);
			f_b = objective(b,SIG_in,m_in,f_in);
			printf("\t\t\tf_a = %f, f_b = %f\n",f_a,f_b);
		}
	}
	
	while (( (b-a)>epsilon )&&(i < max_iters)) {

		ret = objective((a+b)/2.0,SIG_in,m_in,f_in);
		
		
		if ((ret*f_b) < 0) {
			a = (a+b)/2.0;
			f_a = objective(a,SIG_in,m_in,f_in);
		}
		
		if ((f_a*ret) < 0) {
			b = (a+b)/2.0;
			f_b = objective(b,SIG_in,m_in,f_in);
		}
		i++;
	}
	return((a+b)/2.0);
}

double solve_BISECT_THETA_ECMM(double A_in, double B_in, double epsilon, double SIG_in, int m_in, double *f_in, int max_iters) {
	int  i;
	
	double objective(double x, double SIG, int m, double *f) {
		int i;
		double partial;
		
		partial = 0.0;
		for (i = 1; i <= m; i++)
			partial += x / (x - f[i]);
		partial *= SIG/(6.0*((double)m));
		partial -= 1.0;
		
		return(partial);
	}
	
	double a,b,f_a,f_b,ret;
	
	
	a = A_in;
	b = B_in;
	
	f_a = objective(a,SIG_in,m_in,f_in);
	f_b = objective(b,SIG_in,m_in,f_in);
	
	i = 0;
	
	if ( f_a*f_b >= 0 ) {
		printf("\t\tInvalid starting points to bisection method\n");
		
		printf("\t\ta = %e , b = %e\n",a,b);
		printf("\t\tf_a = %e , f_b = %e\n",f_a , f_b);
		
		while ( f_a*f_b >= 0 ) {
			printf("\t\t\tEnter new values for a and b\n");
			scanf("%lf %lf",&a,&b);
			printf("\t\t\tECHO a = %f, b = %f\n",a,b);
			
			f_a = objective(a,SIG_in,m_in,f_in);
			f_b = objective(b,SIG_in,m_in,f_in);
			printf("\t\t\tf_a = %f, f_b = %f\n",f_a,f_b);
		}
	}
	
	while (( (b-a)>epsilon )&&(i < max_iters)) {
		
		ret = objective((a+b)/2.0,SIG_in,m_in,f_in);
		
		
		if ((ret*f_b) < 0) {
			a = (a+b)/2.0;
			f_a = objective(a,SIG_in,m_in,f_in);
		}
		
		if ((f_a*ret) < 0) {
			b = (a+b)/2.0;
			f_b = objective(b,SIG_in,m_in,f_in);
		}
		i++;
	}
	return((a+b)/2.0);
}


/*--------------------------------------------------------------*/
/*	Purpose: Calculate value of potential function				*/
/*  Returns: potential function value							*/
/*  Args   :												    */
/*		(double) TAU : value of TAU								*/
/*		(double*) f	 : scaled edge congestion					*/
/*		(double*) f_bar : scaled edge congestion of current		*/
/*					iteration									*/
/*		(double*) tmp_f	: workspace								*/
/*		(double) SIG : value of SIGMA							*/
/*		(double) EPS : value of EPSILON							*/
/*		(double) m   : number of edges (graph constraints)		*/
/*--------------------------------------------------------------*/
double pot_fun(double TAU , double *f, double *f_bar, double *tmp_f, double SIG, double EPS, double m) {
	double	theta,tmp_LAM,partial,pot_val;
	int		i;
	double	LB,UB,Tol;
	
	for (i = 1; i <= m+1; i++)
		tmp_f[i] = f[i];
	
	for (i = 1; i <= m+1; i++)
		tmp_f[i] = (1.0 - TAU)*f[i] + TAU*f_bar[i];
	
	tmp_LAM = -DBL_MAX;
	for (i = 1; i <= m+1; i++)
		if (tmp_f[i] > tmp_LAM)
			tmp_LAM = tmp_f[i];
	
	
	/* solve for theta */
	LB = tmp_LAM / (1.0 - (SIG/6.0)*(1.0/(m+1.0)));
	UB = tmp_LAM / (1.0 - (SIG/6.0));
	Tol = (EPS*EPS)/(m+1.0);							
	theta = solve_BISECT_THETA(LB, UB , Tol, SIG,m,tmp_f,MAX_BISECT_ITERS);
	
	/* calculate potential value */
	pot_val = log(theta);
	partial = 0.0;
	for (i = 1; i <= m+1; i++)
		partial += log(theta - tmp_f[i]);
	partial *= ((SIG/6.0)/(m+1.0));
	pot_val = pot_val - partial;
	
	return pot_val;
}

/*--------------------------------------------------------------*/
/*	Purpose: Uses secant method to compute an approximation of	*/
/*		the derivitive of the potential function				*/
/*  Returns: ~ potfun'											*/
/*  Args   :												    */
/*		(double) TAU : value of TAU								*/
/*		(double*) f	 : scaled edge congestion					*/
/*		(double*) f_bar : scaled edge congestion of current		*/
/*					iteration									*/
/*		(double*) tmp_f	: workspace								*/
/*		(double) SIG : value of SIGMA							*/
/*		(double) EPS : value of EPSILON							*/
/*		(double) m   : number of edges (graph constraints)		*/
/*		(double) delta : differential							*/
/*--------------------------------------------------------------*/
double pot_fun_deriv(double TAU , double *f, double *f_bar, double *tmp_f, double SIG, double EPS, double m, double delta) {
	return (  (pot_fun(TAU+delta ,f, f_bar, tmp_f, SIG, EPS, m) - pot_fun(TAU ,f, f_bar, tmp_f, SIG, EPS, m))/delta );
}


/*--------------------------------------------------------------*/
/*	Purpose: bisection method for calculating optimal value of	*/
/*		TAU for minimizing potential function					*/
/*  Returns: Optimal value of TAU								*/
/*  Args   :												    */
/*		(double) A_in,B_in : initial interval					*/
/*		(double) epsilon : quits when range < epsilon			*/
/*		(double) SIG_in	 : SIGMA								*/
/*		(double) m_in	 : number of edges						*/
/*		(double*) f_bar  : scaled edge congestion of current	*/
/*					iteration									*/
/*		(double*) f	 : scaled edge congestion					*/
/*		(double*) tmp_f : workspace								*/
/*		(double) EPS : value of EPSILON							*/
/*		(int)	max_iters: maximum number of iterations for		*/
/*			bisection method.   This is a failsafe				*/
/*																*/
/*--------------------------------------------------------------*/
double solve_BISECT_potfun(double A_in, double B_in, double epsilon, double delta, double SIG, int m, double *f_bar, double *f, double *tmp_f, double EPS, int max_iters) {
  int  i;	
  double a,b,f_a,f_b,ret;
  	  
	  
	
	a = A_in;
	b = B_in;
	
	f_a = pot_fun_deriv(a ,f, f_bar, tmp_f, SIG, EPS, m, delta);
	f_b = pot_fun_deriv(b ,f, f_bar, tmp_f, SIG, EPS, m, delta);
	
	i = 0;

	if ( f_a*f_b >= 0 ) {
		return 0.0;
	}
	
	while (( (b-a)>epsilon )&&(i < max_iters)) {
		
		ret = pot_fun_deriv( (a+b)/2.0 ,f, f_bar, tmp_f, SIG, EPS, m, delta);
		
		
		if ((ret*f_b) < 0) {
			a = (a+b)/2.0;
			f_a = pot_fun_deriv( a , f, f_bar, tmp_f, SIG, EPS, m, delta);
		}
		
		if ((f_a*ret) < 0) {
			b = (a+b)/2.0;
			f_b = pot_fun_deriv( b , f, f_bar, tmp_f, SIG, EPS, m, delta);
		}
		i++;
	}
	return((a+b)/2.0);
}				

/*--------------------------------------------------------------*/
/*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   */
/*					DEPRICATED									*/
/*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   */
/*--------------------------------------------------------------*/
double secant(double x_0, double x_1, double epsilon, double delta, int max_iters, 
			  double func(double TAU , double *f, double *f_bar, double *tmp_f, double SIG, double EPS, double m),
			  double *f, double *f_bar, double *tmp_f, double SIG, double EPS, double m) {
	int k;
	double	f_im1, f_i;
	double	x_im1, x_i, x_ip1;
	
	double func_prime(double TAU , double *f, double *f_bar, double *tmp_f, double SIG, double EPS, double m, double delta) {
		return (( func(TAU+delta,f,f_bar,tmp_f,SIG,EPS,m) - func(TAU,f,f_bar,tmp_f,SIG,EPS,m) )/ delta);
	}
	
	
	x_im1 = x_0;		f_im1 = func_prime(x_im1,f,f_bar,tmp_f,SIG,EPS,m,delta);
	x_i   = x_1;		f_i   = func_prime(x_i,f,f_bar,tmp_f,SIG,EPS,m,delta);
	k = 1;
	
	x_ip1 = x_i - func_prime(x_i,f,f_bar,tmp_f,SIG,EPS,m,delta)*( (x_i - x_im1)/(func_prime(x_i,f,f_bar,tmp_f,SIG,EPS,m,delta) - func_prime(x_im1,f,f_bar,tmp_f,SIG,EPS,m,delta)));
	
	while ( (k < max_iters) && (fabs(x_ip1 - x_i) > epsilon) ) {
		x_im1 = x_i;
		x_i = x_ip1;
		
		x_ip1 = x_i - func_prime(x_i,f,f_bar,tmp_f,SIG,EPS,m,delta)*( (x_i - x_im1)/(func_prime(x_i,f,f_bar,tmp_f,SIG,EPS,m,delta) - func_prime(x_im1,f,f_bar,tmp_f,SIG,EPS,m,delta)));
		k++;
	}
	
	return x_ip1;
}
#endif
