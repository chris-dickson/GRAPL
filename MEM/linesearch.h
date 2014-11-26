/*--------------------------------------------------------------------------------------*/
/* linesearch.h																			*/
/*	Header for various linesearch functions												*/
/*																						*/
/*--------------------------------------------------------------------------------------*/

#ifndef LINESEARCH_H_
#define LINESEARCH_H_

#define MAX_BISECT_ITERS 10000

double solve_BISECT_THETA(double A_in, double B_in, double epsilon, double SIG_in, int m_in, double *f_in, int max_iters);
double solve_BISECT_THETA_ECMM(double A_in, double B_in, double epsilon, double SIG_in, int m_in, double *f_in, int max_iters);


double pot_fun(double TAU , double *f, double *f_bar, double *tmp_f, double SIG, double EPS, double m);

double pot_fun_deriv(double TAU , double *f, double *f_bar, double *tmp_f, double SIG, double EPS, double m, double delta);

double solve_BISECT_potfun(double A_in, double B_in, double epsilon, double delta, double SIG, int m, double *f_bar, double *f, double *tmp_f, double EPS, int max_iters);

double secant(double x_0, double x_1, double epsilon, double delta, int max_iters, 
			  double func(double TAU , double *f, double *f_bar, double *tmp_f, double SIG, double EPS, double m),
			  double *f, double *f_bar, double *tmp_f, double SIG, double EPS, double m);


#endif
