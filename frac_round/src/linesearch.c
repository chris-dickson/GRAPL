double solve_BISECT(double A_in, double B_in, double epsilon, double SIG_in, int m_in, double *f_in) {
	FILE *fp;
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
	
	fp = fopen("SolveOut.txt","w");
	
	if (fp == NULL) {
		printf("Can't open filve SolveOut.txt\n"); fflush(stdout);
		exit(1);
	}
	
	a = A_in;
	b = B_in;
	
	f_a = objective(a,SIG_in,m_in,f_in);
	f_b = objective(b,SIG_in,m_in,f_in);

	i = 0;
	fprintf(fp,"Iteration %6d : a = %e , b = %e , f_a = %e , f_b = %e\n",i,a,b,f_a,f_b);
		
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
	
	while ( (b-a)>epsilon ) {

		ret = objective((a+b)/2.0,SIG_in,m_in,f_in);
		
	//	printf("%f\n",ret);
		
		if ((ret*f_b) < 0) {
			a = (a+b)/2.0;
			f_a = objective(a,SIG_in,m_in,f_in);
		}
		
		if ((f_a*ret) < 0) {
			b = (a+b)/2.0;
			f_b = objective(b,SIG_in,m_in,f_in);
		}
		i++;
		fprintf(fp,"Iteration %6d : a = %e , b = %e , f_a = %e , f_b = %e\n",i,a,b,f_a,f_b);	
	}
	fclose(fp);
	return((a+b)/2.0);
}				

