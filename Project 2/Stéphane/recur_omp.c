#include "walltime.h"
#include <math.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char *argv[]) {
	int N = 2000000000;
	double up = 1.00000001;
	double Sn = 1.00000001;
	int n;
	/* allocate memory for the recursion */
	double *opt = (double *)malloc((N + 1) * sizeof(double));
	
	if (opt == NULL)
	die("failed to allocate problem size");
	
	double time_start = wall_time();
	// TODO: YOU NEED TO PARALLELIZE THIS LOOP
	
	#pragma omp parallel firstprivate(Sn)
	{
		int tid = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
		int chunksize = N/nthreads;
		double offset = pow(up, tid*chunksize);
		Sn = Sn*offset;

		for (int i = 0; i < chunksize; i++) {
			opt[i + tid*chunksize] = Sn;
			Sn *= up;
		}
	}
	
	Sn = opt[N-1]*up;
	
	
	//~ #pragma omp parallel firstprivate(Sn)
	//~ {
		//~ int tid = omp_get_thread_num();
		//~ int nthreads = omp_get_num_threads();
		//~ int chunksize = N/nthreads;		
		//~ #pragma omp for schedule (dynamic [chunksize])
			  //~ for (n = 0; n <= N; ++n) {
				  //~ opt[n] = Sn;
				  //~ Sn *= up;
			  //~ }
	//~ }
			
	printf("Parallel RunTime   :  %f seconds\n", wall_time() - time_start);
	printf("Final Result Sn    :  %.17g \n", Sn);
	
	double temp = 0.0;
	for (n = 0; n <= N; ++n) {
		temp += opt[n] * opt[n];
	}
	
	printf("Result ||opt||^2_2 :  %f\n", temp / (double)N);
	printf("\n");
	
	return 0;
}
