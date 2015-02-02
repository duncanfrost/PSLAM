#include "AmoDefines.h"

void InitProcAndGPU()
{
#ifdef USE_OMP
	if(omp_get_max_threads()>4)
		omp_set_num_threads(4);
	else if(omp_get_max_threads()>2)
		omp_set_num_threads(2);
#endif
		
}

