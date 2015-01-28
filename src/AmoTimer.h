#pragma once


#include <sys/time.h>

#define USE_OMP_C
#ifndef USE_OMP_C
class amoTimer
{
  public:
	amoTimer(){};
	void start()
	{			    
			gettimeofday(&t1, NULL);		
	};
	void stop(std::string addPrint="")
	{
			gettimeofday(&t2, NULL);  
			double elapsedTime; 
			elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      
			elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms		
			std::cout << "time "<<addPrint<<": "<<elapsedTime << " ms.\n";
	};
  private:
	timeval t1, t2;     

};
#else

#include <omp.h>
class amoTimer
{
  public:
	amoTimer(){};
	void start()
	{			    
			t1 = omp_get_wtime( );		
	};
	void stop(std::string addPrint="")
	{
			t2 = omp_get_wtime( ); 
			std::cout << "time "<<addPrint<<": "<<(t2-t1)*1000. << " ms.\n";
	};
	float getTimeDiffms(std::string addPrint="")
	{
			t2 = omp_get_wtime( ); 
			return (t2-t1)*1000.;
	};
  private:
	double t1, t2;     

};
#endif



