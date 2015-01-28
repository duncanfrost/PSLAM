//definition of runnable object
#ifndef __RUNNABLE_H___
#define __RUNNABLE_H___


#include<pthread.h>
#include<sys/types.h>

using namespace std;
class ObjectToRun
{
    public:
      virtual void run(){}
};
void* executer(void* param);
class amoRunnable : public ObjectToRun
{
  pthread_t thid;
  int ret;
  public:
   void start();
   void join();
   void kill();
};

#endif