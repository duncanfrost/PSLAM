#include "Runnable.h"

void amoRunnable::start()
{
   ret=pthread_create(&thid,NULL,executer,(void*)this);
}
void amoRunnable::join()
{
   pthread_join(thid,NULL);
}
void* executer(void* param)
{
   ObjectToRun *obj=(ObjectToRun*)param;
   obj->run();
}

/*void amoRunnable::kill()
{
   pthread_cancel(thid);
   join();
}*/
