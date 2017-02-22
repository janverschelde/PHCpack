/* jobqueue.h for managing a queue of jobs */

#ifndef _JOBQUEUE_H_
#define _JOBQUEUE_H_

#include <cstdlib>
#include <iostream>
#include <pthread.h>

class JobQueue
{
   private:

      pthread_mutex_t read_lock;

   public:

      int nbr;       // number of jobs
      int *nextjob;  // index of the next job
      int *work;     // array of nbr jobs

      JobQueue ( int n )
      {
         this->nbr = n;
         this->nextjob = (int*) calloc(1,sizeof(int));
         *(this->nextjob) = -1;
         this->work = (int*) calloc(n,sizeof(int));
         for(int i=0; i<n; i++)
            this->work[i] = 1 + rand() % 5;
         pthread_mutex_init(&read_lock, NULL);
      }

      void write ( void )
      {
         std::cout << "Number of jobs : " << this->nbr << std::endl;
         for(int i=0; i<this->nbr; i++)
            std::cout << "job " << i <<  " : "
                      << this->work[i] << std::endl;
         std::cout << "Current job : " << *(this->nextjob) << std::endl;
      }

      int get_next_job ( int verbose = 1 )
      /*
       * If there is no next job, then returns -1,
       * else returns the index of the next job. */
      {
         int thenextjob = -1;

         pthread_mutex_lock(&read_lock);
         int* ptr = this->nextjob;
         if(verbose > 0) std::cout << "  at job " << *ptr << std::endl; 
         if(*ptr < this->nbr-1) thenextjob = ++(*ptr);
         if(verbose > 0) std::cout << "next job " << thenextjob << std::endl; 
         pthread_mutex_unlock(&read_lock);
 
         return thenextjob;     
      }
};

#endif /* _JOBQUEUE_H_ */
