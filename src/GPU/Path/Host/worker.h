/* The class Worker encapsulates a pthread for shared memory multithreaded
   execution of any given number of functions. */

#ifndef _WORKER_H_
#define _WORKER_H_

#include <cstdlib>
#include <iostream>
#include <pthread.h>

class Worker
{
   private:

      pthread_t thread;      // the thread of execution
      pthread_attr_t attrib; // attributes of the thread

   public:

      int idn;           // identification number

      Worker ( int nbr ) // initialize the identification number
      {
         this->idn = nbr;
      }

      void write ( void ) // writes the identification number
      {
         std::cout << "Hello from thread " << this->idn
                   << "." << std::endl;
      }

      int work ( void* do_job ( void* args ), void* args, int verbose = 1 )
      // launches the thread
      {
         if(verbose > 0)
             std::cout << "worker " << this->idn
                       << " starts to work ..." << std::endl;

         pthread_attr_init(&attrib);
         pthread_create(&(this->thread),&attrib,do_job,args);

         if(verbose > 0)
             std:: cout << "worker " << this->idn
                        << " returns from thread creation" << std::endl;

         return 0;
      }

      int join ( void ) // wait for the thread to finish
      {
         pthread_join(this->thread,NULL);

         return 0;
      }

      ~Worker ( void ) // destroys the thread
      {
         pthread_exit(NULL);
      }
};

#endif /* _WORKER_H_ */
