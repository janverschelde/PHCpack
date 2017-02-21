// Test on the operations in the Job Queue

#include "jobqueue.h"

using namespace std;

int process_jobs ( int nbt, JobQueue* jobs );
/*
 * Processing the jobs with nbt threads. */

int main ( void )
{
   int nbr;
   cout << "Give the number of jobs : ";
   cin >> nbr;

   JobQueue jobs(nbr);

   jobs.write();

   while(1)
   {
      int nextjob = jobs.get_next_job();

      if(nextjob == -1) break;

      cout << "Index of the next job : " << nextjob
           << " with work " << jobs.work[nextjob] << endl;
   }

   *jobs.nextjob = -1; // reset the job queue
   int nbt;
   cout << "Give the number threads : ";
   cin >> nbt;

   int fail = process_jobs(nbt,&jobs);

   return 0;
}

int process_jobs ( int nbt, JobQueue* jobs )
{
   return 0;
}
