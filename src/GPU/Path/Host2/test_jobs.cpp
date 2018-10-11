// Test on the operations in the Job Queue

#include <unistd.h>
#include "jobqueue.h"
#include "worker.h"

using namespace std;

typedef struct // type for the do_job function
{
   int label;  // for the idn of the worker
   void *jobs; // for the job queue
}
WorkItem;      // argument for the do_job function

int process_jobs ( int nbt, JobQueue* jobs );
/*
 * Processing the jobs with nbt threads. */

void* do_job ( void* args );
/*
 * Processes the jobs in the job queue, given in args. */

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
   Worker *crew = (Worker*)calloc(nbt,sizeof(Worker));
   WorkItem *wit = (WorkItem*) calloc(nbt,sizeof(WorkItem));

   cout << "Assigning unique identification to each worker ..." << endl;
   for(int i=0; i<nbt; i++) crew[i].idn = i;
   cout << "Writing the data at each worker ..." << endl;
   for(int i=0; i<nbt; i++) crew[i].write();
   cout << "Starting the work crew ..." << endl;
   for(int i=0; i<nbt; i++)
   {
      wit[i].label = crew[i].idn;
      wit[i].jobs = jobs;

      cout << "Launching worker " << wit[i].label << endl;

      crew[i].work(&do_job,(void*)&wit[i]);
   }
   cout << "Waiting for the workers to finish ..." << endl;
   for(int i=0; i<nbt; i++) crew[i].join();

   return 0;
}

void* do_job ( void* args )
{
   WorkItem *wit = (WorkItem*) args;
   JobQueue *jobs = (JobQueue*) wit->jobs;

   cout << "Worker " << wit->label << " entered job loop ..." << endl;

   while(1)
   {
      int nextjob = jobs->get_next_job();

      if(nextjob == -1) break;

      cout << "Worker " << wit->label
           << " has job : " << nextjob
           << " with work " << jobs->work[nextjob] << endl;

      sleep(jobs->work[nextjob]);
   }

   return NULL;
}
