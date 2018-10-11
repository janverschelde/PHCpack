// Test on the methods in the Worker class.

#include <unistd.h>
#include "worker.h"

using namespace std;

void* write ( void* args );
/*
 * Expects a string in the argument args.
 * The string is written to screen. */

int main ( void )
{
   int nbr;
   cout << "Give an identification number for a worker : ";
   cin >> nbr;

   Worker thread(nbr);

   thread.write();

   string hello = "Hello world!\0";

   int fail;
   fail = thread.work(&write,(void*)&hello);

   cout << "waiting for the thread to join ..." << endl;
   thread.join();

   cout << "giving the worker more work ...";

   fail = thread.work(&write,(void*)&hello);

   cout << "waiting for the thread to join ..." << endl;
   thread.join();

   cout << "done!" << endl;

   return 0;
}

void* write ( void* args )
{
   string* message = (string*) args;

   cout << *message << endl;

   cout << "sleeping for 2 time units ..." << endl;

   sleep(2);

   return NULL;
}
