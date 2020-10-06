// The file parse_run_arguments.cpp contains the defintion of the function
// with prototype in parse_run_arguments.h.

#include <iostream>
#include <cstdlib>
#include "parse_run_arguments.h"

using namespace std;

int parse_run_arguments
 ( int argc, char *argv[], int *blocksize, int *dim, int *freq, int *mode )
{
   if(argc < 5)
   {
      cout << argv[0] << " needs four parameters, for example:" << endl;
      cout << argv[0] << " blocksize dim freq mode" << endl;
      cout << " where blocksize is the number of threads in a block," << endl;
      cout << "       dim is the dimension of the problem," << endl;
      cout << "       freq is the number of repeated runs, and" << endl;
      cout << "       mode is 0, 1, or 2 for execution mode." << endl;
      cout << "Please try again ..." << endl;

      return 1;
   }
   else
   {
      *blocksize = atoi(argv[1]);  // number of threads in a block
      *dim = atoi(argv[2]);        // dimension
      *freq = atoi(argv[3]);       // number of repeated runs
      *mode = atoi(argv[4]);       // execution mode

      return 0;
   }
}
