/* Illustrates the linking of a C program with the shared object libPHCpack.
 * Assuming the code is saved in the file with name 'welcome.c'
 * the executable 'welcome' is made executing the statement
 * gcc -o welcome welcome.c ../lib/libPHCpack.so
 * where the shared object libPHCpack is present in the folder ../lib. */

#include <stdlib.h>
#include <stdio.h>

extern int _ada_use_c2phc ( int job, int *a, int *b, double *c, int v );

int main ( int argc, char* argv[] )
{
   int *a;
   int *b;
   double *c;

 // writes the welcome banner to PHC

   int fail = _ada_use_c2phc(0, a, b, c, 1);

   int len;
   int name[30];

 // retrieves the version string of PHCpack

   fail = _ada_use_c2phc(999, &len, name, c, 1);

   char *version = calloc(30, sizeof(char));

   for(int i=0; i<30; i++) version[i] = (char) name[i];

   printf("%s\n", version);

   return 0;
}
