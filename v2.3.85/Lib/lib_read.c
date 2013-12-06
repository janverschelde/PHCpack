/* prompts the user for a system and writes it to screen */

#include <stdio.h>
#include <stdlib.h>
#include "phcpack.h"

#define v 1  /* verbose flag:
                 0 no output during the computations,
                 1 only one-line message after call to phcpack */

int main ( int argc, char *argv[] )
{
   int fail,n;
   char outfile[80],ch;

   adainit();

   printf("\nGive name of the output file : ");
   scanf("%s",outfile);
   n = (int) strlen(outfile);
   fail = define_output_file_with_string(n,outfile);
   fail = write_string_to_defined_output_file(n,outfile);
   fail = close_output_file();
   scanf("%c",&ch); /* skip newline symbol */

   printf("\nReading and writing a polynomial system ...\n");

   fail = read_target_system();
   if(v==1) printf("-> read_target_system returns %d\n",fail);
   fail = copy_target_system_to_container();
   if(v==1) printf("-> copy_target_system_to_container returns %d\n",fail);
   printf("\nThe system read :\n");
   fail = print_system();
   if(v==1) printf("-> print_system returns %d\n",fail);

   adafinal();

   return 0;
}
