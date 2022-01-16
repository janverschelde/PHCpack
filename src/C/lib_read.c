/* prompts the user for a system and writes it to screen */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phcpack.h"
#include "solcon.h"

#define v 1  /* verbose flag:
                 0 no output during the computations,
                 1 only one-line message after call to phcpack */

int input_argument_free ( void );
/*
 * In this test program, the Ada code prompts the user for the input 
 * file names and the calls to the Ada code is thus argument free,
 * as far as the input data concerns. */

int named_input_files ( void );
/*
 * The names of the input files in this test are read by the C function
 * and passed to the Ada routine who will then use the name to read the
 * input data from file and into the containers. */

int main ( int argc, char *argv[] )
{
   adainit();

   // input_argument_free();
   named_input_files();

   adafinal();

   return 0;
}

int input_argument_free ( void )
{
   int fail,n;
   char outfile[80],ch;
   char somesys[22] = "\nTITLE : some system\n";

   printf("\nReading and writing a polynomial system ...\n");

   fail = read_standard_target_system();
   if(v==1) printf("-> read_target_system returns %d\n",fail);
   fail = copy_target_system_to_container();
   if(v==1) printf("-> copy_target_system_to_container returns %d\n",fail);
   printf("\nThe system read :\n");
   fail = print_system();
   if(v==1) printf("-> print_system returns %d\n",fail);

   printf("\nGive name of the output file : ");
   scanf("%s",outfile);
   n = (int) strlen(outfile);
   fail = define_output_file_with_string(n,outfile);
   if(v==1) printf("-> define_output_file_with_string returns %d\n",fail);
   fail = print_system();
   if(v==1) printf("-> print_system returns %d\n",fail);
 //  fail = write_string_to_defined_output_file(21,"\nTITLE : some system\n");
   fail = write_string_to_defined_output_file(21,somesys);
   if(v==1) printf("-> write_string_to_defined_output_file returns %d\n",fail);
   fail = close_output_file();
   if(v==1) printf("-> close_output_file returns %d\n",fail);

   return 0;
}

int named_input_files ( void )
{
   int n, fail, len;
   char inputfile[80];
   int mode = 7;
 /* reading target system :  0 for standard, 1 for double double,
                             2 for quad double, 3 for multiprecision
    reading start system :  4 for standard, 5 for double double,
                            6 for quad double, 7 for multiprecision */

   if(mode < 4) printf("\nReading the target system from file ...\n");
   if(mode > 3) printf("\nReading the start system from file ...\n");
   printf("Give the file name : "); scanf("%s", inputfile);
   n = (int) strlen(inputfile);

   if(mode == 0)
   {
      printf("\n... test reading in standard doubles ...\n");
      fail = read_standard_target_system_from_file(n, inputfile);
      fail = copy_container_to_target_system();
      printf("\nThe system read :\n"); fail = write_standard_target_system();
   }
   else if(mode == 1)
   {
      printf("\n... test reading in double doubles ...\n");
      fail = read_dobldobl_target_system_from_file(n, inputfile);
      fail = copy_dobldobl_container_to_target_system();
      printf("\nThe system read :\n"); fail = write_dobldobl_target_system();
   }
   else if(mode == 2)
   {
      printf("\n... test reading in quad doubles ...\n");
      fail = read_quaddobl_target_system_from_file(n, inputfile);
      fail = copy_quaddobl_container_to_target_system();
      printf("\nThe system read :\n"); fail = write_quaddobl_target_system();
   }
   else if(mode == 3)
   {
      printf("\n... test reading in multiprecision ...\n");
      fail = read_multprec_target_system_from_file(80, n, inputfile);
      fail = copy_multprec_container_to_target_system();
      printf("\nThe system read :\n"); fail = write_multprec_target_system();
   }
   else if(mode == 4)
   {
      printf("\n... test reading in standard doubles ...\n");
      fail = read_standard_start_system_from_file(n, inputfile);
      fail = copy_container_to_start_system();
      printf("\nThe system read :\n"); fail = write_standard_start_system();
      fail = solcon_number_of_standard_solutions(&len);
      printf("\nNumber of start solutions in container : %d\n", len);
   }
   else if(mode == 5)
   {
      printf("\n... test reading in double doubles ...\n");
      fail = read_dobldobl_start_system_from_file(n, inputfile);
      fail = copy_dobldobl_container_to_start_system();
      printf("\nThe system read :\n"); fail = write_dobldobl_start_system();
      fail = solcon_number_of_dobldobl_solutions(&len);
      printf("\nNumber of start solutions in container : %d\n", len);
   }
   else if(mode == 6)
   {
      printf("\n... test reading in quad doubles ...\n");
      fail = read_quaddobl_start_system_from_file(n, inputfile);
      fail = copy_quaddobl_container_to_start_system();
      printf("\nThe system read :\n"); fail = write_quaddobl_start_system();
      fail = solcon_number_of_quaddobl_solutions(&len);
      printf("\nNumber of start solutions in container : %d\n", len);
   }
   else
   {
      printf("\n... test reading in multiprecision ...\n");
      fail = read_multprec_start_system_from_file(80, n, inputfile);
      fail = copy_multprec_container_to_start_system();
      printf("\nThe system read :\n"); fail = write_multprec_start_system();
      fail = solcon_number_of_multprec_solutions(&len);
      printf("\nNumber of start solutions in container : %d\n", len);
   }
   return 0;
}
