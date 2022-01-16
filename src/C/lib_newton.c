/* reads a target and start system and then solves the target system,
   using the start system in an artificial-parameter homotopy */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phcpack.h"
#include "syscon.h"
#include "solcon.h"

void test_validate ( void );
/* tests the validation with Newton */

void test_standard_Newton_step ( void );
/* runs the Newton step with standard arithmetic */

void test_dobldobl_Newton_step ( void );
/* runs the Newton step with double double arithmetic */

void test_quaddobl_Newton_step ( void );
/* runs the Newton step with quad double arithmetic */

void test_multprec_Newton_step ( void );
/* runs the Newton step with multiprecision arithmetic */

void test_standard_Newton_Laurent_step ( void );
/* runs the Newton step with standard arithmetic
 * on Laurent polynomial systems */

void test_dobldobl_Newton_Laurent_step ( void );
/* runs the Newton step with double double arithmetic
 * on Laurent polynomial systems */

void test_quaddobl_Newton_Laurent_step ( void );
/* runs the Newton step with quad double arithmetic
 * on Laurent polynomial systems */

void test_multprec_Newton_Laurent_step ( void );
/* runs the Newton step with multiprecision arithmetic
 * on Laurent polynomial systems */

void test_deflate ( void );
/* calls the standard double precision deflate */

void test_varbprec_Newton_Laurent_step ( void );
/* test on variable precision Newton step */

int main ( int argc, char *argv[] )
{
   int choice;
   char ch;

   adainit();

   printf("\nMENU to run Newton's method : \n");
   printf("  0. test plain validation on a start system;\n");
   printf("  1. run Newton step with standard double arithmetic;\n");
   printf("  2. run Newton step with double double arithmetic;\n");
   printf("  3. run Newton step with quad double arithmetic;\n");
   printf("  4. run Newton step with multiprecision arithmetic;\n");
   printf("  5. standard double precision deflation with defaults;\n");
   printf("  6. standard double Newton step on Laurent system;\n");
   printf("  7.   double double Newton step on Laurent system;\n");
   printf("  8.     quad double Newton step on Laurent system;\n");
   printf("  9.  multiprecision Newton step on Laurent system;\n");
   printf(" 10. variable precision Newton step on Laurent system.\n");
   printf("Type 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, or 10 to select : ");
   scanf("%d",&choice);
   scanf("%c",&ch); /* skip newline symbol */

   if(choice == 0)
      test_validate();
   else if(choice == 1)
      test_standard_Newton_step();
   else if(choice == 2)
      test_dobldobl_Newton_step();
   else if(choice == 3)
      test_quaddobl_Newton_step();
   else if(choice == 4)
      test_multprec_Newton_step();
   else if(choice == 5)
      test_deflate();
   else if(choice == 6)
      test_standard_Newton_Laurent_step();
   else if(choice == 7)
      test_dobldobl_Newton_Laurent_step();
   else if(choice == 8)
      test_quaddobl_Newton_Laurent_step();
   else if(choice == 9)
      test_multprec_Newton_Laurent_step();
   else if(choice == 10)
      test_varbprec_Newton_Laurent_step();
   else
      printf("invalid selection, please try again\n");

   adafinal();

   return 0;
}

void test_validate ( void )
{
   int fail;

   printf("\nCalling Newton validation in PHCpack...\n");
   fail = read_standard_start_system();
   fail = define_output_file();
   printf("\nSee the output file results ...\n\n");
   fail = copy_start_system_to_container();
   fail = copy_start_solutions_to_container();
   fail = validate_solutions();
}

void test_standard_Newton_step ( void )
{
   int fail,dim,len;

   printf("\nRunning Newton step with standard arithmetic ...\n");
   fail = syscon_read_standard_system();
   fail = syscon_number_of_standard_polynomials(&dim);
   printf("The system container has %d polynomials.\n",dim);
   fail = solcon_read_standard_solutions();
   fail = solcon_number_of_standard_solutions(&len);
   printf("The solution container has size %d.\n",len);
   fail = solcon_dimension_of_standard_solutions(&dim);
   printf("The solutions in the container have dimension %d.\n",dim);
   fail = standard_Newton_step();
   printf("The solutions after the Newton step :\n");
   fail = solcon_write_standard_solutions();
}

void test_standard_Newton_Laurent_step ( void )
{
   int fail,dim,len;

   printf("\nRunning Newton step with standard arithmetic ...\n");
   fail = syscon_read_standard_Laurent_system();
   fail = syscon_number_of_standard_Laurentials(&dim);
   printf("The system container has %d Laurent polynomials.\n",dim);
   fail = solcon_read_standard_solutions();
   fail = solcon_number_of_standard_solutions(&len);
   printf("The solution container has size %d.\n",len);
   fail = solcon_dimension_of_standard_solutions(&dim);
   printf("The solutions in the container have dimension %d.\n",dim);
   fail = standard_Newton_Laurent_step();
   printf("The solutions after the Newton step :\n");
   fail = solcon_write_standard_solutions();
}

void test_dobldobl_Newton_step ( void )
{
   int fail,dim,len;

   printf("\nRunning Newton step with double double arithmetic ...\n");
   fail = syscon_read_dobldobl_system();
   fail = syscon_number_of_dobldobl_polynomials(&dim);
   printf("The system container has %d polynomials.\n",dim);
   fail = solcon_read_dobldobl_solutions();
   fail = solcon_number_of_dobldobl_solutions(&len);
   printf("The solution container has size %d.\n",len);
   fail = solcon_dimension_of_dobldobl_solutions(&dim);
   printf("The solutions in the container have dimension %d.\n",dim);
   fail = dobldobl_Newton_step();
   printf("The solutions after the Newton step :\n");
   fail = solcon_write_dobldobl_solutions();
}

void test_dobldobl_Newton_Laurent_step ( void )
{
   int fail,dim,len;

   printf("\nRunning Newton step with double double arithmetic ...\n");
   fail = syscon_read_dobldobl_Laurent_system();
   fail = syscon_number_of_dobldobl_Laurentials(&dim);
   printf("The system container has %d Laurent polynomials.\n",dim);
   fail = solcon_read_dobldobl_solutions();
   fail = solcon_number_of_dobldobl_solutions(&len);
   printf("The solution container has size %d.\n",len);
   fail = solcon_dimension_of_dobldobl_solutions(&dim);
   printf("The solutions in the container have dimension %d.\n",dim);
   fail = dobldobl_Newton_Laurent_step();
   printf("The solutions after the Newton step :\n");
   fail = solcon_write_dobldobl_solutions();
}

void test_quaddobl_Newton_step ( void )
{
   int fail,dim,len;

   printf("\nRunning Newton step with quad double arithmetic ...\n");
   fail = syscon_read_quaddobl_system();
   fail = syscon_number_of_quaddobl_polynomials(&dim);
   printf("The system container has %d polynomials.\n",dim);
   fail = solcon_read_quaddobl_solutions();
   fail = solcon_number_of_quaddobl_solutions(&len);
   printf("The solution container has size %d.\n",len);
   fail = solcon_dimension_of_quaddobl_solutions(&dim);
   printf("The solutions in the container have dimension %d.\n",dim);
   fail = quaddobl_Newton_step();
   printf("The solutions after the Newton step :\n");
   fail = solcon_write_quaddobl_solutions();
}

void test_quaddobl_Newton_Laurent_step ( void )
{
   int fail,dim,len;

   printf("\nRunning Newton step with quad double arithmetic ...\n");
   fail = syscon_read_quaddobl_Laurent_system();
   fail = syscon_number_of_quaddobl_Laurentials(&dim);
   printf("The system container has %d Laurent polynomials.\n",dim);
   fail = solcon_read_quaddobl_solutions();
   fail = solcon_number_of_quaddobl_solutions(&len);
   printf("The solution container has size %d.\n",len);
   fail = solcon_dimension_of_quaddobl_solutions(&dim);
   printf("The solutions in the container have dimension %d.\n",dim);
   fail = quaddobl_Newton_Laurent_step();
   printf("The solutions after the Newton step :\n");
   fail = solcon_write_quaddobl_solutions();
}

void test_multprec_Newton_step ( void )
{
   int fail,dim,len,deci;

   printf("\nRunning Newton step with multiprecision arithmetic ...\n");
   printf("\ngive the number of decimal places in the working precision : ");
   scanf("%d",&deci);
   fail = syscon_read_multprec_system(deci);
   fail = syscon_number_of_multprec_polynomials(&dim);
   printf("The system container has %d polynomials.\n",dim);
   fail = solcon_read_multprec_solutions();
   fail = solcon_number_of_multprec_solutions(&len);
   printf("The solution container has size %d.\n",len);
   fail = solcon_dimension_of_multprec_solutions(&dim);
   printf("The solutions in the container have dimension %d.\n",dim);
   fail = multprec_Newton_step(deci);
   printf("The solutions after the Newton step :\n");
   fail = solcon_write_multprec_solutions();
}

void test_multprec_Newton_Laurent_step ( void )
{
   int fail,dim,len,deci;

   printf("\nRunning Newton step with multiprecision arithmetic ...\n");
   printf("\ngive the number of decimal places in the working precision : ");
   scanf("%d",&deci);
   fail = syscon_read_multprec_Laurent_system(deci);
   fail = syscon_number_of_multprec_Laurentials(&dim);
   printf("The system container has %d Laurent polynomials.\n",dim);
   fail = solcon_read_multprec_solutions();
   fail = solcon_number_of_multprec_solutions(&len);
   printf("The solution container has size %d.\n",len);
   fail = solcon_dimension_of_multprec_solutions(&dim);
   printf("The solutions in the container have dimension %d.\n",dim);
   fail = multprec_Newton_Laurent_step(deci);
   printf("The solutions after the Newton step :\n");
   fail = solcon_write_multprec_solutions();
}

void test_deflate ( void )
{
   int fail,dim,len;
   const int maxitr = 3;
   const int maxdef = 3;
   const double tolerr = 1.0e-8;
   const double tolres = 1.0e-8;
   const double tolrnk = 1.0e-6;

   printf("\nRunning deflation ...\n");
   fail = read_standard_start_system();
   fail = copy_start_system_to_container();
   fail = copy_start_solutions_to_container();
   fail = syscon_number_of_standard_polynomials(&dim);
   printf("The system container has %d polynomials.\n",dim);
   fail = solcon_number_of_standard_solutions(&len);
   printf("The solution container has size %d.\n",len);
   fail = standard_deflate(maxitr,maxdef,tolerr,tolres,tolrnk);
   printf("The solutions after deflation :\n");
   fail = solcon_write_standard_solutions();
}

void test_varbprec_Newton_Laurent_step ( void )
{
   int fail,dim,lenpols,nf,nq,nv,ns,lensols;
   int wanted,maxitr,maxprc;
   char infile[80],c;
   char *system;

   printf("\nRunning variable precision Newton step ...\n");
   printf("\nGive name of the input file : ");
   scanf("%s",infile);
   nf = (int) strlen(infile);
   system = read_polynomials_from_file(nf,infile,&lenpols,&nq,&nv,&fail);
   printf("The system has %d characters :\n%s\n",lenpols,system);
   scanf("%c",&c); /* skip new line symbol */
   printf("\n");
   fail = solcon_read_multprec_solutions();
   fail = solcon_number_of_multprec_solutions(&lensols);
   printf("The solution container has size %d.\n",lensols);
   fail = solcon_dimension_of_multprec_solutions(&dim);
   printf("The solutions in the container have dimension %d.\n",dim);
   if(dim != nv)
      printf("Dimension of solutions does not match number of variables!\n");
   else
   {
      printf("Calling the variable precision Newton steps method ...\n");
      printf("-> give the wanted number of decimal places : ");
      scanf("%d",&wanted);
      printf("-> give the maximum number of Newton steps : ");
      scanf("%d",&maxitr);
      printf("-> give the maximum number of places in the precision : ");
      scanf("%d",&maxprc);
      fail = varbprec_Newton_Laurent_step
        (dim,wanted,maxitr,maxprc,lenpols,system);
      printf("The solutions after the variable precision Newton steps :\n");
      fail = solcon_write_multprec_solutions();
   }
}
