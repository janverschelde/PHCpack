/* Tests the membership test for a witness set and test point. */

#include <stdio.h>
#include "syscon.h"
#include "solcon.h"
#include "witset.h"

#define verbose 1 /* verbose flag */

void standard_read_point ( int n, double *x );
/*
 * DESCRIPTION :
 *   Prompts the user for 2*n doubles, for the real and imaginary parts
 *   of the standard double precision coordinates of the point x. */

void dobldobl_read_point ( int n, double *x );
/*
 * DESCRIPTION :
 *   Prompts the user for 4*n doubles, for the real and imaginary parts
 *   of the double double precision coordinates of the point x. */

void quaddobl_read_point ( int n, double *x );
/*
 * DESCRIPTION :
 *   Prompts the user for 8*n doubles, for the real and imaginary parts
 *   of the quad double precision coordinates of the point x. */

int standard_membership_test ( void );
/*
 * DESCRIPTION :
 *   Prompts for a witness set defined by an ordinary polynomial system,
 *   and prompts for a test point.
 *   Runs the membership test in standard double precision. */

int dobldobl_membership_test ( void );
/*
 * DESCRIPTION :
 *   Prompts for a witness set defined by an ordinary polynomial system
 *   and prompts for a test point.
 *   Runs the membership test in double double precision. */

int quaddobl_membership_test ( void );
/*
 * DESCRIPTION :
 *   Prompts for a witness set defined by an ordinary polynomial system
 *   and prompts for a test point.
 *   Runs the membership test in quad double precision. */

int standard_Laurent_membership_test ( void );
/*
 * DESCRIPTION :
 *   Prompts for a witness set defined by a Laurent polynomial system,
 *   and prompts for a test point.
 *   Runs the membership test in standard double precision. */

int dobldobl_Laurent_membership_test ( void );
/*
 * DESCRIPTION :
 *   Prompts for a witness set defined by a Laurent polynomial system
 *   and prompts for a test point.
 *   Runs the membership test in double double precision. */

int quaddobl_Laurent_membership_test ( void );
/*
 * DESCRIPTION :
 *   Prompts for a witness set defined by a Laurent polynomial system
 *   and prompts for a test point.
 *   Runs the membership test in quad double precision. */

int main ( int argc, char *argv[] )
{
   char ans;
   int precision,fail;

   adainit();

   printf("\nMENU for the working precision :\n");
   printf("  0. standard double precision;\n");
   printf("  1. double double precision;\n");
   printf("  2. quad double precision.\n");
   printf("Type 0, 1, or 2 to make a choice : ");
   scanf("%d",&precision);
   scanf("%c",&ans); /* skip end of line character */

   printf("Laurent polynomial system ? (y/n) ");
   scanf("%c",&ans);

   if(ans == 'y')
   {
      scanf("%c",&ans); 
      if(precision == 0)
         fail = standard_Laurent_membership_test();
      else if(precision == 1)
         fail = dobldobl_Laurent_membership_test();
      else if(precision == 2)
         fail = quaddobl_Laurent_membership_test();
      else
         printf("Selected precision level is not supported.\n");
   }
   else
   {
      scanf("%c",&ans);
      if(precision == 0)
         fail = standard_membership_test();
      else if(precision == 1)
         fail = dobldobl_membership_test();
      else if(precision == 2)
         fail = quaddobl_membership_test();
      else
         printf("Selected precision level is not supported.\n");
   }
   adafinal();

   return 0;
}

void standard_read_point ( int n, double *x )
{
   int k;

   for(k=0; k<n; k++)
   {
      printf("Give the real part for x[%d] : ",k);
      scanf("%lf",&x[2*k]);
      printf("Give the imaginary part for x[%d] : ",k);
      scanf("%lf",&x[2*k+1]);
   }
   if(verbose>0)
   {
      printf("The coordinates of the point :\n");
      for(k=0; k<n; k++)
         printf("x[%d] = %.15e  %.15e\n",k,x[2*k],x[2*k+1]);
   }
}

void dobldobl_read_point ( int n, double *x )
{
   int k;

   for(k=0; k<n; k++)
   {
      printf("Give the real part for x[%d] : ",k);
      scanf("%lf",&x[4*k]);
      scanf("%lf",&x[4*k+1]);
      printf("Give the imaginary part for x[%d] : ",k);
      scanf("%lf",&x[4*k+2]);
      scanf("%lf",&x[4*k+3]);
   }
   if(verbose>0)
   {
      printf("The coordinates of the point :\n");
      for(k=0; k<n; k++)
         printf("x[%d] = %.15e  %.15e  %.15e  %.15e\n",
                    k,x[4*k],x[4*k+1],x[4*k+2],x[4*k+3]);
   }
}

void quaddobl_read_point ( int n, double *x )
{
   int k;

   for(k=0; k<n; k++)
   {
      printf("Give the real part for x[%d] : ",k);
      scanf("%lf",&x[8*k]); scanf("%lf",&x[8*k+1]);
      scanf("%lf",&x[8*k+2]); scanf("%lf",&x[8*k+3]);
      printf("Give the imaginary part for x[%d] : ",k);
      scanf("%lf",&x[8*k+4]); scanf("%lf",&x[8*k+5]);
      scanf("%lf",&x[8*k+6]); scanf("%lf",&x[8*k+7]);
   }
   if(verbose>0)
   {
      printf("The coordinates of the point :\n");
      for(k=0; k<n; k++)
         printf(
            "x[%d] = %.15e  %.15e  %.15e  %.15e %.15e  %.15e  %.15e  %.15e\n",
                k, x[4*k],x[4*k+1],x[4*k+2],x[4*k+3],
                   x[4*k+4],x[4*k+5],x[4*k+6],x[4*k+7]);
   }
}

int standard_membership_test ( void )
{
   int fail,n,dim,deg,nv;
   char ans;

   printf("\nReading a witness set ...\n");
   fail = read_witness_set(&n,&dim,&deg);
   nv = n - dim;

   if(verbose>0)  /* only in verbose mode */
   {
      printf("\nThe ambient dimension : %d.\n",n);
      printf("The dimension of the solution set : %d.\n",dim);
      printf("The degree of the solution set : %d.\n",deg);
      printf("\nDo you want to see the embedded system ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe system read :\n");
         fail = syscon_write_standard_system();
      }
      scanf("%c",&ans); /* skip end of line character */
      printf("\nDo you want to see the solutions ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe solutions read :\n");
         fail = solcon_write_standard_solutions();
      }
      scanf("%c",&ans); /* skip end of line character */
   }
   {
      int fail,onsys,onset,nbtasks=0;
      double tpt[2*nv];
      const double restol = 1.0e-8;
      const double homtol = 1.0e-6;
      printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
      printf("\nReading the coordinates of the test point x ...\n");
      standard_read_point(nv,tpt);
      fail = standard_homotopy_membership_test
               (1,nv,dim,restol,homtol,tpt,&onsys,&onset,nbtasks);
   }
   return 0;
}

int dobldobl_membership_test ( void )
{
   int fail,n,dim,deg,nv;
   char ans;

   printf("\nReading a witness set ...\n");
   fail = read_dobldobl_witness_set(&n,&dim,&deg);
   nv = n - dim;

   if(verbose>0)  /* only in verbose mode */
   {
      printf("\nThe ambient dimension : %d.\n",n);
      printf("The dimension of the solution set : %d.\n",dim);
      printf("The degree of the solution set : %d.\n",deg);
      printf("\nDo you want to see the embedded system ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe system read :\n");
         fail = syscon_write_dobldobl_system();
      }
      scanf("%c",&ans); /* skip end of line character */
      printf("\nDo you want to see the solutions ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe solutions read :\n");
         fail = solcon_write_dobldobl_solutions();
      }
      scanf("%c",&ans); /* skip end of line character */
   }
   {
      int fail,onsys,onset,nbtasks=0;
      double tpt[4*nv];
      const double restol = 1.0e-8;
      const double homtol = 1.0e-6;
      printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
      printf("\nReading the coordinates of the test point x ...\n");
      dobldobl_read_point(nv,tpt);
      fail = dobldobl_homotopy_membership_test
               (1,nv,dim,restol,homtol,tpt,&onsys,&onset,nbtasks);
   }
   return 0;
}

int quaddobl_membership_test ( void )
{
   int fail,n,dim,deg,nv;
   char ans;

   printf("\nReading a witness set ...\n");
   fail = read_quaddobl_witness_set(&n,&dim,&deg);
   nv = n - dim;

   if(verbose>0)  /* only in verbose mode */
   {
      printf("\nThe ambient dimension : %d.\n",n);
      printf("The dimension of the solution set : %d.\n",dim);
      printf("The degree of the solution set : %d.\n",deg);
      printf("\nDo you want to see the embedded system ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe system read :\n");
         fail = syscon_write_quaddobl_system();
      }
      scanf("%c",&ans); /* skip end of line character */
      printf("\nDo you want to see the solutions ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe solutions read :\n");
         fail = solcon_write_quaddobl_solutions();
      }
      scanf("%c",&ans); /* skip end of line character */
   }
   {
      int fail,onsys,onset,nbtasks=0;
      double tpt[8*nv];
      const double restol = 1.0e-8;
      const double homtol = 1.0e-6;
      printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
      printf("\nReading the coordinates of the test point x ...\n");
      quaddobl_read_point(nv,tpt);
      fail = quaddobl_homotopy_membership_test
               (1,nv,dim,restol,homtol,tpt,&onsys,&onset,nbtasks);
   }
   return 0;
}

int standard_Laurent_membership_test ( void )
{
   int fail,n,dim,deg,nv;
   char ans;

   printf("\nReading a witness set ...\n");
   fail = read_standard_Laurent_witness_set(&n,&dim,&deg);
   nv = n - dim;

   if(verbose>0)  /* only in verbose mode */
   {
      printf("\nThe ambient dimension : %d.\n",n);
      printf("The dimension of the solution set : %d.\n",dim);
      printf("The degree of the solution set : %d.\n",deg);
      printf("\nDo you want to see the embedded system ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe system read :\n");
         fail = syscon_write_standard_Laurent_system();
      }
      scanf("%c",&ans); /* skip end of line character */
      printf("\nDo you want to see the solutions ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe solutions read :\n");
         fail = solcon_write_standard_solutions();
      }
      scanf("%c",&ans); /* skip end of line character */
   }
   {
      int fail,onsys,onset,nbtasks;
      double tpt[2*nv];
      const double restol = 1.0e-8;
      const double homtol = 1.0e-6;
      printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
      printf("\nReading the coordinates of the test point x ...\n");
      standard_read_point(nv,tpt);
      fail = standard_Laurent_homotopy_membership_test
               (1,nv,dim,restol,homtol,tpt,&onsys,&onset,nbtasks);
   }
   return 0;
}

int dobldobl_Laurent_membership_test ( void )
{
   int fail,n,dim,deg,nv;
   char ans;

   printf("\nReading a witness set ...\n");
   fail = read_dobldobl_Laurent_witness_set(&n,&dim,&deg);
   nv = n - dim;

   if(verbose>0)  /* only in verbose mode */
   {
      printf("\nThe ambient dimension : %d.\n",n);
      printf("The dimension of the solution set : %d.\n",dim);
      printf("The degree of the solution set : %d.\n",deg);
      printf("\nDo you want to see the embedded system ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe system read :\n");
         fail = syscon_write_dobldobl_Laurent_system();
      }
      scanf("%c",&ans); /* skip end of line character */
      printf("\nDo you want to see the solutions ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe solutions read :\n");
         fail = solcon_write_dobldobl_solutions();
      }
      scanf("%c",&ans); /* skip end of line character */
   }
   {
      int fail,onsys,onset,nbtasks=0;
      double tpt[4*nv];
      const double restol = 1.0e-8;
      const double homtol = 1.0e-6;
      printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
      printf("\nReading the coordinates of the test point x ...\n");
      dobldobl_read_point(nv,tpt);
      fail = dobldobl_Laurent_homotopy_membership_test
               (1,nv,dim,restol,homtol,tpt,&onsys,&onset,nbtasks);
   }
   return 0;
}

int quaddobl_Laurent_membership_test ( void )
{
   int fail,n,dim,deg,nv;
   char ans;

   printf("\nReading a witness set ...\n");
   fail = read_quaddobl_Laurent_witness_set(&n,&dim,&deg);
   nv = n - dim;

   if(verbose>0)  /* only in verbose mode */
   {
      printf("\nThe ambient dimension : %d.\n",n);
      printf("The dimension of the solution set : %d.\n",dim);
      printf("The degree of the solution set : %d.\n",deg);
      printf("\nDo you want to see the embedded system ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe system read :\n");
         fail = syscon_write_quaddobl_Laurent_system();
      }
      scanf("%c",&ans); /* skip end of line character */
      printf("\nDo you want to see the solutions ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe solutions read :\n");
         fail = solcon_write_quaddobl_solutions();
      }
      scanf("%c",&ans); /* skip end of line character */
   }
   {
      int fail,onsys,onset,nbtasks=0;
      double tpt[8*nv];
      const double restol = 1.0e-8;
      const double homtol = 1.0e-6;
      printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
      printf("\nReading the coordinates of the test point x ...\n");
      quaddobl_read_point(nv,tpt);
      fail = quaddobl_Laurent_homotopy_membership_test
               (1,nv,dim,restol,homtol,tpt,&onsys,&onset,nbtasks);
   }
   return 0;
}
