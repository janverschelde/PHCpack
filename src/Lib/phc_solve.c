/* Prompts the user for a system and solves it calling the blackbox
 * solver in PHCpack. */

#include <stdio.h>
#include <stdlib.h>
#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"

void read_poly ( int *nc, char p[] );
/* 
 * DESCRIPTION :
 *   Reads a polynomial as a string of characters.
 *   The number of characters is returned in nc. */

int input_output_on_files ( void );
/*
 * DESCRIPTION :
 *   This function calls the blackbox solver in PHCpack reading a
 *   polynomial system from file and writing the solutions to file. */

int Laurent_input_output_on_files ( void );
/*
 * DESCRIPTION :
 *   This function calls the blackbox solver in PHCpack reading a Laurent
 *   polynomial system from file and writing the solutions to file. */

int interactive_input_output ( void );
/*
 * DESCRIPTION :
 *   The user is prompted to give a polynomial system which is solved
 *   with the blackbox solver of PHCpack.  The solutions are written
 *   to the screen, so no files are read or created. */

int interactive_Laurent_input_output ( void );
/*
 * DESCRIPTION :
 *   The user is prompted to give a polynomial system which is solved
 *   with the blackbox solver of PHCpack.  The solutions are written
 *   to the screen, so no files are read or created. */

void greetings( void );
/*
 * DESCRIPTION :
 *   Displays the version string of PHCpack. */

int main ( int argc, char *argv[] )
{
   int choice,fail;

   adainit();

   greetings();
 
   printf("\nWelcome to the blackbox solver in PHCpack :\n");
   printf("  1. the input polynomial system will be typed in; or\n");
   printf("  2. the system typed in is a Laurent polynomial system; or\n");
   printf("  3. an input file contains the input polynomial system; or\n");
   printf("  4. the input Laurent system is on an input file.\n");
   printf("Type 1, 2, 3, or 4 to make your choice : ");
   scanf("%d",&choice);
   if(choice == 1)
      fail = interactive_input_output();
   else if (choice == 2)
      fail = interactive_Laurent_input_output();
   else if (choice == 3)
      fail = input_output_on_files();
   else
      fail = Laurent_input_output_on_files();

   adafinal();

   return 0;
}

int input_output_on_files ( void )
{
   int fail,rc;

   fail = syscon_read_system();
   printf("\nThe system in the container : \n");
   fail = syscon_write_system();
   fail = solve_system(&rc);
   printf("\nThe root count : %d\n",rc);
   printf("\nThe solutions :\n");
   fail = solcon_write_solutions();

   return fail;
}

int Laurent_input_output_on_files ( void )
{
   int fail,rc;

   fail = syscon_read_Laurent_system();
   printf("\nThe system in the container : \n");
   fail = syscon_write_Laurent_system();
   fail = solve_Laurent_system(&rc,0); /* not silent by default */
   printf("\nThe root count : %d\n",rc);
   printf("\nThe solutions :\n");
   fail = solcon_write_solutions();

   return fail;
}

int interactive_input_output ( void )
{
   int n,fail,k,nc,i,rc,len;
   char ch,p[800];

   printf("\nGive the number of polynomials in the system : ");
   scanf("%d",&n);

   fail = syscon_initialize_number(n);
   printf("\nReading %d polynomials, ",n);
   printf("terminate each with ; (semicolon)...\n");
   for(k=1; k<=n; k++)
   { 
      printf("-> polynomial %d : ",k);
      ch = getchar();    /* skip previous newline symbol */
      read_poly(&nc,p);
     /* printf("  p[%d] = ",k); 
        for(i=0; i<nc; i++) printf("%c",p[i]);
        printf("\n"); */
      fail = syscon_store_polynomial(nc,n,k,p);
   }
   printf("The system in the container : \n");
   syscon_write_system();
   fail = solve_system(&rc);
   printf("\nThe root count : %d\n",rc);
   /* printf("\nThe solutions :\n");
   fail = solcon_write_solutions(); */
   printf("interactive selection of solutions... \n");
   do
   {
      printf("\nGive an index to a solution (-1 to exit) : ");
      scanf("%d",&k);
      if(k < 0) break;
      fail = solcon_length_solution_string(k,&len);
      {
         char s[len];
         fail = solcon_write_solution_string(k,len,s);
         printf("\nSolution %d : \n%s",k,s);
      } 
   } while (k >= 0);

   return fail;
}

int interactive_Laurent_input_output ( void )
{
   int n,fail,k,nc,i,rc,len;
   char ch,p[800];

   printf("\nGive the number of polynomials in the system : ");
   scanf("%d",&n);

   fail = syscon_initialize_number_of_Laurentials(n);
   printf("\nReading %d polynomials, ",n);
   printf("terminate each with ; (semicolon)...\n");
   for(k=1; k<=n; k++)
   { 
      printf("-> polynomial %d : ",k);
      ch = getchar();    /* skip previous newline symbol */
      read_poly(&nc,p);
     /* printf("  p[%d] = ",k); 
        for(i=0; i<nc; i++) printf("%c",p[i]);
        printf("\n"); */
      fail = syscon_store_Laurential(nc,n,k,p);
   }
   printf("The system in the container : \n");
   syscon_write_Laurent_system();
   fail = solve_Laurent_system(&rc,0); /* not silent by default */
   printf("\nThe root count : %d\n",rc);
   /* printf("\nThe solutions :\n");
   fail = solcon_write_solutions(); */
   printf("interactive selection of solutions... \n");
   do
   {
      printf("\nGive an index to a solution (-1 to exit) : ");
      scanf("%d",&k);
      if(k < 0) break;
      fail = solcon_length_solution_string(k,&len);
      {
         char s[len];
         fail = solcon_write_solution_string(k,len,s);
         printf("\nSolution %d : \n%s",k,s);
      } 
   } while (k >= 0);

   return fail;
}

void read_poly ( int *nc, char p[] )
{
   int c, i = 0;
 
   while ((c = getchar()) != EOF && c != ';') p[i++] = c;
   p[i++] = ';';
   p[i++] = '\0';
   *nc = i;
}

void greetings( void )
{
   int fail,n;
   char s[40];

   fail = version_string(&n,s);

   printf("Running %s ...\n",s);
}
