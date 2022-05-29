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

int input_output_on_files ( int precision );
/*
 * DESCRIPTION :
 *   This function calls the blackbox solver reading a polynomial
 *   system from file and writing the solutions to screen.
 *   Three levels of precision are supported, depending on precision:
 *   0 for standard double, 1 for double double, 2 for quad double. */

int Laurent_input_output_on_files ( int precision );
/*
 * DESCRIPTION :
 *   This function calls the blackbox solver reading a Laurent
 *   polynomial system from file and writing the solutions to file. 
 *   Three levels of precision are supported, depending on precision:
 *   0 for standard double, 1 for double double, 2 for quad double. */

int interactive_input_output ( int precision );
/*
 * DESCRIPTION :
 *   The user is prompted to give a polynomial system which is solved
 *   with the blackbox solver in double, double double, or quad double
 *   precision, depending on the value of precision: 0, 1, or 2.
 *   Solutions are written to the screen, so no files are read or created. */

int interactive_Laurent_input_output ( int precision );
/*
 * DESCRIPTION :
 *   The user is prompted to give a polynomial system which is solved
 *   with the blackbox solver in double, double double, or quad double
 *   precision, depending on the value of precision: 0, 1, or 2. 
 *   Solutions are written to the screen, so no files are read or created. */

int prompt_for_precision ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for the level of precision and returns
 *   0 for standard double precision, 1 for double double precision, or
 *   2 for quad double precision. */

void greetings ( void );
/*
 * DESCRIPTION :
 *   Displays the version string of PHCpack. */

int main ( int argc, char *argv[] )
{
   int choice,fail,prec;

   adainit();

   greetings();

   prec = prompt_for_precision();
 
   printf("\nWelcome to the blackbox solver in PHCpack :\n");
   printf("  1. the input polynomial system will be typed in; or\n");
   printf("  2. the system typed in is a Laurent polynomial system; or\n");
   printf("  3. an input file contains the input polynomial system; or\n");
   printf("  4. the input Laurent system is on an input file.\n");
   printf("Type 1, 2, 3, or 4 to make your choice : ");
   scanf("%d",&choice);
   if(choice == 1)
      fail = interactive_input_output(prec);
   else if (choice == 2)
      fail = interactive_Laurent_input_output(prec);
   else if (choice == 3)
      fail = input_output_on_files(prec);
   else
      fail = Laurent_input_output_on_files(prec);

   adafinal();

   return 0;
}

int input_output_on_files ( int precision )
{
   int fail,rc,nrc,nbtasks,vrb,i,focusmv;
   char rocos[1024];

   if(precision == 0)
   {
      fail = syscon_read_standard_system();
      printf("\nThe system in the container : \n");
      fail = syscon_write_standard_system();
      printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
      printf("\nFocus on mixed volumes (1 = yes, 0 = no) : ");
      scanf("%d",&focusmv);
      printf("\nGive the verbose level : "); scanf("%d",&vrb);
      fail = solve_standard_system(&rc,0,&nrc,rocos,nbtasks,focusmv,vrb);
      printf("\nROOT COUNTS :\n");
      for(i=0; i<nrc; i++) printf("%c",rocos[i]);
      printf("\nThe root count : %d\n",rc);
      printf("\nThe solutions :\n");
      fail = solcon_write_standard_solutions();
   }
   else if(precision == 1)
   {
      fail = syscon_read_dobldobl_system();
      printf("\nThe system in the container : \n");
      fail = syscon_write_dobldobl_system();
      printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
      printf("\nGive the verbose level : "); scanf("%d",&vrb);
      fail = solve_dobldobl_system(&rc,0,&nrc,rocos,nbtasks,vrb);
      printf("\nROOT COUNTS :\n");
      for(i=0; i<nrc; i++) printf("%c",rocos[i]);
      printf("\nThe root count : %d\n",rc);
      printf("\nThe solutions :\n");
      fail = solcon_write_dobldobl_solutions();
   }
   else if(precision == 2)
   {
      fail = syscon_read_quaddobl_system();
      printf("\nThe system in the container : \n");
      fail = syscon_write_quaddobl_system();
      printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
      printf("\nGive the verbose level : "); scanf("%d",&vrb);
      fail = solve_quaddobl_system(&rc,0,&nrc,rocos,nbtasks,vrb);
      printf("\nROOT COUNTS :\n");
      for(i=0; i<nrc; i++) printf("%c",rocos[i]);
      printf("\nThe root count : %d\n",rc);
      printf("\nThe solutions :\n");
      fail = solcon_write_quaddobl_solutions();
   }

   return fail;
}

int Laurent_input_output_on_files ( int precision )
{
   int fail,rc,nrc,nbtasks,focusmv,vrb,i;
   char rocos[1024];

   if(precision == 0)
   {
      fail = syscon_read_standard_Laurent_system();
      printf("\nThe system in the container : \n");
      fail = syscon_write_standard_Laurent_system();
      printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
      printf("\nFocus on mixed volumes (1 = yes, 0 = no) : ");
      scanf("%d",&focusmv);
      printf("\nGive the verbose level : "); scanf("%d",&vrb);
      fail = solve_standard_Laurent_system
               (&rc,0,&nrc,rocos,nbtasks,focusmv,vrb);
      printf("\nROOT COUNTS :\n");
      for(i=0; i<nrc; i++) printf("%c",rocos[i]);
      printf("\nThe root count : %d\n",rc);
      printf("\nThe solutions :\n");
      fail = solcon_write_standard_solutions();
   }
   else if(precision == 1)
   {
      fail = syscon_read_dobldobl_Laurent_system();
      printf("\nThe system in the container : \n");
      fail = syscon_write_dobldobl_Laurent_system();
      printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
      printf("\nGive the verbose level : "); scanf("%d",&vrb);
      fail = solve_dobldobl_Laurent_system(&rc,0,&nrc,rocos,nbtasks,vrb);
      printf("\nROOT COUNTS :\n");
      for(i=0; i<nrc; i++) printf("%c",rocos[i]);
      printf("\nThe root count : %d\n",rc);
      printf("\nThe solutions :\n");
      fail = solcon_write_dobldobl_solutions();
   }
   else if(precision == 2)
   {
      fail = syscon_read_quaddobl_Laurent_system();
      printf("\nThe system in the container : \n");
      fail = syscon_write_quaddobl_Laurent_system();
      printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
      printf("\nGive the verbose level : "); scanf("%d",&vrb);
      fail = solve_quaddobl_Laurent_system(&rc,0,&nrc,rocos,nbtasks,vrb);
      printf("\nROOT COUNTS :\n");
      for(i=0; i<nrc; i++) printf("%c",rocos[i]);
      printf("\nThe root count : %d\n",rc);
      printf("\nThe solutions :\n");
      fail = solcon_write_quaddobl_solutions();
   }
   return fail;
}

int standard_interactive_input_output ( void )
{
   int n,fail,k,nc,i,rc,nrc,len,nbtasks,focusmv,vrb;
   char rocos[1024];
   char ch,p[800];

   printf("\nGive the number of polynomials in the system : ");
   scanf("%d",&n);

   fail = syscon_initialize_number_of_standard_polynomials(n);
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
      fail = syscon_store_standard_polynomial(nc,n,k,p);
   }
   printf("The system in the container : \n");
   syscon_write_standard_system();
   printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
   printf("\nFocus on polynomials (1 = yes, 0 = no) : ");
   scanf("%d",&focusmv);
   printf("\nGive the verbose level : "); scanf("%d",&vrb);
   fail = solve_standard_system(&rc,0,&nrc,rocos,nbtasks,focusmv,vrb);
   printf("\nTHE ROOT COUNTS :\n");
   for(i=0; i<nrc; i++) printf("%c",rocos[i]);
   printf("\nThe root count : %d\n",rc);
   /* printf("\nThe solutions :\n");
   fail = solcon_write_standard_solutions(); */
   printf("interactive selection of solutions... \n");
   do
   {
      printf("\nGive an index to a solution (-1 to exit) : ");
      scanf("%d",&k);
      if(k < 0) break;
      fail = solcon_length_standard_solution_string(k,&len);
      {
         char s[len];
         fail = solcon_write_standard_solution_string(k,len,s);
         printf("\nSolution %d : \n%s",k,s);
      } 
   } while (k >= 0);

   return fail;
}

int dobldobl_interactive_input_output ( void )
{
   int n,fail,k,nc,i,rc,nrc,len,nbtasks,vrb;
   char rocos[1024];
   char ch,p[800];

   printf("\nGive the number of polynomials in the system : ");
   scanf("%d",&n);

   fail = syscon_initialize_number_of_dobldobl_polynomials(n);
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
      fail = syscon_store_dobldobl_polynomial(nc,n,k,p);
   }
   printf("The system in the container : \n");
   syscon_write_dobldobl_system();
   printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
   printf("\nGive the verbose level : "); scanf("%d",&vrb);
   fail = solve_dobldobl_system(&rc,0,&nrc,rocos,nbtasks,vrb);
   printf("\nThe root count : %d\n",rc);
   /* printf("\nThe solutions :\n");
   fail = solcon_write_dobldobl_solutions(); */
   printf("interactive selection of solutions... \n");
   do
   {
      printf("\nGive an index to a solution (-1 to exit) : ");
      scanf("%d",&k);
      if(k < 0) break;
      fail = solcon_length_dobldobl_solution_string(k,&len);
      {
         char s[len];
         fail = solcon_write_dobldobl_solution_string(k,len,s);
         printf("\nSolution %d : \n%s",k,s);
      } 
   } while (k >= 0);

   return fail;
}

int quaddobl_interactive_input_output ( void )
{
   int n,fail,k,nc,i,rc,nrc,len,nbtasks,vrb;
   char rocos[1024];
   char ch,p[800];

   printf("\nGive the number of polynomials in the system : ");
   scanf("%d",&n);

   fail = syscon_initialize_number_of_quaddobl_polynomials(n);
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
      fail = syscon_store_quaddobl_polynomial(nc,n,k,p);
   }
   printf("The system in the container : \n");
   syscon_write_quaddobl_system();
   printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
   printf("\nGive the verbose level : "); scanf("%d",&vrb);
   fail = solve_quaddobl_system(&rc,0,&nrc,rocos,nbtasks,vrb);
   printf("\nThe root count : %d\n",rc);
   /* printf("\nThe solutions :\n");
   fail = solcon_write_quaddobl_solutions(); */
   printf("interactive selection of solutions... \n");
   do
   {
      printf("\nGive an index to a solution (-1 to exit) : ");
      scanf("%d",&k);
      if(k < 0) break;
      fail = solcon_length_quaddobl_solution_string(k,&len);
      {
         char s[len];
         fail = solcon_write_quaddobl_solution_string(k,len,s);
         printf("\nSolution %d : \n%s",k,s);
      } 
   } while (k >= 0);

   return fail;
}

int interactive_input_output ( int precision )
{
   if(precision == 0)
      standard_interactive_input_output();
   else if(precision == 1)
      dobldobl_interactive_input_output();
   else if(precision == 2)
      quaddobl_interactive_input_output();

   return 0;
}

int standard_interactive_Laurent_input_output ( void )
{
   int n,fail,k,nc,nrc,i,rc,len,nbtasks,focusmv,vrb;
   char rocos[1024];
   char ch,p[800];

   printf("\nGive the number of polynomials in the system : ");
   scanf("%d",&n);

   fail = syscon_initialize_number_of_standard_Laurentials(n);
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
      fail = syscon_store_standard_Laurential(nc,n,k,p);
   }
   printf("The system in the container : \n");
   syscon_write_standard_Laurent_system();
   printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
   printf("\nFocus on mixed volumes (1 = yes, 0 = no) : ");
   scanf("%d",&focusmv);
   printf("\nGive the verbose level : "); scanf("%d",&vrb);
   fail = solve_standard_Laurent_system
            (&rc,0,&nrc,rocos,nbtasks,focusmv,vrb);
   printf("\nThe root count : %d\n",rc);
   /* printf("\nThe solutions :\n");
   fail = solcon_write_standard_solutions(); */
   printf("interactive selection of solutions... \n");
   do
   {
      printf("\nGive an index to a solution (-1 to exit) : ");
      scanf("%d",&k);
      if(k < 0) break;
      fail = solcon_length_standard_solution_string(k,&len);
      {
         char s[len];
         fail = solcon_write_standard_solution_string(k,len,s);
         printf("\nSolution %d : \n%s",k,s);
      } 
   } while (k >= 0);

   return fail;
}

int dobldobl_interactive_Laurent_input_output ( void )
{
   int n,fail,k,nc,i,rc,nrc,len,nbtasks,vrb;
   char rocos[1024];
   char ch,p[800];

   printf("\nGive the number of polynomials in the system : ");
   scanf("%d",&n);

   fail = syscon_initialize_number_of_dobldobl_Laurentials(n);
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
      fail = syscon_store_dobldobl_Laurential(nc,n,k,p);
   }
   printf("The system in the container : \n");
   syscon_write_dobldobl_Laurent_system();
   printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
   printf("\nGive the verbose level : "); scanf("%d",&vrb);
   fail = solve_dobldobl_Laurent_system(&rc,0,&nrc,rocos,nbtasks,vrb);
   printf("\nThe root count : %d\n",rc);
   /* printf("\nThe solutions :\n");
   fail = solcon_write_dobldobl_solutions(); */
   printf("interactive selection of solutions... \n");
   do
   {
      printf("\nGive an index to a solution (-1 to exit) : ");
      scanf("%d",&k);
      if(k < 0) break;
      fail = solcon_length_dobldobl_solution_string(k,&len);
      {
         char s[len];
         fail = solcon_write_dobldobl_solution_string(k,len,s);
         printf("\nSolution %d : \n%s",k,s);
      } 
   } while (k >= 0);

   return fail;
}

int quaddobl_interactive_Laurent_input_output ( void )
{
   int n,fail,k,nc,i,rc,nrc,len,nbtasks,vrb;
   char rocos[1024];
   char ch,p[800];

   printf("\nGive the number of polynomials in the system : ");
   scanf("%d",&n);

   fail = syscon_initialize_number_of_quaddobl_Laurentials(n);
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
      fail = syscon_store_quaddobl_Laurential(nc,n,k,p);
   }
   printf("The system in the container : \n");
   syscon_write_quaddobl_Laurent_system();
   printf("\nGive the number of tasks : "); scanf("%d",&nbtasks);
   printf("\nGive the verbose level : "); scanf("%d",&vrb);
   fail = solve_quaddobl_Laurent_system(&rc,0,&nrc,rocos,nbtasks,vrb);
   printf("\nThe root count : %d\n",rc);
   /* printf("\nThe solutions :\n");
   fail = solcon_write_quaddobl_solutions(); */
   printf("interactive selection of solutions... \n");
   do
   {
      printf("\nGive an index to a solution (-1 to exit) : ");
      scanf("%d",&k);
      if(k < 0) break;
      fail = solcon_length_quaddobl_solution_string(k,&len);
      {
         char s[len];
         fail = solcon_write_quaddobl_solution_string(k,len,s);
         printf("\nSolution %d : \n%s",k,s);
      } 
   } while (k >= 0);

   return fail;
}

int interactive_Laurent_input_output ( int precision )
{
   if(precision == 0)
      standard_interactive_Laurent_input_output();
   else if(precision == 1)
      dobldobl_interactive_Laurent_input_output();
   else if(precision == 2)
      quaddobl_interactive_Laurent_input_output();

   return 0;
}

void read_poly ( int *nc, char p[] )
{
   int c, i = 0;
 
   while ((c = getchar()) != EOF && c != ';') p[i++] = c;
   p[i++] = ';';
   p[i++] = '\0';
   *nc = i;
}

int prompt_for_precision ( void )
{
   char c;

   while(1)
   {
      printf("\nMENU for the working precision :\n");
      printf("  0. standard double precision;\n");
      printf("  1. double double precision;\n");
      printf("  2. quad double precision.\n");
      printf("Type 0, 1, or 2 to select the precision : ");
      c = getchar();
      if(c == '0')
         return 0;
      else if(c == '1')
         return 1;
      else if(c == '2')
         return 2;
      else
      {
         printf("Wrong selection, please try again ...\n");
         scanf("%c",&c); /* skip newline symbol */
      }
   }
}

void greetings ( void )
{
   int fail,n;
   char s[40];

   fail = version_string(&n,s);

   printf("Running %s ...\n",s);
}
