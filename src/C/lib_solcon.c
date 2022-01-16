/* simple test program in C on the operations in the solutions container */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "solcon.h"

typedef struct list LIST;
struct list
{
   double rt;   /* real part of complex continuation parameter t */
   double it;   /* imaginary part of t */
   int m;       /* multiplicity of the solution vector */
   double *v;   /* solution is sequence of real and imaginary parts */
   double err;  /* error on the correction term from Newton's method */
   double rco;  /* inverse of estimate of condition number */
   double res;  /* residual of the solution vector */
   LIST *next;  /* next solution in the list */
};

LIST *push ( LIST *l, int n, double rt, double it, int m, double *v,
             double err, double rco, double res );
/*
 * DESCRIPTION :
 *   Pushes the solution to the front of the list l. */

void write_solution ( int n, double rt, double it, int m, double *v,
                      double err, double rco, double res );
/*
 * DESCRIPTION :
 *   Writes the n-dimensional solution to the screen. */

void write_solution_list ( LIST *l, int n );
/*
 * DESCRIPTION :
 *   Writes the list of n-dimensional solution vectors. */

void test_solution_retrieval ( int n );
/*
 * DESCRIPTION :
 *   Prompts the user for some solution number,
 *   retrieves and prints the solution vector,
 *   n is the complex dimension of the solution */

void test_replace_solution ( int n );
/*
 * DESCRIPTION :
 *   Prompts the user for some solution number,
 *   and replaces the solution. */

LIST *test_retrievals ( int len, int n );
/*
 * DESCRIPTION :
 *   Writes the solution list of length len and complex dimension n
 *   on standard output, returns the list of solutions. */

void test_additions ( int n, LIST *l );
/*
 * DESCRIPTION :
 *   Fills up the container with the solutions in the list l. */

void test_solution_container ( void );
/*
 * DESCRIPTION :
 *   tests main operations in the solutions container. */

void test_incremental_read_and_write ( void );
/*
 * DESCRIPTION :
 *   Tests incremental read and write of solutions. */

void test_solution_strings( void );
/*
 * DESCRIPTION :
 *   Tests writing of solutions to strings. */

void parse_solution_strings ( void );
/*
 * DESCRIPTION :
 *   Tests parsing of solutions from strings into container. */

void test_standard_container ( void );
/*
 * DESCRIPTION :
 *   Test operations in the container for standard solutions. */

void test_dobldobl_container ( void );
/*
 * DESCRIPTION :
 *   test operations in the container for double double solutions. */

void test_quaddobl_container ( void );
/*
 * DESCRIPTION :
 *   Test operations in the container for quad double solutions. */

void test_multprec_container ( void );
/*
 * DESCRIPTION :
 *   Test operations in the container for multiprecision solutions. */

void standard_next_retrievals ( void );
/*
 * DESCRIPTION :
 *   Prints the standard solution vectors with the retrieve_next. */

void test_standard_next_retrievals ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for solutions and then calls standard_next_retrievals. */

void retrieve_all_standard_solution_strings ( void );
/*
 * DESCRIPTION :
 *   Prints the strings of all standard solutions. */

void test_retrieve_all_standard_solution_strings ( void );
/* 
 * DESCRIPTION :
 *   Prompts the user for solutions and then calls 
 *   retrieve_all_standard_solution_strings. */

void test_read_solutions_from_file ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and then attempts to read
 *   solutions from file. */

int main(void)
{
   int choice;
   char ch;

   printf("\nMENU to test solution containers :\n");
   printf("  1. test with standard doubles;\n");
   printf("  2. test with double doubles;\n");
   printf("  3. test with quad doubles;\n");
   printf("  4. test with multiprecision numbers;\n");
   printf("  5. test retrieval of standard solutions with next.\n");
   printf("  6. test retrieval of all standard solution strings.\n");
   printf("  7. test reading solutions from a file with given name.\n");
   printf("Type 1, 2, 3, 4, 5, 6, or 7 to choose a test : ");
   scanf("%d",&choice);
   scanf("%c",&ch); /* skip new line character */
   printf("\n");

   adainit();

   if(choice == 1)
      test_standard_container();
   else if(choice == 2)
      test_dobldobl_container();
   else if(choice == 3)
      test_quaddobl_container();
   else if(choice == 4)
      test_multprec_container();
   else if(choice == 5)
      test_standard_next_retrievals();
   else if(choice == 6)
      test_retrieve_all_standard_solution_strings();
   else if(choice == 7)
      test_read_solutions_from_file();
   else
      printf("invalid choice, please try again\n");

   adafinal();

   return 0;
}

void retrieve_and_write ( int n, int k )
{
   int i,m,fail;
   double sol[2*n+5];

   fail = solcon_retrieve_standard_solution(n,k,&m,sol);

   if (fail == 1)
      printf("some failure occurred...\n");
   else
   { 
      printf("The solution vector :\n");
      for(i=0; i<n; i++)
         printf(" %.15e  %.15e\n", sol[2+2*i], sol[2+2*i+1]);
   }
}

void test_solution_retrieval ( int n )
{
   int k;

   printf("Give index to a solution : "); scanf("%d", &k);

   retrieve_and_write(n,k);
}

void test_replace_solution ( int n )
{
   int i,k,m,fail;
   double sol[2*n+5];

   printf("Give k : "); scanf("%d", &k);
   retrieve_and_write(n,k);

   printf("Setting solution vector to 1,2,3,...\n");
   sol[0] = 1.111;
   sol[1] = 2.222;
   m = 12345;
   sol[2*n+2] = 3.333;
   sol[2*n+3] = 4.444;
   sol[2*n+4] = 5.555; 
   for(i=0; i<2*n; i++) sol[2+i] = i;

   fail = solcon_replace_standard_solution(n,k,m,sol);

   if (fail == 1)
      printf("some failure occurred...\n");
   else
   { 
      printf("The changed solution retrieved again :\n");
      retrieve_and_write(n,k);
   }
}

void write_solution ( int n, double rt, double it, int m, double *v,
                      double err, double rco, double res )
{
   int i;

   printf("t : %.15e  %.15e\n", rt, it);
   printf("m : %d\n", m);
   printf("The solution for t :\n");
   for(i=0; i<n; i++)
      printf(" %.15e  %.15e\n", v[2*i], v[2*i+1]);
   printf("== err : %.3e = rco : %.3e = res : %.3e ==\n",err,rco,res);
}

LIST *push ( LIST *l, int n, double rt, double it, int m, double *v,
             double err, double rco, double res )
{
   int i;
   LIST *nl = (LIST*)calloc(1,sizeof(LIST));

   nl->rt = rt;
   nl->it = it;
   nl->m = m;
   nl->v = (double*)calloc(2*n,sizeof(double));
   for(i=0; i<2*n; i++)
      nl->v[i] = v[i];
   nl->err = err;
   nl->rco = rco;
   nl->res = res;
   nl->next = l;

   return nl;
}

void write_solution_list ( LIST *l, int n )
{
   LIST *p;
   int k;

   for(p=l,k=1; p!=NULL; p=p->next,k++)
   {
      printf("Solution %d :\n",k);
      write_solution(n,p->rt,p->it,p->m,p->v,p->err,p->rco,p->res);
   }
}

LIST *test_retrievals ( int len, int n )
{
   int i,k,m,fail;
   double sol[2*n+5],err,rco,res;
   LIST *l = NULL;

   for(k=1; k<=len; k++)
   {
      fail = solcon_retrieve_standard_solution(n,k,&m,sol);
      printf("Solution %d :\n", k);
      err = sol[2*n+2];
      rco = sol[2*n+3];
      res = sol[2*n+4];
      write_solution(n,sol[0],sol[1],m,sol+2,err,rco,res);
      l = push(l,n,sol[0],sol[1],m,sol+2,err,rco,res);
   }

   return l;
}

void test_additions ( int n, LIST *l )
{
   int i,k,m,fail;
   LIST *p;
   double sol[2*n+5];

   for(p=l,k=1; p!=NULL; p=p->next,k++)
   {
      sol[0] = p->rt;
      sol[1] = p->it;
      m = p->m;
      for(i=0;i<2*n;i++) sol[2+i] = p->v[i];
      sol[2*n+2] = p->err;
      sol[2*n+3] = p->rco;
      sol[2*n+4] = p->res;
      fail = solcon_append_standard_solution(n,m,sol);
   }
}

void test_solution_container ( void )
{
   int k,len,n,fail;
   double *c;
   LIST *l;

   fail = solcon_read_standard_solutions();
   printf("The solution list :\n");
   fail = solcon_write_standard_solutions();
   fail = solcon_number_of_standard_solutions(&len);
   printf("Number of solutions : %d\n",len);
   fail = solcon_dimension_of_standard_solutions(&n);
   printf("Dimension of vectors : %d\n",n);
   test_solution_retrieval(n);
   test_replace_solution(n);
   l = test_retrievals(len,n);
   printf("The solution list is reverse order : \n");
   write_solution_list(l,n);
   fail = solcon_clear_standard_solutions();
   test_additions(n,l);   
   printf("The reconstructed list in reverse order :\n");
   fail = solcon_write_standard_solutions(); 
}

void test_incremental_read_and_write ( void )
{
   int fail,*a,*b,len,dim,k,m,i;
   double *c,err,rco,res;
   char ans;

   fail = solcon_open_solution_input_file();
   printf("\n");
   fail = solcon_create_solution_output_file();
   printf("\nShould the file be scanned for a banner ? (y/n) ");
   scanf("%c",&ans);
   if (ans == 'y') fail = solcon_scan_solution_banner();
   fail = solcon_read_solution_dimensions(&len,&dim);
   printf("#solutions : %d  dimension : %d\n",len,dim);
   fail = solcon_write_solution_dimensions(len,dim);
   c = (double*)calloc(2*dim+5,sizeof(double));
   k = 0;
   for(i=0; i<len; i++)
   {
      fail = solcon_read_next_solution(dim,&m,c);
      printf("Solution %d :\n", i+1);
      err = c[2*dim+2]; rco = c[2*dim+3]; res = c[2*dim+4];
      write_solution(dim,c[0],c[1],m,c+2,err,rco,res);
      fail = solcon_write_next_solution(&k,dim,m,c);
      printf("number of solutions written to file : %d\n",k);
   }
   fail = solcon_close_solution_input_file(0);
   fail = solcon_close_solution_output_file();
}

void test_solution_strings ( void )
{
   int len,fail,k,n,n0,n1,n2;

   fail = solcon_read_standard_solutions();
   fail = solcon_number_of_standard_solutions(&len);
   printf("Number of solutions : %d\n",len);
   printf("Give solution number : "); scanf("%d",&k);
   fail = solcon_length_standard_solution_string(k,&n);
   fail = solcon_length_solution_intro(k,&n0);
   fail = solcon_length_solution_vector(k,&n1);
   fail = solcon_length_solution_diagnostics(k,&n2);
   printf("Length of solution %d : %d = %d + %d + %d\n",k,n,n0,n1,n2);
   {
      char s[n+1],s0[n0+1],s1[n1+1],s2[n2+1];
      fail = solcon_write_standard_solution_string(k,n,s);
      printf("Solution %d as string :\n%s\n",k,s);
      fail = solcon_write_solution_intro(k,n0,s0);
      printf("introduction to solution %d :\n%s",k,s0);
      fail = solcon_write_solution_vector(k,n1,s1);
      printf("vector in the solution %d :\n%s",k,s1);
      fail = solcon_write_solution_diagnostics(k,n2,s2);
      printf("diagnostics of the solution %d :\n%s",k,s2);
   }
}

void parse_solution_strings ( void )
/*
 * Reads solutions from file into the container,
 * stores all solutions as strings, clears the container,
 * and fills the container using the solution strings. */
{
   int fail,len,dim;

   fail = solcon_read_standard_solutions();
   fail = solcon_number_of_standard_solutions(&len);
   fail = solcon_dimension_of_standard_solutions(&dim);
   {
      char* solutions[len];
      int nc[len];
      int k;

      for(k=0; k<len; k++)
      {
         fail = solcon_length_standard_solution_string(k+1,&nc[k]);
         solutions[k] = (char*)calloc(nc[k]+1,sizeof(char));
         fail = solcon_write_standard_solution_string
                  (k+1,nc[k]+1,solutions[k]);
      }
      printf("\nCheck for all solution strings ...\n");
      do
      {
         printf("Give a number less than %d (0 to quit) : ",len+1);
         scanf("%d",&k);
         if((k <= 0) || (k > len)) break;
         printf("Solution %d :\n%s\n",k,solutions[k-1]);
      } while (k > 0);
      printf("\nClearing the solution container ...\n");
      solcon_clear_standard_solutions();
      printf("\nAppending the solution strings to the container ...\n");
      for(k=0; k<len; k++)
         fail = solcon_append_standard_solution_string
                   (dim,nc[k]+1,solutions[k]);
      printf("\nThe solutions in the container :\n");
      solcon_write_standard_solutions();
   }
}

void test_standard_container ( void )
{
   int option;
   char ch;

   printf("\nMENU to test operations in solutions container :\n");
   printf("  1. test main operations in the container;\n");
   printf("  2. incremental read/write of solutions from/to file;\n");
   printf("  3. test the writing of solutions to strings;\n");
   printf("  4. test parsing of solutions from strings into container.\n");
   printf("Type 1, 2, 3, or 4 to make your choice : "); scanf("%d",&option);
   scanf("%c",&ch); /* skip new line character */
   printf("\n");

   if(option == 1)
      test_solution_container();
   else if(option == 2)
      test_incremental_read_and_write();
   else if(option == 3)
      test_solution_strings();
   else if(option == 4)
      parse_solution_strings();
   else
      printf("%d is wrong choice.  Please try again...\n", option);
}

void test_dobldobl_container ( void )
{
   int fail,len,dim,k,nc;
   char *solution;
   char ans;

   fail = solcon_read_dobldobl_solutions();
   fail = solcon_write_dobldobl_solutions();
   fail = solcon_number_of_dobldobl_solutions(&len);
   printf("Number of solutions : %d\n",len);
   fail = solcon_dimension_of_dobldobl_solutions(&dim);
   printf("Dimension of vectors : %d\n",dim);

   printf("\nTesting solution string retrieval ...\n\n");
   printf("Give an index to a solution : ");
   scanf("%d",&k); scanf("%c",&ans); /* skip new line */
   fail = solcon_length_dobldobl_solution_string(k,&nc);
   solution = (char*)calloc(nc+1,sizeof(char));
   fail = solcon_write_dobldobl_solution_string(k,nc+1,solution);
   printf("solution %d :\n%s\n",k,solution);

   printf("hit enter to continue with append ...\n");
   scanf("%c",&ans);
   fail = solcon_append_dobldobl_solution_string(dim,nc+1,solution);
   printf("solutions after the append :\n");
   fail = solcon_write_dobldobl_solutions();
}

void test_quaddobl_container ( void )
{
   int fail,len,dim,k,nc;
   char *solution;
   char ans;

   fail = solcon_read_quaddobl_solutions();
   fail = solcon_write_quaddobl_solutions();
   fail = solcon_number_of_quaddobl_solutions(&len);
   printf("Number of solutions : %d\n",len);
   fail = solcon_dimension_of_quaddobl_solutions(&dim);
   printf("Dimension of vectors : %d\n",dim);

   printf("\nTesting solution string retrieval ...\n\n");
   printf("Give an index to a solution : ");
   scanf("%d",&k); scanf("%c",&ans);
   fail = solcon_length_quaddobl_solution_string(k,&nc);
   solution = (char*)calloc(nc+1,sizeof(char));
   fail = solcon_write_quaddobl_solution_string(k,nc+1,solution);
   printf("solution %d :\n%s\n",k,solution);

   printf("hit enter to continue with append ...\n");
   scanf("%c",&ans);
   fail = solcon_append_quaddobl_solution_string(dim,nc+1,solution);
   printf("solutions after the append :\n");
   fail = solcon_write_quaddobl_solutions();
}

void test_multprec_container ( void )
{
   int fail,len,dim,k,nc;
   char *solution;
   char ans;

   fail = solcon_read_multprec_solutions();
   fail = solcon_write_multprec_solutions();
   fail = solcon_number_of_multprec_solutions(&len);
   printf("Number of solutions : %d\n",len);
   fail = solcon_dimension_of_multprec_solutions(&dim);
   printf("Dimension of vectors : %d\n",dim);

   printf("\nTesting solution string retrieval ...\n\n");
   printf("Give an index to a solution : ");
   scanf("%d",&k); scanf("%c",&ans);
   fail = solcon_length_multprec_solution_string(k,&nc);
   solution = (char*)calloc(nc+1,sizeof(char));
   fail = solcon_write_multprec_solution_string(k,nc+1,solution);
   printf("solution %d :\n%s\n",k,solution);

   printf("hit enter to continue with append ...\n");
   scanf("%c",&ans);
   fail = solcon_append_multprec_solution_string(dim,nc+1,solution);
   printf("solutions after the append :\n");
   fail = solcon_write_multprec_solutions();
}

void standard_next_retrievals ( void )
{
   int dim,len,fail;
   fail = solcon_dimension_of_standard_solutions(&dim);
   fail = solcon_number_of_standard_solutions(&len);
   printf("The container has %d solutions of dimension %d.\n",len,dim);
   {
      int i,k,m;
      double sol[2*dim+5];
      do
      {
         fail = solcon_retrieve_next_standard_solution(dim,&k,&m,sol);
         if(fail != 0)
            printf("return value of k is %d\n",k);
         else
         { 
            printf("The solution vector %d :\n",k);
            for(i=0; i<dim; i++)
               printf(" %.15e  %.15e\n", sol[2+2*i], sol[2+2*i+1]);
         }
      }
      while(k > 0);
   }
}

void test_standard_next_retrievals ( void )
{
   int fail = solcon_read_standard_solutions();
   standard_next_retrievals();
}

void retrieve_all_standard_solution_strings ( void )
{
   int dim,len,fail;
   fail = solcon_dimension_of_standard_solutions(&dim);
   fail = solcon_number_of_standard_solutions(&len);
   printf("The container has %d solutions of dimension %d.\n",len,dim);
   {
      int k,len;
      do
      {
         fail = solcon_length_current_standard_solution_string(&k,&len);
         if(fail != 0)
            printf("return value of k is %d\n",k);
         else
         { 
            char sol[len+1];
            fail = solcon_write_current_standard_solution_string(&k,len,sol);
            printf("The solution vector %d :\n%s\n",k, sol);
            fail = solcon_move_current_standard_to_next(&k);
         }
      }
      while(k > 0);
   }
}

void test_retrieve_all_standard_solution_strings ( void )
{
   int fail = solcon_read_standard_solutions();
   retrieve_all_standard_solution_strings();
}

void test_read_solutions_from_file ( void )
{
   char name[80];
   int nbc,fail,len;

   printf("Give the name of the file : ");
   scanf("%s",name);
   printf("Reading solutions from file with name \"%s\" ...\n",name);

   nbc = strlen(name);
   printf("-> number of characters in name : %d\n",nbc);

   printf("\ncalling solcon_read_standard_solutions_from_file ...\n");
   fail = solcon_read_standard_solutions_from_file(nbc,name);
   fail = solcon_number_of_standard_solutions(&len);
   printf("Number of solutions : %d\n",len);
}
