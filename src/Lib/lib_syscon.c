/* simple test program in C on the operations in the systems container */

#include <stdio.h>
#include <stdlib.h>
#include "syscon.h"

typedef struct list LIST;
struct list
{
   double rc;  /* real part of coefficient */
   double ic;  /* imaginary part of coefficient */
   int  *exp;  /* exponent vector */
   LIST *next; /* next monomial in the list */
};

LIST *push ( LIST *l, int n, double *c, int *e );
/* pushes the monomial to the front of the list l */

void write_monomial_list ( LIST *l, int n );
/*
 * DESCRIPTION :
 *   Writes the list of monomials,
 *   the number of variables is in n. */

void test_retrievals ( int neq, int nvr, LIST *p[] );
/*
 * DESCRIPTION :
 *   Prints the system in the container,
 *   and returns the lists of terms in p.
 *
 * ON ENTRY :
 *   neq     number of equations in the system;
 *   nvr     number of variables in the system.
 *
 * ON RETURN :
 *   p       an array of lists of terms.  */

void test_additions ( int neq, int nvr, LIST *p[] );
/* 
 * DESCRIPTION :
 *   Fills up the container with the terms in p.
 *   The number of equations is in neq. */

void test_symbol_table ( int *nvr );
/*
 * DESCRIPTION :
 *   Tests operations in the symbol table.
 *   On return in nvr is the number of symbols. */

void test_standard_container ( void );
/* test operations of standard system container */

void test_dobldobl_container ( void );
/* test operations of the double double system container */

void test_quaddobl_container ( void );
/* test operations of the quad double system container */

void test_multprec_container ( void );
/* test operations of the quad double system container */

void show_random_system ( void );
/* prompts the user for the dimensions of a random system
 * generates it and writes it to screen */

int main(void)
{
   int choice;

   printf("\nMENU for testing the systems containers :\n");
   printf("  0. show a random polynomial system;\n");
   printf("  1. test with standard double coefficients;\n");
   printf("  2. test with double double coefficients;\n");
   printf("  3. test with quad double coefficients;\n");
   printf("  4. test with multiprecision coefficients.\n");
   printf("Type 1, 2, 3, or 4 to select the precision : ");
   scanf("%d",&choice);

   adainit();
   if(choice == 0)
      show_random_system();
   else if(choice == 1)
      test_standard_container();
   else if(choice == 2)
      test_dobldobl_container();
   else if(choice == 3)
      test_quaddobl_container();
   else if(choice == 4)
      test_multprec_container();
   else
      printf("invalid choice, please try again...\n");
   adafinal();

   return 0;
}

void test_retrievals ( int neq, int nvr, LIST *p[] )
{
   int i,j,k,nt,fail;
   double c[2];
   int d[nvr];

   for(i=1; i<=neq; i++)
   {
      fail = syscon_number_of_standard_terms(i, &nt);
      printf("  #terms in polynomial %d : %d\n", i, nt);

      p[i-1] = NULL;
      for(j=1; j<=nt; j++)
      {
         fail = syscon_retrieve_standard_term(i,j,nvr,d,c);
         printf("%22.15e  %22.15e", c[0], c[1]);
         for (k=0; k<nvr; k++) printf(" %d", d[k]);
         printf("\n");
         p[i-1] = push(p[i-1],nvr,c,d);
      }
   }
}

void test_additions ( int neq, int nvr, LIST *p[] )
{
   int i,fail;
   double *c;
   LIST *l;

   fail = syscon_initialize_number_of_standard_polynomials(neq);
   for(i=0; i<neq; i++)
      for(l=p[i]; l!=NULL; l=l->next)
      {
         double cf[2];
         cf[0] = l->rc;
         cf[1] = l->ic;
         fail = syscon_add_standard_term(i+1,nvr,l->exp,cf);
      }
   printf("\nThe reconstructed system :\n");
   fail = syscon_write_standard_system();
}

LIST *push ( LIST *l, int n, double *c, int *e )
{
   int i;
   LIST *nl = (LIST*)calloc(1,sizeof(LIST));

   nl->rc = c[0];
   nl->ic = c[1];
   nl->exp = (int*)calloc(n,sizeof(int));
   for(i=0; i<n; i++)
      nl->exp[i] = e[i];

   nl->next = l;

   return nl;
}

void write_monomial_list ( LIST *l, int n )
{
   LIST *p;
   int k;

   for(p=l; p!= NULL; p=p->next)
   {
      printf("%22.15e  %22.15e", p->rc, p->ic);
      for(k=0; k<n; k++) printf(" %d", p->exp[k]);
      printf("\n");
   }
}

void test_symbol_table ( int *nvr )
{
   int fail;
   int sbsize = 80;
   int size;
   char *s;

   fail = syscon_number_of_symbols(nvr);
   printf("number of symbols in the table : %d\n", *nvr);
   printf("the symbols :");
   fail = syscon_write_symbols();
   printf("\n");
   size = (*nvr)*sbsize;
   s = (char*)calloc(size, sizeof(char));
   fail = syscon_string_of_symbols(&size, s);
   printf("the string of symbols : %s\n", s);
   printf("number of characters in the string : %d\n", size);
}

void test_standard_container ( void )
{
   int neq,nvr,fail,*d;
   double *c;

   fail = syscon_read_standard_system();
   fail = syscon_write_standard_system();
   fail = syscon_number_of_standard_polynomials(&neq);

   test_symbol_table(&nvr);

   printf("\nThe number of equations in the system : %d\n", neq);
   printf("The number of variables in the system : %d\n\n", nvr);
   {
      LIST *p[neq];
      int i;

      test_retrievals(neq, nvr, p);

      fail = syscon_clear_standard_system();

      for(i=0; i<neq; i++)
      {
         printf("The terms in polynomial %d : \n", i+1);
         write_monomial_list(p[i], nvr);
      }

      test_additions(neq, nvr, p);
   }
}

void test_dobldobl_container ( void )
{
   int fail,n,i,t,deg,nvr;

   fail = syscon_read_dobldobl_system();
   fail = syscon_write_dobldobl_system();
   fail = syscon_number_of_dobldobl_polynomials(&n);
   printf("number of polynomials : %d\n",n);
   printf("number of terms in each polynomial :");
   for(i=1; i<=n; i++)
   {
      fail = syscon_number_of_dobldobl_terms(i,&t);
      printf(" %d",t);
   }
   printf("\ndegree of each polynomial :");
   for(i=1; i<=n; i++)
   {
      fail = syscon_degree_of_dobldobl_polynomial(i,&deg);
      printf(" %d",deg);
   }
   printf("\nthe polynomials as strings :\n");
   for(i=1; i<=n; i++)
   {
      char buffer[2000];
      int nc;
      fail = syscon_load_dobldobl_polynomial(i,&nc,buffer);
      printf("polynomial %d : %s\n",i,buffer);
   }
   printf("\n");
   test_symbol_table(&nvr);
}

void test_quaddobl_container ( void )
{
   int fail,n,i,t,deg,nvr;

   fail = syscon_read_quaddobl_system();
   fail = syscon_write_quaddobl_system();
   fail = syscon_number_of_quaddobl_polynomials(&n);
   printf("number of polynomials : %d\n",n);
   printf("number of terms in each polynomial :");
   for(i=1; i<=n; i++)
   {
      fail = syscon_number_of_quaddobl_terms(i,&t);
      printf(" %d",t);
   }
   printf("\ndegree of each polynomial :");
   for(i=1; i<=n; i++)
   {
      fail = syscon_degree_of_quaddobl_polynomial(i,&deg);
      printf(" %d",deg);
   }
   printf("\nthe polynomials as strings :\n");
   for(i=1; i<=n; i++)
   {
      char buffer[4000];
      int nc;
      fail = syscon_load_quaddobl_polynomial(i,&nc,buffer);
      printf("polynomial %d : %s\n",i,buffer);
   }
   printf("\n");
   test_symbol_table(&nvr);
}

void test_multprec_container ( void )
{
   int deci,fail,n,i,t,deg,nvr;

   printf("\nGive the number of decimal places in the working precision : ");
   scanf("%d",&deci);
   fail = syscon_read_multprec_system(deci);
   fail = syscon_write_multprec_system();
   fail = syscon_number_of_multprec_polynomials(&n);
   printf("number of polynomials : %d\n",n);
   printf("number of terms in each polynomial :");
   for(i=1; i<=n; i++)
   {
      fail = syscon_number_of_multprec_terms(i,&t);
      printf(" %d",t);
   }
   printf("\ndegree of each polynomial :");
   for(i=1; i<=n; i++)
   {
      fail = syscon_degree_of_multprec_polynomial(i,&deg);
      printf(" %d",deg);
   }
   printf("\nthe polynomials as strings :\n");
   for(i=1; i<=n; i++)
   {
      char buffer[4000];
      int nc;
      fail = syscon_load_multprec_polynomial(i,&nc,buffer);
      printf("polynomial %d : %s\n",i,buffer);
   }
   printf("\n");
   test_symbol_table(&nvr);
}

void show_random_system ( void )
{
   int n,m,d,c,neq,fail;

   printf("\nGenerating a random system ...\n");
   printf("-> enter the number of variables : "); scanf("%d", &n);
   printf("-> enter the number of equations : "); scanf("%d", &neq);
   printf("-> enter the number of monomials : "); scanf("%d", &m);
   printf("-> enter the degree bound : "); scanf("%d", &d);
   printf("-> enter the coefficient type : "); scanf("%d", &c);

   fail = syscon_random_system(n,m,d,c,neq);
   if(fail == 0) printf("Making of a random system succeeded.\n");
   syscon_write_standard_system();
}
