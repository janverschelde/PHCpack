/* simple test program in C on the Ada procedure use_syscon */

#include<stdio.h>
#include<stdlib.h>

extern void adainit();
extern int _ada_use_syscon ( int task, int *a, int *b, double *c );
extern void adafinal();

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

void write ( LIST *l, int n );
/* writes the list of monomials */

void test_retrievals ( int n, LIST *p[n] );
/* prints the system in the container,
 * and returns the lists of terms in p */

void test_additions ( int n, LIST *p[n] );
/* fills up the container with the terms in p */

int main(void)
{
   printf("\nTesting Systems Container...\n");

   adainit();
   {
      int n,fail,*d;
      double *c;

      fail = _ada_use_syscon(0,&n,d,c);    /* read system */
      fail = _ada_use_syscon(1,&n,d,c);    /* write system */
      fail = _ada_use_syscon(2,&n,d,c);    /* get dimension */

      printf("\nThe dimension of the system : %d\n",n);
      {
         LIST *p[n];
         int i;

         test_retrievals(n,p);

         fail = _ada_use_syscon(7,&n,d,c); /* clear container */

         for(i=0; i<n; i++)
         {
            printf("The terms in polynomial %d : \n", i+1);
            write(p[i],n);
         }

         test_additions(n,p);
      }
   }
   adafinal();

   return 0;
}

void test_retrievals ( int n, LIST *p[n] )
{
   int i,j,k,fail;
   double c[2];
   int m[3],d[n];

   m[0] = n;
   for(i=1; i<=n; i++)
   {
      m[1] = i;
      fail = _ada_use_syscon(4,m,d,c);
      printf("  #terms in polynomial %d : %d\n", i, m[0]);

      p[i-1] = NULL;
      for(j=1; j<=m[0]; j++)
      {
         m[2] = j;
         fail = _ada_use_syscon(5,m,d,c);
         printf(" %.15f  %.15f", c[0], c[1]);
         for (k=0; k<n; k++) printf(" %d", d[k]);
         printf("\n");
         p[i-1] = push(p[i-1],n,c,d);
      }
   }
}

void test_additions ( int n, LIST *p[n] )
{
   int i,*d,fail,m[3];
   double *c;
   LIST *l;

   fail = _ada_use_syscon(3,&n,d,c);
   m[0] = n;
   for(i=0; i<n; i++)
      for(l=p[i]; l!=NULL; l=l->next)
      {
         double cf[2];
         cf[0] = l->rc;
         cf[1] = l->ic;
         m[1] = i+1;
         fail = _ada_use_syscon(6,m,l->exp,cf);
      }
   fail = _ada_use_syscon(1,m,d,c);
}

LIST *push ( LIST *l, int n, double *c, int *e )
{
   int i;
   LIST *nl = (LIST*)calloc(1,sizeof(LIST));

   nl->rc = c[0];
   nl->ic = c[1];
   nl->exp = calloc(n,sizeof(int));
   for (i=0; i<n; i++)
     nl->exp[i] = e[i];

   nl->next = l;

   return nl;
}

void write ( LIST *l, int n )
{
   LIST *p;
   int k;

   for(p=l; p!= NULL; p=p->next)
   {
      printf(" %.15f  %.15f", p->rc, p->ic);
      for (k=0; k<n; k++) printf(" %d", p->exp[k]);
      printf("\n");
   }
}
