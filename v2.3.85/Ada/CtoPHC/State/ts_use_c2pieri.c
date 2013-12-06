/* simple test program to the interactive operations of use_c2pieri */

#include<stdio.h>
#include<stdlib.h>
#include<time.h>

extern void adainit();
extern int _ada_use_c2pieri ( int job, int *a, int *b, double *c );
extern void adafinal();

int initialize_dimensions ( int *m, int *p, int *q );
/* prompts user to give (m,p,q) and initializes the machine,
 * returns 0 and the values for (m,p,q) if all went well,
 * returns 1 if a failure happened */

int initialize_input_planes ( int m, int p, int q );
/* initializes the machine with m*p + q*(m+p) random input m-planes,
 * returns 0 if all went well, returns 1 if failed */

int initialize_interpolation_points ( int m, int p, int q );
/* initializes the machine with m*p + q*(m+p) random interpolation points,
 * returns 0 if all went well, returns 1 if failed */

int store_start_pivots ( int p );
/* prompts the user for top and bottom pivots of a localization pattern
 * for a curve at the start and stores those pivots in the machine,
 * returns 0 if all went well, returns 1 if failed */

int store_target_pivots ( int p );
/* prompts the user for top and bottom pivots of a localization pattern
 * for a target curve and stores those pivots in the machine,
 * returns 0 if all went well, returns 1 if failed */

int store_start_curve ( void );
/* stores the coefficients of the start solution curve,
 * returns 0 if all went well, returns 1 if failed */

int retrieve_target_curve ( void );
/* retrieves the coefficients of the target solution curve,
 * returns 0 if all went well, returns 1 if failed */

int main(void)
{
   printf("\nTesting Pieri homotopies from C ...\n");

   srand(time(NULL));

   adainit();
   {
      int m,p,q,fail,choice = 0;
      int *a,*b;
      double res,*c;
      char skip_newline;

      m = 0; q = 0; p = 0;
      do
      { 
         fail = _ada_use_c2pieri(choice,a,b,c);
         printf("Give a number (0 for menu, -1 to exit) : ");
         scanf("%d",&choice);

         switch (choice)
	 {
           case -1: { printf("bye bye\n");
                      break; }
           case  0:   break;
           case  1: { fail = initialize_dimensions(&m,&p,&q);
                      printf("C: m = %d, p = %d, q = %d\n", m,p,q);
                      break; }
           case  2: { fail = initialize_input_planes(m,p,q);
                      break; }
           case  3: { fail = initialize_interpolation_points(m,p,q);
                      break; }
           case  4 : { fail = store_start_pivots(p); break; }
           case  5 : { fail = store_target_pivots(p); break; }
           case  6 : { fail = store_start_curve(); break; }
           case  7 : { fail = retrieve_target_curve(); break; }
           case  8 : { fail = _ada_use_c2pieri(8,a,b,c); break; }
           case  9 : { fail = _ada_use_c2pieri(9,a,b,c); break; }
           case 10 : { fail = _ada_use_c2pieri(10,a,b,&res);
                       printf("Sum of residuals : %.15le\n",res);
                       break; }
           case 11 : { fail = _ada_use_c2pieri(11,a,b,&res);
                       printf("Sum of residuals : %.15le\n",res);
                       break; }
           case 12 : { fail = _ada_use_c2pieri(12,a,b,c); break; }
           default:   fail = -1;
         }

         if (choice < 0) break;

         if (fail == -1)
            printf("Invalid choice.  Please try again...\n");
	 else if (fail == 0)
            printf("operation succeeded\n");
         else
            printf("OPERATION FAILED!!!\n");

         skip_newline = getchar(); choice = 0;

      } while (choice >= 0);
   }
   adafinal();

   return 0;
}

int initialize_dimensions ( int *m, int *p, int *q )
{
   int fail,mpq[3];
   int *b;
   double *c;

   printf("\nReading the dimensions of the problem...\n");
   printf("  Give m (dimension of input planes)  : "); scanf("%d",&mpq[0]);
   printf("  Give p (dimension of output planes) : "); scanf("%d",&mpq[1]);
   printf("  Give q (degree of solution curves)  : "); scanf("%d",&mpq[2]);

   fail = _ada_use_c2pieri(1,mpq,b,c);
   *m = mpq[0];
   *p = mpq[1];
   *q = mpq[2];

   return fail;
}

int initialize_input_planes ( int m, int p, int q )
{
   int n = m*p + q*(m+p);
   int nbcff,i;
   int fail,mpq[3];
   double *c;

   if (n > 0)
   {
      mpq[0] = m; mpq[1] = p; mpq[2] = q;
      nbcff = 2*(m+p)*m*n;    
      printf("Generating %d random doubles... \n",nbcff);
      c = (double*)calloc(nbcff,sizeof(double));
      for (i=0; i<nbcff;i++) c[i] = ((double) rand())/RAND_MAX;
      fail = _ada_use_c2pieri(2,mpq,&n,c);
   }
   else
   {
      printf("m = %d, p = %d, q = %d => invalid n = %d\n", m,p,q,n);
      printf("Please initialize first (m,p,q).\n");
      fail = 1;
   }

   return fail;
}

int initialize_interpolation_points ( int m, int p, int q )
{
   int n = m*p + q*(m+p);
   int nbcff,i;
   int fail,mpq[3];
   double *c;

   if (n > 0)
   {
      mpq[0] = m; mpq[1] = p; mpq[2] = q;
      nbcff = 2*n;
      printf("Generating %d random doubles... \n",nbcff);
      c = (double*)calloc(nbcff,sizeof(double));
      for (i=0; i<nbcff;i++) c[i] = ((double) rand())/RAND_MAX;
      fail = _ada_use_c2pieri(3,mpq,&n,c);
   }
   else
   {
      printf("m = %d, p = %d, q = %d => invalid n = %d\n", m,p,q,n);
      printf("Please initialize (m,p,q) first.\n");
      fail = 1;
   }

   return fail;
}

void read_pivots ( int p, int pv[] )
{
   int i;

   printf("  give %d top pivots : ", p);
   for (i=0; i<p; i++) scanf("%d",&pv[i]);

   printf("  give %d bottom pivots : ", p);
   for (i=p; i<2*p; i++) scanf("%d",&pv[i]);
}

int store_start_pivots ( int p )
{
   int fail;
   double *c;

   if (p > 0)
   {
      int pivots[2*p];

      printf("\nReading pivots of start curve...\n");
      read_pivots(p,pivots);
      fail = _ada_use_c2pieri(4,&p,pivots,c); 
   }
   else
   {
      printf("p = %d, please initialize p first.\n", p);
      fail = 1;
   }

   return fail;
}

int store_target_pivots ( int p )
{
   int fail;
   double *c;

   if (p > 0)
   {
      int pivots[2*p];

      printf("\nReading pivots of target curve...\n");
      read_pivots(p,pivots);
      fail = _ada_use_c2pieri(5,&p,pivots,c);
   }
   else
   {
      printf("p = %d, please initialize p first.\n", p);
      fail = 1;
   }

   return fail;
}

int store_start_curve ( void )
{
   int fail,n;
   int *b;

   printf("\nGive degree of freedom at start curve : ");
   scanf("%d", &n);

   printf("Generating %d complex numbers...\n",n);
   {
      double c[2*n];
      int i;

      for (i=0; i<2*n; i++) c[i] = (double) i;

      fail = _ada_use_c2pieri(6,&n,b,c);
   }

   return fail;
}

int retrieve_target_curve ( void )
{
   int fail,n;
   int *b;

   printf("\nGive degree of freedom of target curve : ");
   scanf("%d", &n);

   {
      double c[2*n];
      int i;

      fail = _ada_use_c2pieri(7,&n,b,c);

      printf("the coefficients of the target curve :\n");
      for (i=0; i<n; i = i+2)
         printf("  %.15le  %.15le \n", c[i], c[i+1]);
   }

   return fail;
}
