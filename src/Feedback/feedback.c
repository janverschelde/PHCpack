/* This program calls the Ada function pieri_solver */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <string.h>

#include "timer.h"
#include "c2ada_dc_matrix.h"
#include "dc_matrix.h"
#include "dc_inverse.h"

#include "dcmplx.h"

#include "pieri_sols.h"

extern void adainit();
extern int _ada_pieri_solver
 ( int m, int p, int q, int nb, int output_level,
   double *pts, double *input, char *filename );
extern void adafinal();

void skip_info(FILE *ifp);
/* skips the dimensions information */

void feedback
 ( int n, int m, int p, int q, int nb, int output_level, 
   int nn, int input, char *ifn, char *ofn );

int main_feedback ( char *inputfile, char *outputfile )
{
   int n, m,p,q, nb, level, nn, i, input;
   double *points,*planes;
   FILE *ifp, *ofp;

   printf("the output put file name : %s\n", outputfile);

   ifp=fopen(inputfile, "r");  /* open for reading */
   if(ifp == NULL)
   {
      printf("Opening %s for reading failed!\n",inputfile);
      exit(1);
   }
   fscanf(ifp, "%d", &n);
   printf("The number of the internal states n = %d\n", n);
   fscanf(ifp, "%d", &m);  
   printf("The system's input dimension m = %d.\n", m);
   fscanf(ifp, "%d", &p);
   printf("The system's output dimension p = %d.\n", p);
   fscanf(ifp, "%d", &q);
   printf("The number of internal states for compensators q = %d.\n", q);
   fscanf(ifp, "%d", &nb);
   printf("Give the number of maps wanted (<0 for all) : %d ", nb);

   fscanf(ifp, "%d", &level);
   printf( "\nType 0, 1, 2, or 3 to select output level :\n");
   printf( "  0. no intermediate output;\n");
   printf( "  1. only final determinant validation;\n");
   printf( "  2. + validation of all intermediate determinants;\n");
   printf( "  3. + intermediate output of all path trackers.\n");
   printf( "The amount of the intermediate output: %d\n", level);

   fscanf(ifp, "%d", &input);
   printf("Type 0, 1, 2 to select input :\n"); 
   printf("  0. interactive input of real numbers\n");
   printf("  1. interactive input of complex numbers\n");
   printf("  2. random input of complex numbers\n");  
   printf("  3. random input of real numbers\n");
   printf("  4. interactive input of matrices and poles\n");
   printf("The selected input is: %d\n", input);    

   nn = m*p + q*(m+p);
   srand(time(NULL));

   fclose(ifp); 

   feedback(n, m, p, q, nb, level, nn, input, inputfile, outputfile);

   return 0;
}

void skip_info(FILE *ifp)
{
  int i;
  char s[30];
  
  i=0;
  while(1)            /*skips all the space before the first number*/
  {
    if(!isspace(fgetc(ifp))) 
       break;
  } 
  fgets(s, 30, ifp);  /*skips first line */
}

void feedback
 ( int n, int m, int p, int q, int nb, int output_level,
   int nn, int input, char *ifn, char *ofn)
{
   dcmplx A[n][n], B[n][m], C[p][n];
   dcmplx s[nn], Is[n][n], Is_A[n][n], tmp[p][n], M[p][m];  
   double a[nn*2], b[nn*p*m*2], *points, *planes, r;
   int i, j, k, start, nbsols;
   FILE *ifp, *ofp;
   timer t_phc;

   ofp=fopen(ofn, "w"); /*open for writing*/
   fprintf(ofp, "n=%d\n", n);
   fprintf(ofp, "m=%d\n", m);
   fprintf(ofp, "p=%d\n", p);
   fprintf(ofp, "q=%d\n", q);
 
   if(input == 0)
   {
      ifp=fopen(ifn, "r"); /* open for reading */
      skip_info(ifp);
      read_dcmatrix0(n, n, A, ifp);
      printf("The given matrix A(%d*%d) is:\n", n, n);
      print_dcmatrix(n, n, A);
      read_dcmatrix0(n, m, B, ifp);
      printf("The given matrix B(%d*%d) is:\n", n, m);
      print_dcmatrix(n, m, B);
      read_dcmatrix0(p, n, C, ifp);
      printf("The given matrix C(%d*%d) is:\n", p, n);
      print_dcmatrix(p, n, C); 
      for(i=0; i<n+q; i++)
         read_dcmplx0(&s[i], ifp);
      for(i=n+q; i<nn; i++)  /* generate more poles as interpolation points */
      {
        s[i] = create1(cos(rand()));
        if(s[i].re>0) s[i] = min_dcmplx(s[i]);  
      } 
      fclose(ifp);
   }
   if(input == 1)
   {
      ifp=fopen(ifn, "r"); /*open for reading*/
      skip_info(ifp);
      read_dcmatrix2(n, n, A, ifp);
      printf("The given matrix A(%d*%d) is:\n", n, n);
      print_dcmatrix(n, n, A);
      read_dcmatrix2(n, m, B, ifp);
      printf("The given matrix B(%d*%d) is:\n", n, m);
      print_dcmatrix(n, m, B);
      read_dcmatrix2(p, n, C, ifp);
      printf("The given matrix C(%d*%d) is:\n", p, n);
      print_dcmatrix(p, n, C); 
      for(i=0; i<n+q; i++)
         read_dcmplx1(&s[i], ifp);
      for(i=n+q; i<nn; i++)  /* generate more poles as interpolation points */
      {
         s[i] = create1(cos(rand()));
         if(s[i].re>0) s[i] = min_dcmplx(s[i]);  
      } 
      fclose(ifp);
   }
   if(input==2)
   {
      random_dcmatrix ( n, n, A);
      printf("\nThe random generated matrix A is:\n");
      print_dcmatrix(n, n, A);
      random_dcmatrix ( n, m, B);
      printf("\nThe random generated matrix B is:\n");
      print_dcmatrix(n, m, B);
      random_dcmatrix ( p, n, C);
      printf("\nThe random generated matrix C is:\n");
      print_dcmatrix(p, n, C);
      s[0] = create2(-0.23423423423, 0); /* fix one pole to test realization */
      for(i=1; i<nn; i++)
      {
         r = rand();
         s[i] = create2(cos(r), sin(r));
         if(s[i].re>0) s[i] = min_dcmplx(s[i]); 
         s[++i] = conjugate(s[i]);
         if(i==(nn-2))
         {
            if((nn%2)==0) s[++i] = create1(-1.0); 
         }  
      }
      printf("\nThe random generated poles are:\n");
      for(i=0; i<nn; i++) writeln_dcmplx(s[i]);
   }
   if(input==3)
   {
      random_dcmatrix0 ( n, n, A);
      printf("\nThe random generated matrix A is:\n");
      print_dcmatrix(n, n, A);
      random_dcmatrix0 ( n, m, B);
      printf("\nThe random generated matrix B is:\n");
      print_dcmatrix(n, m, B);
      random_dcmatrix0 ( p, n, C);
      printf("\nThe random generated matrix C is:\n");
      print_dcmatrix(p, n, C);
      s[0] = create2(-0.23423423423, 0);  /* fix one pole for test */
      for(i=1; i<nn; i++)
      {
         r = rand();
         s[i] = create2(cos(r), sin(r));
         if(s[i].re>0) s[i] = min_dcmplx(s[i]);
         s[++i] = conjugate(s[i]);
         if(i==(nn-2))
         {
            if((nn%2)==0) s[++i] = create1(-1.0);
         }
      }
      printf("\nThe random generated poles are:\n");
      for(i=0; i<nn; i++) writeln_dcmplx(s[i]);
   }
   if(input == 4)
   {
      ifp=fopen(ifn, "r"); /*open for reading*/
      skip_info(ifp);
      read_dcmatrix0(n, n, A, ifp);
      printf( "The given matrix A(%d*%d) is:\n", n, n);
      print_dcmatrix(n, n, A);

      read_dcmatrix0(n, m, B, ifp);
      printf("The given matrix B(%d*%d) is:\n", n, m);
      print_dcmatrix(n, m, B);

      read_dcmatrix0(p, n, C, ifp);
      printf("The given matrix C(%d*%d) is:\n", p, n);
      print_dcmatrix(p, n, C);

      for(i=0; i<n+q; i++) read_dcmplx1(&s[i], ifp);
      for(i=n+q; i<nn; i++)  /*generate more poles as interpolation points */
      {
         s[i] = create1(cos(rand()));
         if(s[i].re>0) s[i] = min_dcmplx(s[i]);
      }
      fclose(ifp);
   }
   /* print the input matrices in matlab format for further study */
   fprintf(ofp,"A=[\n");  print_dcmatrix1(n, n, A, ofp); fprintf(ofp,"]\n");
   fprintf(ofp,"B=[\n");  print_dcmatrix1(n, m, B, ofp); fprintf(ofp,"]\n");
   fprintf(ofp,"C=[\n");  print_dcmatrix1(p, n, C, ofp); fprintf(ofp,"]\n");

   fprintf(ofp, "\nPoles=[");
   for(i=0; i<nn; i++)
   {
       write_dcmplx1(s[i], ofp);
       if(i!=(nn-1)) fprintf(ofp, ",");
   }
   fprintf(ofp, "]\n");
   /* end of input */
   j = 0;
   for(i=0; i<nn; i++)
   {
      a[j++] = s[i].re;
      a[j++] = s[i].im;
   }
   start = 0;
   for(k=0; k<n+q; k++)
   { 
      for(i=0; i<n; i++)
         for(j=0; j<n; j++)
         {
           if(i==j) Is[i][j] = s[k];
           else Is[i][j] = zero;
         }
      sub_dcmatrix(n, n, Is, A, Is_A);
      dcinverse(n, Is_A);

      multiply_dcmatrix(p, n, n, C, Is_A, tmp);
      multiply_dcmatrix(p, n, m, tmp, B, M);
      c2ada_dc_matrix( p, m, M, nn*p*m*2, b, start);
      start = start + p*m*2;  
   }
   /* generate some random planes */
   for( i=start ; i<nn*p*m*2; i++ )
   {
      b[i++] = cos(rand()); 
      b[i] = 0.0;
   }
   fflush(stdout);
   fflush(ofp); 

   printf("\nComputing the feedback laws with Pieri homotopies ...\n");
   tstart(&t_phc);
   /* adainit(); */ /* not needed as main program is in Ada */
   nbsols = _ada_pieri_solver(m, p, q, nb, output_level, a, b, ofn);
   /* adafinal(); */
   tstop(&t_phc);
   printf("Total "); /* This subroutine spends almost all the time */
   tprint(t_phc);

   printf("\nSee %s for the realization of the output feedbacks.\n", ofn);
   /*  fclose(ofp); */
}
