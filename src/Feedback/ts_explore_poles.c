/* This program calls the Ada function Pieri_Solver */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "c2ada_dc_matrix.h"
/* #include "poly_matrix.h" */
#include "dc_matrix.h"
#include "dc_inverse.h"

#include "dcmplx.h"
#include "pieri_sols.h"

extern int _ada_pieri_solver(int m, int p, int q, int nb, int output_level,
                             double *pts, double *input, char *filename);
extern void adainit();
extern void adafinal();

void feedback(int n, int m, int p, int q, int nb, int output_level, int nn, int input, char *ifn, char *ofn);

int main(int argc, char **argv)
{
  int n, m,p,q, nb, level, nn,dimpts,dimpla,i, input;
  double *points,*planes;
  FILE *ifp, *ofp;

  if(argc!=3)
  {
     printf("\n%s\n\n",
	    "Usage: ts_feedback input_file output_file");
     exit(1);
  }
    ifp=fopen(argv[1], "r"); /*open for reading*/
    ofp=fopen(argv[2], "w"); /*open for writing*/

  fscanf(ifp, "%d", &n);
  printf("The number of the internal states for the given plant (A, B, C) n = %d\n", n);
  fscanf(ifp, "%d", &m);  
  printf( "The system's input dimension m = %d.\n", m);
  fscanf(ifp, "%d", &p);
  printf( "The system's output dimension p = %d.\n", p);
  fscanf(ifp, "%d", &q);
  printf( "The number of the internal states for the dynamic compensators q = %d.\n", q);
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
  printf("\nType 1, 2 to select input\n"); 
  printf("Type 1 for the interactive input\n");
  printf("Type 2 for the random input\n");  
  printf("The selected input is: %d\n", input);    

  nn = m*p + q*(m+p);
  srand(time(NULL));
  feedback(n, m, p, q, nb, level, nn, input, argv[1], argv[2]);

  return 0;
}

void feedback(int n, int m, int p, int q, int nb, int output_level, int nn, int input, char *ifn, char *ofn)
{
  dcmplx A[n][n], B[n][m], C[p][n];
  dcmplx s[nn], Is[n][n], Is_A[n][n], tmp[p][n], M[p][m];  
  double a[nn*2], b[nn*p*m*2], *points, *planes;
  int i, j, k, start, nbsols;
  FILE *ifp, *ofp;
  double pole[4], grid;
int dimpts, dimpla;
  

 ofp=fopen(ofn, "w"); /*open for writing*/
 if(input == 1)
 {
  ifp=fopen(ifn, "r"); /*open for reading*/
  skip(ifp);

  read_dcmatrix2(n, n, A, ifp);
  /* fprintf(ofp, "The given matrix A(%d*%d) is:\n", n, n);
     print_dcmatrix1(n, n, A, ofp); */

  read_dcmatrix2(n, m, B, ifp);
  /* fprintf(ofp, "The given matrix B(%d*%d) is:\n", n, m);
     print_dcmatrix1(n, m, B, ofp); */

  read_dcmatrix2(p, n, C, ifp);
  /* fprintf(ofp, "The given matrix C(%d*%d) is:\n", p, n);
     print_dcmatrix1(p, n, C, ofp);  */
  
  for(i=0; i<n+q; i++)
    read_dcmplx1(&s[i], ifp);
  for(i=n+q; i<nn; i++)     /*generate more poles as interpolation points */
  {
    s[i] = create1(cos(rand()));
    if(s[i].re>0)
      s[i] = min_dcmplx(s[i]);  
  } 
  /* fseek(ifp, 0, SEEK_SET); */
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

  s[0] = create2(-0.23423423423, 0);  /* fix one pole for test purpose */
  for(i=1; i<nn; i++)
  {
    s[i] = random_dcmplx1();
  }
  printf("\nThe given poles are:\n");
  for(i=0; i<nn; i++)
  {
    writeln_dcmplx(s[i]);
  }
 }

  adainit();
 /* Explore all the possible poles */
 grid = 1; 
for(pole[0]=-50.0; pole[0]<0; pole[0]=pole[0]+grid)
 for(pole[1]=pole[0]+grid; pole[1]<0; pole[1]=pole[1]+grid)
  for(pole[2]=pole[1]+grid; pole[2]<0; pole[2]=pole[2]+grid)
   for(pole[3]=pole[2]+grid; pole[3]<0; pole[3]=pole[3]+grid)
{
  printf("The poles are: ");
  for(i=0; i<4; i++)
  {
    printf("%.1lf ", pole[i]);
    s[i].re = pole[i];
    s[i].im = 0.0;  
  }
  printf("\n");
  j = 0;
  for(i=0; i<nn; i++)
  {
    a[j++] = s[i].re;
    a[j++] = s[i].im;
  }

  start = 0;
  for(k=0; k<n+q; k++)
  /* for(k=0; k<6; k++) */
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
    
   /*  printf("The M%d matrix is:\n", k);
       print_dcmatrix(p, m, M); */
  }

  /* generate some random planes */
  for( i=start ; i<nn*p*m*2; i++ )
  {
    b[i++] = cos(rand()); 
    b[i] = 0.0;
  }
 
  fflush(stdout);
  fflush(ofp); 

  /*  printf("\nComputing the feedback law with PHC ...\n"); */

  nbsols = _ada_pieri_solver(m, p, q, nb, output_level, a, b, ofn);


  /*  printf("See %s for the output feedback result.\n", ofn); 
      printf("The number of solutions : %d\n", nbsols);   */ 
  /*  fclose(ofp); */
} 
  adafinal(); 
  /*  
  printf("the double array for s is:\n");
  for(i=0; i<nn*2; i++)
  {
     printf("%.15le ", a[i]);
     if((i+1)%4==0)
        printf("\n");
  }
  
  printf("the double array for M is:\n");
   for(i=0; i<nn*p*m*2; i++)
   {
      printf("%.15le ", b[i]);
      if((i+1)%4==0)
        printf("\n");
   }
  */
}






