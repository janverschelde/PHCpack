/* file "pieri_sols.c" contains definitions of prototypes in pieri_sols.h" */
#include <stdlib.h>
#include <ctype.h>
#include "pieri_sols.h"
#include "poly_matrix.h"
#include "poly_dcmplx.h"
#include "realization.h"
#include "dc_matrix.h"
#include "dc_inverse.h"
#include "timer.h"
#include "append_dcmatrix.h"

/* skips the dimensions information */
void skip(FILE *ifp)
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

void ada2c_poly_matrix( int n, int m, POLY c[n][m], int l, int start_a, int a[l], int q, int start_b, double b[q] )
{

  int i, j, k, r, s;
  
  k=start_a;
  s=start_b;
  for(i=0; i<n; i++)
     for(j=0; j<m; j++)
     {
       c[i][j].d = a[k] ;
       c[i][j].p = (dcmplx*) calloc(c[i][j].d+1, sizeof(dcmplx));
       for(r=0; r<=c[i][j].d; r++)
       {
         c[i][j].p[r].re=b[s++];
         c[i][j].p[r].im=b[s++];
       }
       k++;
     }

}

int pieri_sols( int m, int p, int q, int nbsols, int da, int a[da], int db, double b[db],
                 int npt, double pt[npt], int npl, double pl[npl], char *filename )
{
  POLY M1[p][p], M2[m][p];
  int nn = m*p + q*(m+p);
  dcmplx  F[q][q], G[q][p], H[m][q], K[m][p], tmp1[p+m][p], tmp2[p+m][m], whole[p+m][p+m];
  dcmplx  s[nn], M[nn][p][m], dc_M1[p][p], dc_M2[m][p], Im[m][m], det, det1;
  int i, j, k, t, length_a1, length_b1, c_a, c_b, length_suba, length_subb;
  int sub_a[da/nbsols];
  double sub_b[db/nbsols];
  FILE *ifp, *ofp;  
  dcmplx Iq[q][q], Iq_F[q][q], Ip[p][p], t_tmp[m][q],t_tmp1[m][p], t_tmp2[m][p];
  timer t_realization;
  int n, counter=0;
  
  tstart(&t_realization); 
  k=0;
  for(i=0; i<nn; i++)
  {
    s[i].re = pt[k++];
    s[i].im = pt[k++];
  }
  
  t=0;
  for(k=0; k<nn; k++)
  { 
    for(i=0; i<p; i++)
      for(j=0; j<m; j++)
      {
	M[k][i][j].re = pl[t++];
        M[k][i][j].im = pl[t++];
      }
  }

  length_a1 = p*p;
  length_b1 = 0;
  c_a=0;
  c_b=0;
  length_suba = da/nbsols;
  length_subb = db/nbsols;

  for(i=0; i<length_a1; i++ )
     length_b1 = length_b1 +a[i]+1;
  length_b1 = length_b1*2; 

  ofp=fopen(filename, "a"); /*open for writing*/ 

  I_dcmatrix(m, Im);
  for(i=0; i<nbsols; i++)
  {
    for(j=0; j<length_suba; j++)
    {  
      sub_a[j] = a[c_a];
      c_a++;
    }
    for(k=0; k<length_subb; k++)
    {  
      sub_b[k] = b[c_b];
      c_b++;
    }
    ada2c_poly_matrix(p, p, M1, length_suba, 0, sub_a, length_subb, 0, sub_b);
    ada2c_poly_matrix(m, p, M2, length_suba, length_a1, sub_a, length_subb, length_b1, sub_b); 
    fprintf(ofp,"\n%cSolution No. %d\n",37, i+1);
 
    printf("Solution No. %d\n", i);    
    printf("M1=\n");
    print1(p, p, M1);
    printf("M2=\n");
    print1(m, p, M2);
   
    /* test the results with all the input plane  
       for(t=0; t<nn; t++)  */

    /* test the first input plane only */
    t=0;
    {
      evaluate_matrix(p, p, M1, dc_M1, s[t]);
      det = determinant(p, dc_M1);
      /* printf("The determinant of M1 (at s[0]) is: "); */
      /* writeln_dcmplx(det); */

      evaluate_matrix(m, p, M2, dc_M2, s[t]);
      v_append(p, m, p, dc_M1, dc_M2, tmp1);
      v_append(p, m, m, M[t], Im, tmp2);
      h_append(m+p, p, m, tmp1, tmp2, whole);

      det = determinant(m+p, whole);
     /* printf("The determinant of 5.9 with phc output (before multiply M1_inverse) at s[0] is:\n");
        writeln_dcmplx(det); */
    }
  
    realization( p, m, q, M1, M2, F, G, H, K);
  
    /*test the result with the first input plane and evaluate the matrix at s[0] */
    for(k=0; k<q; k++)
      for(j=0; j<q; j++)
      {
        if(k==j) Iq[k][j] = s[0];
        else Iq[k][j] = zero;
      }
    sub_dcmatrix(q, q, Iq, F, Iq_F);
    dcinverse(q, Iq_F);
    multiply_dcmatrix(m, q, q, H, Iq_F, t_tmp);
    multiply_dcmatrix(m, q, p, t_tmp, G, t_tmp1);
    add_dcmatrix(m, p, t_tmp1, K, t_tmp2);
   /* printf("The calculated H matrix after realization at point s[0] is:\n");
      print_dcmatrix(p, m, t_tmp2); */

   printf("The determinant of the characteristic equation of the closed-loop system is:\n");
   I_dcmatrix(p, Ip);
   v_append(p, m, p, Ip, t_tmp2, tmp1);
   v_append(p, m, m, M[0], Im, tmp2);
   h_append(m+p, p, m, tmp1, tmp2, whole);
   det = determinant(m+p, whole);
   printf("det=");
   writeln_dcmplx(det);
   printf("\n");
    
    /* count how many solutions are complex */
    if(dabs(K[0][0].im) > 1e-8) counter++;    

    fprintf(ofp,"F{%d}=[\n", i+1);
    print_dcmatrix1(q, q, F, ofp);
    fprintf(ofp,"];\n");
  
    fprintf(ofp,"G{%d}=[\n", i+1);
    print_dcmatrix1(q, p, G, ofp);
    fprintf(ofp,"];\n");

    fprintf(ofp,"H{%d}=[\n", i+1);
    print_dcmatrix1(m, q, H, ofp);
    fprintf(ofp,"];\n");

    fprintf(ofp,"K{%d}=[\n", i+1);
    print_dcmatrix1(m, p, K, ofp);
    fprintf(ofp,"];\n");
    
    free_matrix(p, p, M1);
    free_matrix(m, p, M2);
  }
  tstop(&t_realization);
  /* fflush(ofp); */
  fclose(ofp);

  printf("The number of real solutions is: %d\n", nbsols-counter );
  printf("The number of complex solutions is: %d\n", counter);
  printf("The number of total solutions is: %d\n", nbsols);
  printf("\nRealization and verification ");
  tprint(t_realization);

  return 0;
}

