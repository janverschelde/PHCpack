/* file realization.h contains the prototype for functions used for  realizations */
#include "poly_dcmplx.h"

void limit(int n, int m, POLY dd, POLY MM[n][m], dcmplx D[n][m] );
/* get the limit of the transfer function when s go to infinity 
   if the result is not finite, there will be no realization. The 
   result will be stored in D matrix, sum is the order of the realization */

int get_deg(int n, POLY M1[n][n], int a[n], int * sum);
/* get the column degree of the matrix. if a column degree is 0, returns 0; otherwise returns 1 */


void factorize_M1(int p, POLY M1[p][p], int a[p], int sum, dcmplx DH[p][p],
         POLY Dia[p][p], dcmplx DL[p][sum], POLY S[sum][p] );
/* factorize M1(s) = DH*Dia(s) + DL*S(s) */

POLY Get_Hs (int p, int m,  POLY M1[p][p], POLY M2[m][p], POLY H[m][p]);
/* get Hs function and return the least common denominator of Hs matrix, the
   M1 on return will be the inverse matrix (without ds) of the given M1 */ 

void Get_Bc ( int sum, int p, int d[p], dcmplx DH[p][p], dcmplx Bc[sum][p]);
/* given DH and the degree information, return Bc */

void Get_Bc1(int sum, int new_p, int p, int b[new_p],dcmplx DH[new_p][new_p], dcmplx Bc[sum][p]);
/* special case when the column degree of transfer function is zero */

void Get_Ac ( int sum, int p, int d[p], dcmplx DH[p][p], dcmplx DL[p][sum], dcmplx Ac[sum][sum]);
/* given DH, DL and the degree information, return Ac */

void Get_Cc ( int sum, int m, int p, POLY N[m][p], POLY D[p][p], int d[p], dcmplx Dc[m][p], dcmplx Cc[m][sum]);
/* given N, D, Dc and the degree information, return Cc */

void Get_Cc1 ( int sum, int m, int p, POLY N[m][p], int d[p], dcmplx Cc[m][sum]);
/* special case when the column degree of transfer function is zero */

void realization(int p, int m, int q, POLY M1[p][p], POLY M2[m][p], dcmplx Ac[q][q], dcmplx Bc[q][p],
                 dcmplx Cc[m][q], dcmplx Dc[m][p]);

void realization_report(int p, int m, int q, POLY M1[p][p], POLY M2[m][p], dcmplx Ac[q][q], dcmplx Bc[q][p],
			dcmplx Cc[m][q], dcmplx Dc[m][p]);
/* for test purpose only, it will print out the intermediate output.*/

POLY** eliminate_columns(int m, int p, POLY H_s[m][p], int a[p], int new_p);
/* eliminate the columns whose degree is 0 and returns a new transfer function */

void factorize_H1(int m, int p, POLY ds, POLY H[m][p], POLY N[m][p], POLY D[p][p]);
/* call factorize_H function with POLY H[m][n] instead of POLY **H */

void factorize_H(int m, int p, POLY ds, POLY** H, POLY N[m][p], POLY D[p][p]);
/* factorizes transfer function H = N * D_inverse */

void add_zero_columns(int q, int p, int a[p], dcmplx Bc[q][p]);
/* add zero columns to Bc matrix where the column degree of H is 0 */ 

void right_coprime(int n1, int n2, int m, POLY D[n1][m], POLY N[n2][m], POLY Gr[m][m]);
/* use hermite form to make N and D matrices right coprime, and save the coprime matrices in N and D */
