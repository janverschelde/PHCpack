/* This program calls the Ada function Pieri_Solver */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern int _ada_pieri_solver(int m, int p, int q, int nb, int output_level,
                             double *pts, double *input, char *filename);
extern void adainit();
extern void adafinal();

int main()
{
  int m,p,q,nb,level,n,dimpts,dimpla,i,nbsols;
  double *points,*planes;
  char filename[80];

  printf("\nTesting the call to Pieri homotopies...\n\n");

  printf("Give the name of the output file : ");
  scanf("%s", filename);
  printf("The output file name is %s.\n", filename);

  printf("Give the dimension of the input planes : "); scanf("%d",&m);
  printf("Give the dimension of the output planes : "); scanf("%d",&p);
  printf("Give the degree of the maps : "); scanf("%d",&q);
  printf("Give the number of maps wanted (<0 for all) : "); scanf("%d",&nb);

  printf("Give the amount of intermediate output : \n");
  printf("  0. no intermediate output;\n");
  printf("  1. only final determinant validation;\n");
  printf("  2. + validation of all intermediate determinants;\n");
  printf("  3. + intermediate output of all path trackers.\n");
  printf("Type 0, 1, 2, or 3 to select output level : ");
  scanf("%d", &level);

  n = m*p + q*(m+p);

  dimpts = 2*n;
  points = (double*) calloc(dimpts,sizeof(double));
  for(i=0; i<dimpts; i++)
     points[i] = cos(rand());

  dimpla = 2*n*m*p;
  planes = (double*) calloc(dimpla,sizeof(double)); 
  for(i=0; i<dimpla; i++)
     planes[i] = cos(rand());

  adainit();
  nbsols = _ada_pieri_solver(m,p,q,nb,level,points,planes,filename);
  adafinal();

  printf("The number of solutions : %d\n", nbsols);

  return 0;
}
