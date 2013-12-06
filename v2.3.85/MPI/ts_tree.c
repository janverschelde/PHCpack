#include <stdio.h>
#include "tree.h"

int main()
{
  int m, p, q, sum;

  printf("Give p, dimension of the solution planes :");
  scanf("%d", &p);
  printf("Give m, the co-dimension so that n = m+p :");
  scanf("%d", &m);
  printf("Give q, the degree of the maps :");
  scanf("%d", &q);

  if(q==0)
    sum = Build_Bottom_Tree(m,p);
  else
    sum = Q_Build_Bottom_Tree(m,p,q); 
  printf("The number of roots : %d\n", sum);

  return 0;
}
