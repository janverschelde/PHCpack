/* This program calls the Ada function Pieri_Root_Count */

#include <stdio.h>

extern int _ada_pieri_root_count(int m, int p, int q);
extern void adainit();
extern void adafinal();

int main()
{
  int m,p,q;
  int root_number;

  printf("Give the dimension of the input planes : "); scanf("%d",&m);
  printf("Give the dimension of the output planes : "); scanf("%d",&p);
  printf("Give the degree of the maps : "); scanf("%d",&q);

  adainit();
  root_number = _ada_pieri_root_count(m,p,q);
  adafinal();

  printf("The number of the roots : %d\n", root_number);

  return 0;
}
