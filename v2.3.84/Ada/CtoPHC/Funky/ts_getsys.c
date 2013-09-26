#include<stdio.h>
#include<stdlib.h>

extern void adainit();
extern void _ada_getsys1 ( void );
extern void adafinal();

int main(void)
{
   printf("\nTesting Reading of Polynomial Systems\n");

   adainit();
   _ada_getsys1();
   adafinal();

   return 0;
}
