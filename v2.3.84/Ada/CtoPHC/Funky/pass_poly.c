/* This program passes a polynomial in a string to phc for printing. */

#include <stdio.h>

extern void _ada_print_poly(int n, char *s);
extern void adainit();
extern void adafinal();

int main()
{
   int i,n;
   char ch,s[80];

   printf("Give number of variables : ");
   scanf("%d", &n);
   printf("Give a polynomial in %d variables, terminate with ';' :\n", n);
   i = 0;
   ch = ' ';
   while (scanf("%c", &ch)==1)
   {
      s[i++] = ch;
      if (ch == ';') break;
   }
   s[i] = '\0';
   printf("The polynomial read : %s\n", s);

   printf("Calling Ada... \n");
   adainit();
   _ada_print_poly(n,s);
   adafinal();
   printf("... done with the call.\n");

   return 0;
}
