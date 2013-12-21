/* This program passes a polynomial in a string to phc for factoring;
   the factors return in a string, separated by semicolons. */

#include <stdio.h>

extern char *_ada_phc_factor(int n, char *s);
extern void adainit();
extern void adafinal();

int main()
{
   int n,i=0;
   char ch,s[80],*f;

   printf("Give number of variables : ");
   scanf("%d", &n);
   printf("Give a polynomial in %d variables, terminate with ';' :\n", n);
   ch = ' ';
   while (scanf("%c", &ch)==1)
   {
      s[i++] = ch;
      if (ch == ';') break;
   }
   s[i] = '\0';
   printf("The polynomial read : %s\n", s);

   printf("Calling PHC... \n");
   adainit();
   f = _ada_phc_factor(n,s);
   adafinal();
   printf("... done with the call.\n");

   printf("The factors are %s \n", f);

   return 0;
}
