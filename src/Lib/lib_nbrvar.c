/* simple test program in C to scan strings for the number of symbols */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phcpack.h"

/*  Prompts the user for one string with all polynomials in a system.
 *  Each polynomial is terminated by one semicolon.
 *  Then scans those strings for the number of symbols used as
 *  variables in the string representation of the system. */

#define  MAXLEN  80    /* maximum #characters in string */

int read_line ( char *s );
/*
 * DESCRIPTION :
 *   Reads one line of characters, which are stored in s.
 *   On return is the number of characters. */

int main(void)
{
   char pols[MAXLEN];
   int nb,nbc,k,dim,fail;

   adainit();
 
   printf("Give one string with all polyomials in the system,\n");
   printf("each polynomial is terminated by a semi colon;\n");
   printf("Enter your system : ");

   nb = read_line(pols);
   nbc = strlen(pols);

   printf("Number of characters : %d\n",nb);
   printf("Number of characters : %d\n",nbc);

   printf("Your system : ");
   for(k=0; k<nb; k++) printf("%c", pols[k]);
   printf("\n");

   fail = scan_number_of_variables(nb,pols,&dim);

   printf("fail : %d, dimension : %d\n",fail,dim);

   adafinal();

   return 0;
}

int read_line ( char *s )
{
   int c,i;

   for (i=0; (c=getchar()) != EOF && c != '\n' && i < MAXLEN; ++i) s[i] = c;
   s[i] = '\0';

   return i;
}
