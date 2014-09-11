/* test on conversion between strings and arrays of numbers */

#include <stdio.h>
#include "lists_and_strings.h"

void test_intlist2str ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a sequence of integers,
 *   writes the integers to a string and then reads
 *   the string back into a sequence of integers. */
 
void test_dbllist2str ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a sequence of doubles,
 *   writes the doubles to a string and then reads
 *   the string back into a sequence of doubles. */
 
int main ( int argc, char *argv[] )
{
   test_intlist2str();
   test_dbllist2str();
 
   return 0;
}

void test_intlist2str ( void )
{
   int i,n,nbc;
   char s[80];

   printf("give the number of integers : ");
   scanf("%d",&n);
   {
      int data[n];
      for(i=0; i<n; i++)
      {
         printf("-> give element %d : ",i);
         scanf("%d",&data[i]);
      }
      printf("The %d numbers : ",n);
      for(i=0; i<n; i++) printf(" %d",data[i]);
      printf("\n");
      nbc = intlist2str(n,data,s);
   }   
   printf("The string : %s\n",s);
   {
      int ict = itemcount(s);
      int nbr[ict];
      int i;

      printf("number of items : %d\n",ict);
      str2intlist(n,s,nbr);
      printf("the %d numbers :",ict);
      for(i=0; i<ict; i++) printf(" %d",nbr[i]);
      printf("\n");
   }
}

void test_dbllist2str ( void )
{
   int i,n,nbc;
   char s[80];

   printf("give the number of doubles : ");
   scanf("%d",&n);
   {
      double data[n];
      for(i=0; i<n; i++)
      {
         printf("-> give element %d : ",i);
         scanf("%lf",&data[i]);
      }
      printf("The %d numbers : ",n);
      for(i=0; i<n; i++) printf(" %.16le",data[i]);
      printf("\n");
      nbc = dbllist2str(n,data,s);
   }   
   printf("The string : %s\n",s);
   {
      int ict = itemcount(s);
      double nbr[ict];
      int i;

      printf("number of items : %d\n",ict);
      str2dbllist(n,s,nbr);
      printf("the %d numbers :",ict);
      for(i=0; i<ict; i++) printf(" %.16le",nbr[i]);
      printf("\n");
   }
}
