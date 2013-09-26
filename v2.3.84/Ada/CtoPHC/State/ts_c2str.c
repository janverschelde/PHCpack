/* simple test to have phc write C strings */

#include<stdio.h>
#include<stdlib.h>

extern void adainit();
extern int _ada_use_c2phc ( int job, int *a, int *b, double *c );
extern void adafinal();

int main ( int argc, char *argv[] )
{
   double *c;
   int fail,choice,i;
   char skip_newline,ch;

   printf("\nTesting if phc can C strings...\n");

   adainit();
   {
      int n,*s;
      do
      { 
         printf("Give the number of elements : "); scanf("%d",&n);
         s = (int*)calloc(n,sizeof(int));
         skip_newline = getchar();
         printf("Give %d characters : ",n);
         for(i=0; i<n; i++)
         {
            scanf("%c",&ch);
            s[i] = (int) ch;
         }
         printf("The array of integers : ");
         for(i=0; i<n; i++) printf(" %d",s[i]); printf("\n");
         printf("The array of characters : ");
         for(i=0; i<n; i++)
         {
            ch = (char) s[i];
            printf(" %c",ch);
         }
         printf("\n");
         fail = _ada_use_c2phc(158,&n,s,c);
         printf("use_c2phc returned %d as fail value\n",fail);
         printf("Give a number (1 to continue, 0 to exit) : ");
         scanf("%d",&choice);
         skip_newline = getchar();
         free(s);
      } while (choice > 0);
   }
   adafinal();

   return 0;
}
