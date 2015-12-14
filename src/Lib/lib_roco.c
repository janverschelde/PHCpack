/* testing linear-product root counts and random linear-product systems */

#include <stdio.h>
#include <stdlib.h>
#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"
#include "product.h"

int compute_root_count ( int *r );
/*
 * DESCRIPTION :
 *   Prompts the user for a polynomial system, constructs a supporting
 *   set structure, and returns the corresponding root count in r. */ 

int compute_m_homogeneous_Bezout_number ( int *b );
/*
 * DESCRIPTION :
 *   Prompts the user for a polynomial system and computes a m-homogeneous
 *   Bezout number *b based on a partition. */

int construct_start_system ( int r );
/*
 * DESCRPTION :
 *   Constructs a random linear-product system based on the supporting
 *   set structure and with root count in r.
 *
 * REQUIRED : compute_root_count() was executed. */

int main ( void )
{
   int r;

   printf("\nTesting linear-product root counts and systems...\n");

   adainit();

   printf("MENU for computing Bezout bounds :\n");
   printf("  1. test m-homogeneous Bezout numbers;\n");
   printf("  2. test general linear-product root counts.\n");
   printf("Type 1 or 2 to make your choice : ");
   scanf("%d",&r);

   if(r == 1) 
   {
      compute_m_homogeneous_Bezout_number(&r);
   }
   else if(r == 2)
   {
      compute_root_count(&r);
      construct_start_system(r);
      clear_set_structure();
   }
   else
      printf("Wrong choice, please try again...\n");

   adafinal();

   return 0;
}

int compute_root_count ( int *r )
{
   int fail,inout;

   fail = syscon_read_standard_system();
   printf("\nThe system in the container : \n");
   fail = syscon_write_standard_system();
   fail = supporting_set_structure();
   printf("\nA supporting set structure : \n");
   fail = write_set_structure();
   {
      int nc;
      char ss[200];
      fail = set_structure_string(&nc,ss);
      printf("\nthe set structure string : \n%s\n",ss);
      printf("clearing the set structure ...\n");
      clear_set_structure();     
      printf("parsing the set structure string ...\n");
      parse_set_structure(nc,ss);
      printf("after parsing the set structure string :\n");
      fail = write_set_structure();
   }
   inout = 1;
   fail = is_set_structure_supporting(&inout);
   if(inout == 1)
      printf("The set structure supports the polynomial system.\n");
   else
      printf("The set structure does not support the polynomial system.\n");
   fail = linear_product_root_count(r);
   printf("\nThe linear-product root count : %d\n",*r);

   return fail;
}

int compute_m_homogeneous_Bezout_number ( int *b )
{
   int fail,mbz,ncp,nbsols;
   int size=256;
   char partition[size];
   char answer;

   fail = syscon_read_standard_system();
   printf("\nThe system in the container : \n");
   fail = syscon_write_standard_system();

   fail = m_homogeneous_Bezout_number(&mbz,&ncp,partition);
   printf("\nA m-homogeneous Bezout number : %d\n",mbz);
   printf("with partition %s.\n",partition);

   printf("\nDo you wnat to give a partition ? (y/n) ");
   scanf("%c",&answer);
   if(answer == 'y')
   {
      scanf("%c",&answer); /* skip newline */
      printf("\nReading a partition ...\nEnter a string : ");
      {
         int c,i;
         for(i=0; (c=getchar()) != EOF && c != '\n' && i < size; ++i)
            partition[i] = c;
         partition[i] = '\0';
         ncp = i+1;
      }
      printf("-> your partition : %s\n",partition);
      fail = m_partition_Bezout_number(&mbz,ncp,partition);
      printf("the m-homogeneous Bezout number : %d\n",mbz);
   }
   printf("\nConstructing a m-homogeneous start system ...\n");
   fail = m_homogeneous_start_system(ncp,partition);
   printf("An m-homogeneous linear-product start system :\n");
   fail = syscon_write_standard_system();

   fail = solve_linear_product_system();
   printf("\nThe solutions : \n");
   fail = solcon_write_standard_solutions();
   fail = solcon_number_of_standard_solutions(&nbsols);
   if(mbz == nbsols)
      printf("\nComputed %d solutions, as many as the root count.\n",nbsols);
   else
      printf("\nNumber of solutions computed %d /= %d, the root count ?!!\n",
             nbsols,mbz);

   return fail;
}

int construct_start_system ( int r )
{
   int fail;
   int nbsols;

   fail = random_linear_product_system();
   printf("\nA random linear-product system :\n");
   fail = syscon_write_standard_system();
   fail = solve_linear_product_system();
   printf("\nThe solutions : \n");
   fail = solcon_write_standard_solutions();
   fail = solcon_number_of_standard_solutions(&nbsols);
   if(r == nbsols)
      printf("\nComputed %d solutions, as many as the root count.\n",nbsols);
   else
      printf("\nNumber of solutions computed %d /= %d, the root count ?!!\n",
             nbsols,r);

   return fail;
}
