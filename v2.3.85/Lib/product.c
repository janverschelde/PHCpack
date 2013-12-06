/* file product.c contains definitions of the prototypes of product.h */

#include <stdio.h>

extern void adainit( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal( void );

int supporting_set_structure ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc(110,a,b,c);

   return fail;
}

int write_set_structure ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc(111,a,b,c);

   return fail;
}

int linear_product_root_count ( int *r )
{
   int fail,*b;
   double *c;

   fail = _ada_use_c2phc(112,r,b,c);

   return fail;
}

int random_linear_product_system ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc(113,a,b,c);

   return fail;
}

int solve_linear_product_system ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc(114,a,b,c);

   return fail;
}

int clear_set_structure ( void )
{
   int fail,*a,*b;
   double *c;

   fail = _ada_use_c2phc(115,a,b,c);

   return fail;
}
