/* simple test on running a Pade continuation */

#include <stdio.h>
#include "padcon.h"

int main ( void )
{
   adainit();

   printf("Tuning the homotopy continuation parameters ...\n");

   padcon_set_default_parameters();

   adafinal();

   return 0;
}
