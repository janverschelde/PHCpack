/* This file "outputData.cpp" contains the definitions of the operations
 * to process the output of DEMiCs. */

#include "outputData.h"

int demics_allocate_lifting ( int nbrsup, int* crdsup )
{
   int fail;
   double *c;

   fail = _ada_use_c2phc(834,&nbrsup,crdsup,c);

   return fail;
}

int demics_assign_lifting ( int idxsup, int idxpnt, double value )
{
   int fail;

   fail = _ada_use_c2phc(835,&idxsup,&idxpnt,&value);

   return fail;
}

int demics_retrieve_lifting ( int idxsup, int idxpnt, double* value )
{
   int fail;

   fail = _ada_use_c2phc(836,&idxsup,&idxpnt,value);

   return fail;
}

int demics_clear_lifting ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(837,a,b,c);

   return fail;
}

int demics_append_cell_indices ( std::string strcell )
{
   int fail;
   double *c;

   int lenstrcell = strcell.length();

   int cellchars[lenstrcell];

   for(int k=0; k<lenstrcell; k++)
      cellchars[k] = (int) strcell[k];

   fail = _ada_use_c2phc(838,&lenstrcell,cellchars,c);

   return fail;
}

int demics_retrieve_cell_indices ( int idx, char* strcell )
{
   int fail;
   int buf[256];
   int nbr = idx;
   double *c;

   fail = _ada_use_c2phc(839,&nbr,buf,c);

   for(int k=0; k<nbr; k++) strcell[k] = (char) buf[k];
   strcell[nbr] = '\0';

   return fail;
}

int demics_clear_cell_indices ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(840,a,b,c);

   return fail;
}

int demics_store_mixed_volume ( int mv )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc(841,&mv,b,c);

   return fail;
}

int demics_retrieve_mixed_volume ( int* mv )
{
   int *b,fail;
   double *c;

   fail = _ada_use_c2phc(842,mv,b,c);

   return fail;
}
