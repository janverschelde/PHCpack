#include "idle_queue.h"
#include <stdlib.h>
#include <stdio.h>


IDLE_ELEMENT* addslv (IDLE_ELEMENT * listp, int data) 
{
   IDLE_ELEMENT *lp = listp;

   if (listp != NULL) 
   {
     while (listp -> link != NULL)
       listp = (IDLE_ELEMENT *)listp -> link;
     listp -> link = (IDLE_ELEMENT  *) malloc (sizeof (IDLE_ELEMENT));
     listp = (IDLE_ELEMENT *)listp -> link;
     listp -> link = NULL;
     listp -> data = data;
     return lp;
   }
   else 
   {
     listp = (IDLE_ELEMENT  *) malloc (sizeof(IDLE_ELEMENT));
     listp -> link = NULL;
     listp -> data = data;
     return listp;
   }
}

IDLE_ELEMENT* removeslv (IDLE_ELEMENT *lp) 
{
   IDLE_ELEMENT * tempp;
   int temp = lp->data;
   /* printf ("Element removed is %d\n", temp);  */
   tempp = (IDLE_ELEMENT *)lp -> link;
   free (lp);
   return tempp;
}

int num_idle(IDLE_ELEMENT * lp) 
{
   if (lp != NULL)
     return 1+num_idle((IDLE_ELEMENT *)lp->link);
   else
     return 0;
}
