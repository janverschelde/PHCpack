// plain C interface to DEMiCs, to call from Ada

void write_data
 ( int dimension, int nsupports,
   int* mixtype, int* cardsup, int *coordinates );
/*
 * DESCRIPTION :
 *   Writes the input data to screen.
 *
 * ON ENTRY :
 *   verbose    if 1, then data is printed, if 0, remains silent;
 *   dimension  length of the points in each support set;
 *   nsupports  number of distinct support sets;
 *   mixtype    number of occurrences of each set;
 *   cardsup    cardinalities of the support sets. */

extern "C" int demicsrun
 ( int verbose, int dimension, int nsupports,
   int* mixtype, int* cardsup, int *coordinates );
/*
 * DESCRIPTION :
 *   Calls DEMiCs to enumerate all mixed cells.
 *
 * ON ENTRY :
 *   verbose    if 1, then data is printed, if 0, remains silent;
 *   dimension  length of the points in each support set;
 *   nsupports  number of distinct support sets;
 *   mixtype    number of occurrences of each set;
 *   cardsup    cardinalities of the support sets. */
