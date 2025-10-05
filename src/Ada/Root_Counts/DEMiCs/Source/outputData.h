/* This file "outputData.h" contains the prototypes of the operations
 * to process the output of DEMiCs.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __OUTPUTDATA_H__
#define __OUTPUTDATA_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal( void );
#endif

#include <string>

int demics_allocate_lifting ( int nbrsup, int* crdsup );
/*
 * DESCRIPTION :
 *   Allocates memory for the lifting values.
 *   On return is the failure code, which equals 0 if all went well.
 *
 * ON ENTRY :
 *   nbrsup   the number of distinct supports;
 *   crdsup   cardinality of each support set. */

int demics_assign_lifting ( int idxsup, int idxpnt, double value );
/*
 * DESCRIPTION :
 *   Assigns a lifting value for a point defined by two indices.
 *   On return is the failure code, which equals 0 if all went well.
 *
 * ON ENTRY :
 *   idxsup   index to a support set;
 *   idxpnt   index to the point in support set with index idxsup;
 *   value    lifting value for the point defined by the
 *            two indices idxsup and idxpnt. */

int demics_retrieve_lifting ( int idxsup, int idxpnt, double* value );
/*
 * DESCRIPTION :
 *   Retrieves the lifting value of a point defined by two indices.
 *   On return is the failure code, which equals 0 if all went well.
 *
 * ON ENTRY :
 *   idxsup   index to a support set;
 *   idxpnt   index to the point in support set with index idxsup.
 *
 * ON RETURN :
 *   value    lifting value for the point defined by the
 *            two indices idxsup and idxpnt. */

int demics_clear_lifting ( void );
/*
 * DESCRIPTION :
 *   Deallocates the memory used for the lifting of the supports.
 *   On return is the failure code, which equals 0 if all went well. */

int demics_append_cell_indices ( std::string strcell );
/*
 * DESCRIPTION :
 *   Appends the string representation of the cell indices. */

int demics_retrieve_cell_indices ( int idx, char* strcell );
/*
 * DESCRIPTION :
 *   Returns the string at position idx.
 *   The strcell should have allocated sufficient space for the result. */

int demics_clear_cell_indices ( void );
/*
 * DESCRIPTION :
 *   Deallocates memory occupied by the cell indices. */

int demics_store_mixed_volume ( int mv );
/*
 * DESCRIPTION :
 *   Stores the mixed volume mv. */

int demics_retrieve_mixed_volume ( int* mv );
/*
 * DESCRIPTION :
 *   Returns the stored mixed volume in mv. */

#endif
