/* The file "one_level_lp1.h" collects the prototypes to perform tests at
 * one level, as called by the main program "mixed_volume1.c". */
 
#ifndef _OneLevelLP1_
#define _OneLevelLP1_
 
#include "index_tree_lp1.h"

void one_level_LP 
 ( int Strt1Pt, int End1Pt, int *PtIn, int LPdim, double **A,
   int nVar, double *x, double **Binv, int *Bidx, IT_LP *ItLp );
/*
 * DESCRIPTION :
 *   This routine is called in the main program in "mixed_volume.c"
 *   and uses the other functions below in its implementation.
 *
 * ON ENTRY :
 *   Strt1Pt   index where to start the 1-point tests;
 *   End1Pt    index where to end the 1-point tests;
 *   PtIn      one row of a matrix with as many columns as the
 *             total number of points in the supports;
 *   LPDim     current dimension of the LP problem,
 *             its value is used as "na" in the routines below;
 *   A         matrix of constraints at the current level;
 *   nVar      number of variables;
 *   x         solution vector of dimension nVar;
 *   Binv      basis inverse, matrix of dimension nVar;
 *   Bidx      index vector of dimension nVar;
 *   ItLp      index tree with data for LP problems.
 *
 * ON RETURN :
 *   PtIn      updated row of indices;
 *   x         updated solution vector;
 *   Binv      updated basis inverse;
 *   Bidx      updated index vector;
 *   ItLp      updated index tree with data for LP problems. */

void dnulp2_a
 ( int ma, int na, double **a, int nVar, double *c, int *Bidx,
   double *x, double **Binv, int *info );
/*
 * DESCRIPTION :
 *   Auxiliary routine called by one_level_LP,
 *   invoked the first time when performing a 1-point test.
 *
 * ON ENTRY :
 *  ma         #constraints involved in the LP problem;
 *  na         dimension of the LP problem;
 *  a          matrix with coefficients of the inequality constraints;
 *  nVar       number of variables;
 *  c          vector of dimension nVar+2,
 *             used to determine the outgoing constraint;
 *  Bidx       index vector of dimension nVar;
 *  x          solution vector of dimension nVar;
 *  Binv       basis inverse, matrix of dimension nVar.
 *
 * ON RETURN :
 *  Bidx       updated index vector;
 *  x          updated solution vector;
 *  Binv       updated basis inverse;
 *  info       information about the dimension of the problem. */

void dlp2_1pts
 ( int ma, int na, double **a, int nVar, double *c, int TstPt,
   int *Bidx, double *x, double **Binv, int *PtIn, IT_LP *ItLp );
/*
 * DESCRIPTION :
 *   Auxiliary routine called by one_level_LP,
 *   invoked after the first time, when dnulp2_a returns info >= 0.
 *
 * ON ENTRY :
 *  ma         #constraints involved in the LP problem;
 *  na         dimension of the LP problem;
 *  a          matrix with coefficients of the inequality constraints;
 *  nVar       number of variables;
 *  c          vector of dimension nVar+2,
 *             used to determine the outgoing constraint;
 *  TstPt      index of a point ranging between Strt1Pt+1 and End1Pt
 *             of the routine one_level_LP from above;
 *  Bidx       index vector of dimension nVar;
 *  x          solution vector of dimension nVar;
 *  Binv       basis inverse, matrix of dimension nVar;
 *  PtIn       one row of a matrix with as many columns as the
 *             total number of points in the supports;
 *  ItLp       index tree with data for LP problems.
 *
 * ON RETURN :
 *  Bidx       updated index vector;
 *  x          updated solution vector;
 *  Binv       updated basis inverse;
 *  PtIn       vector with one updated entry;
 *  ItLp       updated index tree with data for LP problems. */

void dlp1_1pts
 ( int ma, int na, double **a, int nVar, double *c, int TstPt,
   int *Bidx, double *x, double **Binv, int *PtIn, IT_LP *ItLp );
/*
 * DESCRIPTION :
 *   Auxiliary routine called by one_level_LP,
 *   invoked after the first time, when dnulp2_a returns info < 0.
 *
 * ON ENTRY :
 *  ma         #constraints involved in the LP problem;
 *  na         dimension of the LP problem;
 *  a          matrix with coefficients of the inequality constraints;
 *  nVar       number of variables;
 *  c          vector of dimension nVar+2,
 *             used to determine the outgoing constraint;
 *  TstPt      index of a point ranging between Strt1Pt+1 and End1Pt
 *             of the routine one_level_LP from above;
 *  Bidx       index vector of dimension nVar;
 *  x          solution vector of dimension nVar;
 *  Binv       basis inverse, matrix of dimension nVar;
 *  PtIn       one row of a matrix with as many columns as the
 *             total number of points in the supports;
 *  ItLp       index tree with data for LP problems.
 *
 * ON RETURN :
 *  Bidx       updated index vector;
 *  x          updated solution vector;
 *  Binv       updated basis inverse;
 *  PtIn       vector with one updated entry;
 *  ItLp       updated index tree with data for LP problems. */

void Sort ( int n, int *a );
/*
 * DESCRIPTION :
 *   Sorts the first n numbers in the array a in ascending order.
 *
 * ON ENTRY :
 *   n         dimension of the array a;
 *   a         sequence of n integer numbers.
 *
 * ON RETURN :
 *   a         sequence sorted in ascending order. */

#endif 
