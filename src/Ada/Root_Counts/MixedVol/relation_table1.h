/* The file "relation_tabl1e.h" collects the prototypes of routines to
 * create a so-called "relation table" which records all edges on the
 * lower hulls of all pairwise Minkowski sums. */

#ifndef _RelationTable1_
#define _RelationTable1_  

#include "zero_index_tree1.h"
#include "index_tree_lp1.h"

void RelTable ( int nVar, int nSpt, int **Spt, int *SptIdx,
                double *lft, int **RelTab, L0_IML *L0 );
/*
 * DESCRIPTION :
 *   Sets up the relation table for the given lifted supports.
 *   This is the main routine which calls the other functions
 *   with prototypes listed below.
 *
 * ON ENTRY :
 *   nVar       ambient dimension;
 *   nSpt       number of different supports;
 *   Spt        supports, in array of range 0..nSpt-1;
 *   SptIdx     the i-th supports starts at SptIdx[i] in Spt,
 *              SptIdx[nSpt] equals the total #points in Spt;
 *   lft        lifting for each point in Spt;
 *   RelTab     matrix of dimension equal to the total #points;
 *   L0         initialized level zero tree.
 *
 * ON RETURN :
 *   RelTab     RelTab[i][j] == 1 if points i and j span an edge
 *              on the lower hull of the polytope (if i and j are
 *              from the same support set) or on the lower hull
 *              of the Minkowski sum of the polytopes spanned by
 *              the support set which contains point i and the
 *              other support set which contains point j,
 *              RelTab[i][j] == 0 otherwise,
 *              note that this matrix is symmetric;
 *   L0         updated tree with stored LP information. */

void RlTbLP2_a
 ( int ma, int na, double **a, int nv1, double *c,
   int *LPidx, int *Bidx, double *x, double **Binv, int *info );
/*
 * DESCRIPTION :
 *   Applies LP to find a nondegenerate solutions to perform a 1-point test.
 *
 * ON ENTRY :
 *   ma         number of involved constraints;
 *   na         number of variables (= nVar) plus one;
 *   a          matrix of inequality constraints,
 *              as many rows as total #points plus one,
 *              as many columns as nVar + 2;
 *   nv1        #variables plus one;
 *   c          vector of range 0..nVar+1,
 *              used to determine outgoing constraint;
 *   LPidx      vector of dimension total #points plus one,
 *              index to the involved constaints;
 *   Bidx       index vector of range 0..nVar+1;
 *   x          solution vector of range 0..nVar+1;
 *   Binv       square matrix of dimension nVar+2.
 *
 * ON RETURN :
 *   Bidx       updated index vector;
 *   x          updated solution vector;
 *   Binv       updated matrix;
 *   info       information about the rank. */

void RlTbLP2_e
 ( int ma, int na, int NumCol, double **a, int nv1, int *LPidx,
   int *Bidx, double *x, double **Binv, int *info );
/*
 * DESCRIPTION :
 *   Applies LP to extend the point using a 1-point test,
 *   returns a nondegenerate solution for 1-point tests.
 *
 * ON ENTRY :
 *   ma         number of involved constraints;
 *   na         number of variables (= nVar), original dimension;
 *   NumCol     total number of points in the support;
 *   a          matrix of inequality constraints,
 *              as many rows as total #points plus one,
 *              as many columns as nVar + 2;
 *   nv1        #variables plus one;
 *   LPidx      vector of dimension total #points plus one,
 *              index to the involved constraints;
 *   Bidx       index vector of range 0..nVar+1;
 *   x          solution vector of range 0..nVar+1;
 *   Binv       square matrix of dimension nVar+2.
 *
 * ON RETURN :
 *   a          extended matrix with the "epsilon variable";
 *   Bidx       updated index vector;
 *   x          updated solution vector;
 *   Binv       updated matrix;
 *   info       information about the rank. */

void dlp2_1pt_i
 ( int ma, int na, double **a, int nv1, double *c, int *LPidx, int FixPt,
   int TstPt, int *Bidx, double *x, double **Binv, int **RelTab );
/*
 * DESCRIPTION :
 *   Applies LP to perform a 1-point test in case the current support
 *   is not the first one, for the variable info < 0.
 *
 * ON ENTRY :
 *   ma         number of involved constraints;
 *   na         number of variables (= nVar);
 *   a          matrix of inequality constraints,
 *              as many rows as total #points plus one,
 *              as many columns as nVar + 2;
 *   nv1        #variables plus one;
 *   c          vector of range 0..nVar+1,
 *              used to determine outgoing constraint;
 *   LPidx      vector of dimension total #points plus one;
 *   FixPt      index to a point in the supports;
 *   TstPt      index to an involved constraint;
 *   Bidx       index vector of range 0..nVar+1;
 *   x          solution vector of range 0..nVar+1;
 *   Binv       square matrix of dimension nVar+2;
 *   RelTab     matrix of dimension equal to the total #points.
 *
 * ON RETURN :
 *   Bidx       updated index vector;
 *   x          updated solution vector;
 *   Binv       updated matrix;
 *   RelTab     updated relation table. */

void dlp2_1pt_s
 ( int ma, int na, double **a, int nv1, double *c, int *LPidx, int FixPt,
   int TstPt, int *Bidx, double *x, double **Binv, int **RelTab, L0_IML *L0 );
/*
 * DESCRIPTION :
 *   Applies LP to perform a 1-point test for the first support,
 *   in case the variable info < 0.
 *
 * ON ENTRY :
 *   ma         number of involved constraints;
 *   na         number of variables (= nVar);
 *   a          matrix of inequality constraints,
 *              as many rows as total #points plus one,
 *              as many columns as nVar + 2;
 *   nv1        #variables plus one;
 *   c          vector of range 0..nVar+1,
 *              used to determine outgoing constraint;
 *   LPidx      vector of dimension total #points plus one,
 *              index to the involved constraints;
 *   FixPt      index to a point in the supports;
 *   TstPt      index to an involved constraint;
 *   Bidx       index vector of range 0..nVar+1;
 *   x          solution vector of range 0..nVar+1;
 *   Binv       square matrix of dimension nVar+2;
 *   RelTab     matrix of dimension equal to the total #points;
 *   L0         current level zero tree.
 *
 * ON RETURN :
 *   Bidx       updated index vector;
 *   x          updated solution vector;
 *   Binv       updated matrix;
 *   RelTab     updated relation table;
 *   L0         LP data with Bidx, x, and Binv. */
   
void dlp1_1pt_i
 ( int ma, int na, double **a, int nv1, double *c, int *LPidx, int FixPt,
   int TstPt, int *Bidx, double *x, double **Binv, int **RelTab );
/*
 * DESCRIPTION :
 *   Applies LP to perform a 1-point test when the current support
 *   is not the first one, and for variable info >= 0.
 *
 * ON ENTRY :
 *   ma         number of involved constraints;
 *   na         number of variables;
 *   a          matrix of inequality constraints,
 *              as many rows as total #points plus one,
 *              as many columns as nVar + 2;
 *   nv1        #variables plus one;
 *   c          vector of range 0..nVar+1,
 *              used to determine outgoing constraint;
 *   LPidx      vector of dimension total #points plus one,
 *              index to the involved constraints;
 *   FixPt      index to a point in the supports;
 *   TstPt      index to an involved constraint;
 *   Bidx       index vector of range 0..nVar+1;
 *   x          solution vector of range 0..nVar+1;
 *   Binv       square matrix of dimension nVar+2;
 *   RelTab     matrix of dimension equal to the total #points.
 *
 * ON RETURN :
 *   Bidx       updated index vector;
 *   x          updated solution vector;
 *   Binv       updated matrix;
 *   RelTab     updated relation table. */

void dlp1_1pt_s
 ( int ma, int na, double **a, int nv1, double *c, int *LPidx, int FixPt,
   int TstPt, int *Bidx, double *x, double **Binv, int **RelTab, L0_IML *L0 );
/*
 * DESCRIPTION :
 *   Applies LP to perform a 1-point test for the first support
 *   and for variable info >= 0.
 *
 * ON ENTRY :
 *   ma         number of involved constraints;
 *   na         number of variables;
 *   a          matrix of inequality constraints,
 *              as many rows as total #points plus one,
 *              as many columns as nVar + 2;
 *   nv1        #variables plus one;
 *   c          vector of range 0..nVar+1,
 *              used to determine outgoing constraint;
 *   LPidx      vector of dimension total #points plus one,
 *              index to the involved constraints;
 *   FixPt      index to a point in the supports;
 *   TstPt      index to an involved constraint;
 *   Bidx       index vector of range 0..nVar+1;
 *   x          solution vector of range 0..nVar+1;
 *   Binv       square matrix of dimension nVar+2;
 *   RelTab     matrix of dimension equal to the total #points;
 *   L0         current level zero tree.
 *
 * ON RETURN :
 *   x          updated vector;
 *   Binv       updated matrix;
 *   RelTab     updated relation table;
 *   L0         saved LP data Bidx, x, and Binv. */

#endif
