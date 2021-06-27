/* This file constains the prototypes of the routines called by the main
 * program in "mixed_volume1.h", used to define the LP problems. */

#ifndef _FORMLP1_
#define _FORMLP1_

#include "index_tree_lp1.h"

void form_LP
 ( int nVar, int nSpt, int *SptType, int *SptIdx, int **RelTab, 
   double **ElmMtrx, int *NonZero, int *Lvl2CoDim, double ***A,
   int Lvl, LPdata *ptr, int *Cell, int *Cell_Orig, int *Lvl2LPdim,
   int *Lvl2Spt, int *FixFrstPt, int *FixLstPt, int *MinNumPt, 
   int **PtIn, int *Strt1Pt, int *End1Pt, int **ToOrig,
   int **Pre2Cur, int *LPdim, int *FrstPt, int *LstPt,
   int *Cur2Pre, double *x, double **Binv, int *Bidx, int *info );
/*
 * DESCRIPTION :
 *   Defines the LP problems when the value of Lvl > 1.
 *
 * ON ENTRY :
 *   nVar        number of variables, ambient dimension;
 *   nSpt        number of different supports;
 *   SptType     array of range 0..nSpt with #occurrences of each support;
 *   SptIdx      SptIdx[i] indicates the start of the i-th support;
 *   RelTab      the relation table records edges on lower hull of
 *               each support or on pairwise Minkowski sums;
 *   ElmMtrx     elimination matrix has dimension the size of a cell;
 *   NonZero     array of dimension equal to the size of a mixed cell,
 *               but only one entry will be updated;
 *   Lvl2CoDim   array of dimension equal to the size of a mixed cell,
 *               counts the number of variables eliminated,
 *               only the element at position Lvl is used;
 *   A           3-dimensional matrix, first dimension is the cell size,
 *               then every A[i] is a coefficient matrix of constraints,
 *               of dimension #constraints times nVar + 1;
 *   Lvl         index of the current level;
 *   ptr         pointer to LP data, used only for consultation;
 *   Cell        indices of the current mixed cell;
 *   Cell_Orig   indices to the original mixed cell;
 *   Lvl2LPdim   array of cell size, only entries at Lvl and Lvl-1 are used;
 *   Lvl2Spt     array of cell size, only entry at Lvl is used;
 *   FixFrstPt   array of cell size, only entry at Lvl is used;
 *   FixLstPt    array of cell size, only entry at Lvl is used;
 *   MinNumPt    array of cell size, only entry at Lvl is used;
 *   PtIn        2-dimensional matrix with as many rows as the cell size
 *               and as many columns as the total number of points,
 *               only the points in PtIn[Lvl] and PtIn[Lvl-1] are accessed;
 *   ToOrig      2-dimensional integer matrix, as many rows as cell size,
 *               and as many columns as the total number of points;
 *   Pre2Cur     2-dimensional integer matrix, as many rows as cell size,
 *               and as many columns as the total number of points;
 *   FrstPt      array of size of a mixed cell,
 *               serves as index vector to first point in current support;
 *   LstPt       array of size of mixed cell,
 *               serves as index vector to last point in current support;
 *   Cur2Pre     array of dimension equal to total number of points;
 *   x           solution vector of dimension nVar;
 *   Binv        matrix of dimension nVar, basis inverse;
 *   Bidx        index vector for points in Pre2Cur, of dimension nVar.
 *
 * ON RETURN :
 *   ElmMtrx     updated elimination matrix;
 *   NonZero     entry at Lvl2CoDim[Lvl] is updated;
 *   Cell_Orig   updated indices of the original mixed cell;
 *   PtIn        updated entries at row Lvl;
 *   Strt1Pt     index for where to start the 1-point tests;
 *   End1Pt      index for where to end the 1-point tests;
 *   ToOrig      updated entries at rows Lvl-1 and Lvl;
 *   Pre2Cur     updated entries at row 1 and row Lvl;
 *   LPdim       equals Lvl2LPdim[Lvl];
 *   Cur2Pre     updated index vector;
 *   x           updated solution vector;
 *   Binv        updated basis inverse;
 *   Bidx        updated index vector;
 *   info        0 if normal return,
 *               1 if not enough points to extend, must backtrack. */

void form_LP1
 ( int nVar, int nSpt, int *SptType, int *SptIdx, int **RelTab,
   double ***A, int Lvl, int *Cell, int *Lvl2LPdim, int *FixLstPt,
   int *MinNumPt, int **PtIn, int **ToOrig, int **Pre2Cur,
   int *FrstPt, int *LstPt, int *info );
/*
 * DESCRIPTION :
 *   Defines the LP problems when the value of Lvl <= 1.
 *
 * ON ENTRY :
 *   nVar        number of variables, ambient dimension;
 *   nSpt        number of different supports;
 *   SptType     array of range 0..nSpt with #occurrences of each support;
 *   SptIdx      SptIdx[i] indicates the start of the i-th support;
 *   RelTab      the relation table records edges on lower hull of
 *               each support or on pairwise Minkowski sums;
 *   A           3-dimensional matrix, first dimension is the cell size,
 *               then every A[i] is a coefficient matrix of constraints,
 *               of dimension #constraints times nVar + 1;
 *   Lvl         index of the current level;
 *   Cell        indices of the current mixed cell;
 *   Lvl2LPdim   array of cell size, only entry at Lvl-1 is used;
 *   FixLstPt    array of cell size, only entry at Lvl is used;
 *   MinNumPt    array of cell size, only entry at Lvl is used;
 *   PtIn        2-dimensional matrix with as many rows as the cell size
 *               and as many columns as the total number of points,
 *               only the points in PtIn[Lvl] are accessed;
 *   ToOrig      2-dimensional integer matrix, as many rows as cell size,
 *               and as many columns as the total number of points;
 *   Pre2Cur     2-dimensional integer matrix, as many rows as cell size,
 *               and as many columns as the total number of points;
 *   FrstPt      array of size of a mixed cell,
 *               serves as index vector to first point in current support;
 *   LstPt       array of size of mixed cell,
 *               serves as index vector to last point in current support.
 *
 * ON RETURN :
 *   A           constraints updated in A[Lvl], using A[Lvl-1];
 *   PtIn        updated entries in PtIn[Lvl];
 *   ToOrig      updated entries in ToOrig[Lvl];
 *   Pre2Cur     updated entries in Pre2Cur[Lvl];
 *   FrstPt      set FrstPt[Lvl] to zero;
 *   LstPt       set LstPt[Lvl] to #constraints involved;
 *   info        0 if normal return,
 *               1 if not enough points to extend, must backtrack. */

#endif
