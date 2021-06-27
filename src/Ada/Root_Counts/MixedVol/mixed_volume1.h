/* The file "mixed_volume1.h" contains the main program to compute the mixed
 * volume of a tuple of supports, usually invoked after some preprocessing. */

#ifndef _Mixed_Vol1_
#define _Mixed_Vol1_

#include "cell_stack1.h"

void MixedVol 
 ( int nVar, int nSpt, int CellSize, int *SptType, int *SptIdx, int **Spt,
   double *lft, int *nbCells, CellStack *MCells, int *MVol );
/*
 * DESCRIPTION :
 *   Computes the mixed volume of the given lifted supports.
 *
 * ON ENTRY :
 *   nVar          ambient dimension, number of variables;
 *   nSpt          number of different supports;
 *   CellSize      number of indices in one mixed cell;
 *   SptType       type of each support;
 *   SptIdx        SptIdx[i] indicates the start of i-th support in Spt;
 *   Spt           coordinates of the points in the supports;
 *   lft           lifting of the polynomial system.
 *
 * ON RETURN :
 *   nbCells       number of mixed cells;
 *   MCells        mixed cells of a regular mixed-cell configuration;
 *   MVol          the mixed volume of the polynomial system.   */

int gcd ( int r1, int r2, int *k, int *l );
/*
 * DESCRIPTION :
 *   Returns the gcd of the numbers r1 and r2
 *   with gcd(r1,r2) = k*r1 + l*r2 on return. */

int cell_size ( int nSpt, int *SptType );
/*
 * DESCRIPTION :
 *   Returns the number of indices it takes to represent a mixed cell.
 *
 * ON ENTRY :
 *   nSpt          number of different supports;
 *   SptType       #occurrences of each support. */

void CellVol
 ( int nVar, int nSpt, int **Spt, int *SptType, int *Cell, int *Vol );
/*
 * DESCRIPTION :
 *   Returns the mixed volume of the current cell.
 *
 * ON ENTRY :
 *   nVar         dimension of the ambient space;
 *   nSpt         number of different supports;
 *   Spt          coordinates of the points in the supports;
 *   Cell         a mixed cell.
 *
 * ON RETURN :
 *   Vol          accumulated mixed volume.   */

int solve_linear_system ( int n, double A[n][n], double b[n] );
/*
 * DESCRIPTION :
 *   Solves the linear system A*x = b.
 *
 * NOTE :
 *   This routine is used only to compute the inner normal
 *   when writing the mixed-cell configuration to file.
 *
 * ON ENTRY :   
 *   n        dimension of the linear system;
 *   A        an n-by-n matrix;
 *   b        right hand side vector of size n.
 *
 * ON RETURN :
 *   b        solution x to A*x = b,
 *   if the integer returned by this function is 1,
 *   otherwise if A is singular, the function returns 0. */

void write_cells
 ( int ncfn, char output_file[ncfn],
   int nVar, int nSpt, int *SptType, int **Spt, double *lft,
   int CellSize, int nMCells, CellStack *MCells );
/*
 * DESCRIPTION :
 *   Writes the regular mixed-cell configuration to file,
 *   in a format suitable for PHCpack.
 *
 * ON ENTRY :
 *   ncfn          number of characters in the name of the output file;
 *   output_file   name of the output file for the mixed-cell configuration;
 *   nVar          ambient dimension, number of variables;
 *   nSpt          number of different supports;
 *   SptType       type of each support;
 *   Spt           coordinates of the points in the supports;
 *   lft           lifting of the polynomial system;
 *   CellSize      number of points in a mixed cell;
 *   MCells        mixed cells of a regular mixed-cell configuration.
 * 
 * ON RETURN :
 *   MCells        mixed cells are popped off the stack ... */

#endif
