/* The file "prepare_for_mv1.h" collects the prototypes for the routines used
 * to preprocess a tuple of supports for the mixed-volume computation. */

#ifndef _Preprocess_for_Mixed_Volume_Supports1_
#define _Preprocess_for_Mixed_Volume_Supports1_

void Pre4MV ( int nVar, int nSpt, int *nS, int *SptType,
              int **Spt, int *SptIdx, int **Vtx, int *VtxIdx, int *OldIdx );
/*
 * DESCRIPTION :
 *   As a preprocessing to the mixed-volume computation, this function
 *   extracts the non-vertex points and determines the type of mixture.
 *
 * ON ENTRY :
 *   nVar     the number of variables, embient dimension;
 *   nSpt     the number of supports;
 *   SptType  allocated memory for an array of nVar integers;
 *   Spt      array of points, of range 0..nSpt-1,
 *            with each point having nVar integer coordinates;
 *   SptIdx   array of range 0..nSpt, SptIdx[i] indicates the start
 *            of the i-th support, SptIdx[nVar] contains the total #points;
 *   Vtx      memory allocated for as many points as in Spt;
 *   VtxIdx   copy of SptIdx;
 *   OldIdx   memory allocated for total number of points.
 *
 * ON RETURN :
 *   nS       the number of different supports;
 *   SptType  type of each support, defined for range 0..nS-1;
 *   Vtx      vertex points for each support;
 *   VtxIdx   array of range 0..nS, VtxIdx[i] indicates the start of the
 *            i-th support in Vtx, VtxIdx[nVar] equals total #points;
 *   OldIdx   the position of each point in Vtx
 *            with respect to the original supports in Spt. */

void SortSpt ( int nVar, int nSpt, int **Spt, int *SptIdx );
/*
 * DECRIPTION :
 *   Sorts the point in each support in descending lexicographic order.
 *
 * ON ENTRY :
 *   nVar     the ambient dimension;
 *   nSpt     number of supports;
 *   Spt      the coordinates of the supports are in the rows;
 *   SptIdx   is the index of the 1st point of each support in Spt,
 *            SptIdx[nSpt]-1 is the index of the last point of the supports.
 *
 * ON RETURN :
 *   Spt      sorted supports. */

void LowerTriangular
 ( double **A, int jStrt, int jEnd, int na, int *rnk, int *ib );
/*
 * DESCRIPTION :
 *   Makes the matrix A lower triangular, called by NonVertex below.
 *
 * ON ENTRY :
 *   A        matrix to be triangulated;
 *   jStrt    index of start row;
 *   jEnd     index of end row;
 *   na       number of columns in A.
 *
 * ON RETURN :
 *   A        lower triangular matrix;
 *   rnk      rank of the matrix;
 *   ib       pivoting information. */

void RSimplex
 ( double **A, double *b, int jStrt, int jEnd, int m,
   double *xb, double *sb, int *ib, double **binv, int *info );
/*
 * DESCRIPTION :
 *   Applies the revised simplex method to solve
 *
 *     min e   subject to   -e <= 0 and -e + B y <= b
 *
 *   where (e,y) are variables.
 *   The constraint -e <= 0 corresponds to the jStrt-th row of A.
 *
 * ON ENTRY :
 *   A        active constraints range from row jStrt to jEnd;      
 *   b        right hand side vector;
 *   jStrt    index to the row in A marking the start of the active part;
 *   JEnd     index to the column in A marking the end of the active part;
 *   m        columns in A go up to m;
 *   xb       initial feasible solution;
 *   sb       workspace, initially declared as global variable,
 *            managed by the function NonVertex;
 *   ib       indices of the basis;
 *   binv     inverse of the initial base matrix,
 *            initialized to the identity matrix.
 *
 * ON RETURN :
 *   xb       updated basis vector;
 *   sb       updated workspace vector;
 *   ib       updated index vector;
 *   binv     updated basis inverse matrix;
 *   info     position of the variable epsilon when info>=0. */

void ReArrangeSpt
 ( int nVar, int **Spt, int *SptIdx, int **Vtx, int *VtxIdx,
   int *SptType, int *OldIdx );
/*
 * DESCRIPTION :
 *   Rearranges the supports according to their type.
 *   For example: if A1, A2, A3 are the same and A4, A5 are the same,
 *   then the supports are rearranged into A1, A4, A1, A1, A4.
 *   This routine is called in Pre4MV, after the vertex sets are computed.
 *
 * ON ENTRY :
 *   nVar     the number of variables, embient dimension;
 *   Spt      array of points, of range 0..nSpt-1,
 *            with each point having nVar integer coordinates;
 *   SptIdx   array of range 0..nSpt, SptIdx[i] indicates the start
 *            of the i-th support, SptIdx[nVar] contains the total #points;
 *   Vtx      vertex points of the supports;
 *   VtxIdx   indicates the start of the i-th vertex set in Vtx;
 *   SptType  SptType[i] == -1 means the i-th Spt is the base support,
 *            SptType[i] ==  j means the i-th Spt is the j-th Spt;
 *   OldIdx   the position of each point in Vtx
 *            with respect to the original supports in Spt.
 *
 * ON RETURN :
 *   Spt      support sets, rearranged along their type;
 *   SptIdx   rearranged index vector to the support sets;
 *   Vtx      rearranged vertex points for each support;
 *   VtxIdx   rearranged index vector the vertex sets;
 *   SptType  type of each support;
 *   OldIdx   rearranged position index vector of Vtx w.r.t. Spt. */

void NonVertex
 ( int nVar, int nSpt, int *SptIdx, int **Spt, int *ynVtx );
/*
 * DESCRIPTION :
 *   Applies the revised simplex method to decide which points in the
 *   supports are vertices.
 *
 * ON ENTRY :
 *   nVar     the number of variables;
 *   nSpt     the number of supports;
 *   SptIdx   array of range 0..nSpt, SptIdx[i] indicates the start
 *            of the i-th support, SptIdx[nVar] contains the total #points;
 *   Spt      array of points, of range 0..nSpt-1,
 *            with each point having nVar integer coordinates;
 *   ynVtx    array of dimension equal to the total number of points.
 *
 * ON RETURN :
 *   ynVtx    ynVtx[i] == 1 if the i-th point is a vertex,
 *            otherwise the i-th value in the array is zero. */

#endif
