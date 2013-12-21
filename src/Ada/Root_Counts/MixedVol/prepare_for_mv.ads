with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;

package Prepare_for_MV is

-- DESCRIPTION :
--   This packages defines the specifications of the routines used
--   to preprocess a tuple of supports for the mixed-volume computation.

  procedure Pre4MV
               ( nVar,nSpt : in integer32; nS : out integer32;
                 SptType : in out Standard_Integer_Vectors.Link_to_Vector;
                 Spt : in out Standard_Integer_VecVecs.Link_to_VecVec;
                 SptIdx : in out Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : in out Standard_Integer_VecVecs.Link_to_VecVec;
                 VtxIdx : in out Standard_Integer_Vectors.Link_to_Vector;
                 OldIdx : in out Standard_Integer_Vectors.Link_to_Vector;
                 perm : out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   As a preprocessing to the mixed-volume computation, this function
  --   extracts the non-vertex points and determines the type of mixture.
  --
  -- ON ENTRY :
  --   nVar      the number of variables, ambient dimension;
  --   nSpt      the number of supports, equals nVar;
  --   SptType   allocated memory for an array of nVar integers;
  --   Spt       array of points, of range 0..nSpt-1,
  --             with each point having nVar integer coordinates;
  --   SptIdx    array of range 0..nSpt, SptIdx[i] indicates the start
  --             of the i-th support, SptIdx[nVar] contains the total #points;
  --   Vtx       memory allocated for as many points as in Spt;
  --   VtxIdx    copy of SptIdx;
  --   OldIdx    memory allocated for total number of points.
  --
  -- ON RETURN :
  --   nS        the number of different supports;
  --   SptType   type of each support, defined for range 0..nS-1;
  --   Vtx       vertex points for each support;
  --   VtxIdx    array of range 0..nS, VtxIdx[i] indicates the start of the
  --             i-th support in Vtx, VtxIdx[nVar] equals total #points;
  --   OldIdx    the position of each point in Vtx
  --             with respect to the original supports in Spt;
  --   perm      vector of range 0..nVar-1 with permutation of the sets
  --             in Spt, making smallest support set to go first.

  procedure SortSpt
               ( nVar,nSpt : in integer32;
                 Spt : in out Standard_Integer_VecVecs.Link_to_VecVec;
		 SptIdx : in Standard_Integer_Vectors.Link_to_Vector );

  -- DECRIPTION :
  --   Sorts the point in each support in descending lexicographic order.
  --
  -- ON ENTRY :
  --   nVar      the ambient dimension;
  --   nSpt      number of supports;
  --   Spt       the coordinates of the supports are in the rows;
  --   SptIdx    is the index of the 1st point of each support in Spt,
  --             SptIdx[nSpt]-1 is the index of the last point of the supports.
  --
  -- ON RETURN :
  --   Spt       sorted supports.

  procedure LowerTriangular
               ( A : in out Standard_Floating_Matrices.Link_to_Matrix;
                 jStrt,jEnd,na : in integer32; rnk : out integer32;
                 ib : in out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Makes the matrix A lower triangular, called by NonVertex below.
  --
  -- ON ENTRY :
  --   A         matrix to be triangulated;
  --   jStrt     index of start row;
  --   jEnd      index of end row;
  --   na        number of columns in A.
  -- 
  -- ON RETURN :
  --   A         lower triangular matrix;
  --   rnk       rank of the matrix;
  --   ib        pivoting information.

  procedure RSimplex
               ( A : in Standard_Floating_Matrices.Link_to_Matrix;
                 b : in Standard_Floating_Vectors.Link_to_Vector;
                 jStrt,jEnd,m : in integer32;
                 xb,sb : in out Standard_Floating_Vectors.Link_to_Vector;
                 ib : in out Standard_Integer_Vectors.Link_to_Vector;
                 binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                 info : out integer32 );

  -- DESCRIPTION :
  --   Applies the revised simplex method to solve
  --
  --     min e   subject to   -e <= 0 and -e + B y <= b
  --
  --   where (e,y) are variables.
  --   The constraint -e <= 0 corresponds to the jStrt-th row of A.
  --
  -- ON ENTRY :
  --   A         active constraints range from row jStrt to jEnd;      
  --   b         right hand side vector;
  --   jStrt     index to the row in A marking the start of the active part;
  --   JEnd      index to the column in A marking the end of the active part;
  --   m         columns in A go up to m;
  --   xb        initial feasible solution;
  --   sb        workspace, initially declared as global variable,
  --             managed by the function NonVertex;
  --   ib        indices of the basis;
  --   binv      inverse of the initial base matrix,
  --             initialized to the identity matrix.
  --
  -- ON RETURN :
  --   xb        updated basis vector;
  --   sb        updated workspace vector;
  --   ib        updated index vector;
  --   binv      updated basis inverse matrix;
  --   info      position of the variable epsilon when info>=0.

  procedure ReArrangeSpt
               ( nVar : in integer32;
                 Spt : in out Standard_Integer_VecVecs.Link_to_VecVec;
                 SptIdx : in out Standard_Integer_Vectors.Link_to_Vector;
                 Vtx : in out Standard_Integer_VecVecs.Link_to_VecVec;
                 VtxIdx : in out Standard_Integer_Vectors.Link_to_Vector;
                 SptType : in out Standard_Integer_Vectors.Link_to_Vector;
	         OldIdx : in out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Rearranges the supports according to their type.
  --   For example: if A1, A2, A3 are the same and A4, A5 are the same,
  --   then the supports are rearranged into A1, A4, A1, A1, A4.
  --   This routine is called in Pre4MV, after the vertex sets are computed.
  --
  -- ON ENTRY :
  --   nVar      the number of variables, embient dimension;
  --   Spt       array of points, of range 0..nSpt-1,
  --             with each point having nVar integer coordinates;
  --   SptIdx    array of range 0..nSpt, SptIdx[i] indicates the start
  --             of the i-th support, SptIdx[nVar] contains the total #points;
  --   Vtx       vertex points of the supports;
  --   VtxIdx    indicates the start of the i-th vertex set in Vtx;
  --   SptType   SptType[i] == -1 means the i-th Spt is the base support,
  --             SptType[i] ==  j means the i-th Spt is the j-th Spt;
  --   OldIdx    the position of each point in Vtx
  --             with respect to the original supports in Spt.
  --
  -- ON RETURN :
  --   Spt       support sets, rearranged along their type;
  --   SptIdx    rearranged index vector to the support sets;
  --   Vtx       rearranged vertex points for each support;
  --   VtxIdx    rearranged index vector the vertex sets;
  --   SptType   type of each support;
  --   OldIdx    rearranged position index vector of Vtx w.r.t. Spt.

  procedure NonVertex
               ( nVar,nSpt : in integer32; 
                 SptIdx : in Standard_Integer_Vectors.Link_to_Vector;
                 Spt : in Standard_Integer_VecVecs.Link_to_VecVec;
                 ynVtx : in out Standard_Integer_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Applies the revised simplex method to decide which points in the
  --   supports are vertices.
  -- 
  -- ON ENTRY :
  --   nVar      the number of variables;
  --   nSpt      the number of supports;
  --   SptIdx    array of range 0..nSpt, SptIdx[i] indicates the start
  --             of the i-th support, SptIdx[nVar] contains the total #points;
  --   Spt       array of points, of range 0..nSpt-1,
  --             with each point having nVar integer coordinates;
  --   ynVtx     array of dimension equal to the total number of points.
  --
  -- ON RETURN :
  --   ynVtx     ynVtx[i] == 1 if the i-th point is a vertex,
  --             otherwise the i-th value in the array is zero.

end Prepare_for_MV;
