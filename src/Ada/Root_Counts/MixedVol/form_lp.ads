with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;
with Standard_Floating_VecMats;
with Index_Tree_LP;                     use Index_Tree_LP;
with Relation_Table;                    use Relation_Table;

package Form_LP is

-- DESCRIPTION :
--    The routines in this package define the main LP problems,
--    and are called in the main program "mixed_volume".

  procedure Form_LP
                ( nVar,nSpt : in integer32;
                  SptType,SptIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  RelTab : in Boolean_Matrix;
                  ElmMtrx : in Standard_Floating_Matrices.Link_to_Matrix;
                  NonZero : in out Standard_Integer_Vectors.Link_to_Vector;
                  Lvl2CoDim : in Standard_Integer_Vectors.Link_to_Vector;
                  A : in Standard_Floating_VecMats.Link_to_VecMat;
                  Lvl : in integer32; ptr : in Link_to_LPdata;
                  Cell : in Standard_Integer_Vectors.Link_to_Vector;
                  Cell_Orig : in out Standard_Integer_Vectors.Link_to_Vector;
                  Lvl2LPdim : in Standard_Integer_Vectors.Link_to_Vector;
                  Lvl2Spt : in Standard_Integer_Vectors.Link_to_Vector;
                  FixFrstPt : in Standard_Integer_Vectors.Link_to_Vector;
                  FixLstPt : in Standard_Integer_Vectors.Link_to_Vector;
                  MinNumPt : in Standard_Integer_Vectors.Link_to_Vector;
                  PtIn : in out Standard_Integer_VecVecs.Link_to_VecVec;
                  Strt1Pt,End1Pt : out integer32;
                  ToOrig : in out Standard_Integer_VecVecs.Link_to_VecVec;
                  Pre2Cur : in out Standard_Integer_VecVecs.Link_to_VecVec; 
                  LPdim : out integer32;
                  FrstPt : in out Standard_Integer_Vectors.Link_to_Vector;
                  LstPt : in out Standard_Integer_Vectors.Link_to_Vector;
                  Cur2Pre : in out Standard_Integer_Vectors.Link_to_Vector;
                  x : in out Standard_Floating_Vectors.Link_to_Vector;
                  Binv : in out Standard_Floating_Matrices.Link_to_Matrix;
                  Bidx : in out Standard_Integer_Vectors.Link_to_Vector;
                  info : out integer32 );

  -- DESCRIPTION :
  --   Defines the LP problems when the value of Lvl > 1.
  --
  -- ON ENTRY :
  --   nVar       number of variables, ambient dimension;
  --   nSpt       number of different supports;
  --   SptType    array of range 0..nSpt with #occurrences of each support;
  --   SptIdx     SptIdx[i] indicates the start of the i-th support;
  --   RelTab     the relation table records edges on lower hull of
  --              each support or on pairwise Minkowski sums;
  --   ElmMtrx    elimination matrix has dimension the size of a cell;
  --   NonZero    array of dimension equal to the size of a mixed cell,
  --              but only one entry will be updated;
  --   Lvl2CoDim  array of dimension equal to the size of a mixed cell,
  --              counts the number of variables eliminated,
  --              only the element at position Lvl is used;
  --   A          3-dimensional matrix, first dimension is the cell size,
  --              then every A[i] is a coefficient matrix of constraints,
  --              of dimension #constraints times nVar + 1;
  --   Lvl        index of the current level;
  --   ptr        pointer to LP data, used only for consultation;
  --   Cell       indices of the current mixed cell;
  --   Cell_Orig  indices to the original mixed cell;
  --   Lvl2LPdim  array of cell size, only entries at Lvl and Lvl-1 are used;
  --   Lvl2Spt    array of cell size, only entry at Lvl is used;
  --   FixFrstPt  array of cell size, only entry at Lvl is used;
  --   FixLstPt   array of cell size, only entry at Lvl is used;
  --   MinNumPt   array of cell size, only entry at Lvl is used;
  --   PtIn       2-dimensional matrix with as many rows as the cell size
  --              and as many columns as the total number of points,
  --              only the points in PtIn[Lvl] and PtIn[Lvl-1] are accessed;
  --   ToOrig     2-dimensional integer matrix, as many rows as cell size,
  --              and as many columns as the total number of points;
  --   Pre2Cur    2-dimensional integer matrix, as many rows as cell size,
  --              and as many columns as the total number of points;
  --   FrstPt     array of size of a mixed cell,
  --              serves as index vector to first point in current support;
  --   LstPt      array of size of mixed cell,
  --              serves as index vector to last point in current support;
  --   Cur2Pre    array of dimension equal to total number of points;
  --   x          solution vector of dimension nVar;
  --   Binv       matrix of dimension nVar, basis inverse;
  --   Bidx       index vector for points in Pre2Cur, of dimension nVar.
  --
  -- ON RETURN :
  --   ElmMtrx    updated elimination matrix;
  --   NonZero    entry at Lvl2CoDim[Lvl] is updated;
  --   Cell_Orig  updated indices of the original mixed cell;
  --   PtIn       updated entries at row Lvl;
  --   Strt1Pt    index for where to start the 1-point tests;
  --   End1Pt     index for where to end the 1-point tests;
  --   ToOrig     updated entries at rows Lvl-1 and Lvl;
  --   Pre2Cur    updated entries at row 1 and row Lvl;
  --   LPdim      equals Lvl2LPdim[Lvl];
  --   Cur2Pre    updated index vector;
  --   x          updated solution vector;
  --   Binv       updated basis inverse;
  --   Bidx       updated index vector;
  --   info       0 if normal return,
  --              1 if not enough points to extend, must backtrack.

  procedure Form_LP1
                ( nVar,nSpt : in integer32;
                  SptType,SptIdx : in Standard_Integer_Vectors.Link_to_Vector;
                  RelTab : in Boolean_Matrix;
                  A : in Standard_Floating_VecMats.Link_to_VecMat;
                  Lvl : in integer32;
                  Cell : in Standard_Integer_Vectors.Link_to_Vector;
                  Lvl2LPdim : in Standard_Integer_Vectors.Link_to_Vector;
                  FixLstPt : in Standard_Integer_Vectors.Link_to_Vector;
                  MinNumPt : in Standard_Integer_Vectors.Link_to_Vector;
                  PtIn : in out Standard_Integer_VecVecs.Link_to_VecVec;
                  ToOrig : in out Standard_Integer_VecVecs.Link_to_VecVec;
                  Pre2Cur : in out Standard_Integer_VecVecs.Link_to_VecVec;
                  FrstPt : in out Standard_Integer_Vectors.Link_to_Vector;
                  LstPt : in out Standard_Integer_Vectors.Link_to_Vector;
                  info : out integer32 );

  -- DESCRIPTION :
  --   Defines the LP problems when the value of Lvl <= 1.
  --
  -- ON ENTRY :
  --   nVar       number of variables, ambient dimension;
  --   nSpt       number of different supports;
  --   SptType    array of range 0..nSpt with #occurrences of each support;
  --   SptIdx     SptIdx[i] indicates the start of the i-th support;
  --   RelTab     the relation table records edges on lower hull of
  --              each support or on pairwise Minkowski sums;
  --   A          3-dimensional matrix, first dimension is the cell size,
  --              then every A[i] is a coefficient matrix of constraints,
  --              of dimension #constraints times nVar + 1;
  --   Lvl        index of the current level;
  --   Cell       indices of the current mixed cell;
  --   Lvl2LPdim  array of cell size, only entry at Lvl-1 is used;
  --   FixLstPt   array of cell size, only entry at Lvl is used;
  --   MinNumPt   array of cell size, only entry at Lvl is used;
  --   PtIn       2-dimensional matrix with as many rows as the cell size
  --              and as many columns as the total number of points,
  --              only the points in PtIn[Lvl] are accessed;
  --   ToOrig     2-dimensional integer matrix, as many rows as cell size,
  --              and as many columns as the total number of points;
  --   Pre2Cur    2-dimensional integer matrix, as many rows as cell size,
  --              and as many columns as the total number of points;
  --   FrstPt     array of size of a mixed cell,
  --              serves as index vector to first point in current support;
  --   LstPt      array of size of mixed cell,
  --              serves as index vector to last point in current support.
  --
  -- ON RETURN :
  --   A          constraints updated in A[Lvl], using A[Lvl-1];
  --   PtIn       updated entries in PtIn[Lvl];
  --   ToOrig     updated entries in ToOrig[Lvl];
  --   Pre2Cur    updated entries in Pre2Cur[Lvl];
  --   FrstPt     set FrstPt[Lvl] to zero;
  --   LstPt      set LstPt[Lvl] to #constraints involved;
  --   info       0 if normal return,
  --              1 if not enough points to extend, must backtrack.
  
end Form_LP;
