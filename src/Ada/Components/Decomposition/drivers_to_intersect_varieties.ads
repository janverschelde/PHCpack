with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Drivers_to_Intersect_Varieties is

-- DESCRIPTION :
--   The drivers in this package are responsible to convert the
--   extrinsic formulation into a complete intersection with
--   system and solutions in intrinsic coordinates.

  procedure Intrinsic_Diagonal_Homotopy
              ( file : in file_type; report : in boolean;
                ep1,ep2 : in Poly_Sys; esols1,esols2 : in Solution_List;
                a,b : in natural32 );

  -- DESCRIPTION :
  --   Executes the diagonal homotopy in intrinsic coordinates to
  --   intersect an a-dimensional component given by (ep1,esols1)
  --   with a b-dimensional component given by (ep2,esols2).
  --   All output is written to a file.

  -- REQUIRED : a >= b.

  -- ON ENTRY :
  --   file     file to write output and diagnostics on;
  --   report   indicates whether path trackers report diagnostics,
  --            if false, then path trackers are silent;
  --   ep1      polynomials for the 1st witness set of dimension a;
  --   ep2      polynomials for the 2nd witness set of dimension b;
  --   esols1   solutions in the a-dimensional witness set;
  --   esols2   solutions in the b-dimensional witness set;
  --   a        dimension of the 1st component;
  --   b        dimension of the 2nd component.

  procedure Intrinsic_Diagonal_Homotopy
              ( file : in file_type; report : in boolean;
                ep1,ep2 : in Poly_Sys; esols1,esols2 : in Solution_List;
                a,b : in natural32; f : out Link_to_Poly_Sys;
                p : out Link_to_Matrix; s : out Solution_List );

  -- DESCRIPTION :
  --   Executes the diagonal homotopy in intrinsic coordinates to
  --   intersect an a-dimensional component given by (ep1,esols1)
  --   with a b-dimensional component given by (ep2,esols2).
  --   All output is written to a file.

  -- REQUIRED : a >= b.

  -- ON ENTRY :
  --   file     file to write output and diagnostics on;
  --   report   indicates whether path trackers report diagnostics,
  --            if false, then path trackers are silent;
  --   ep1      polynomials for the 1st witness set of dimension a;
  --   ep2      polynomials for the 2nd witness set of dimension b;
  --   esols1   solutions in the a-dimensional witness set;
  --   esols2   solutions in the b-dimensional witness set;
  --   a        dimension of the 1st component;
  --   b        dimension of the 2nd component.

  -- ON RETURN :
  --   f        polynomial system stacked with ep1 and ep2,
  --            stripped from their original embedding;
  --   p        generators of the intersecting plane cutting out s;
  --   s        intrinsic coordinates of lowest dimensional piece
  --            of the intersection of the two components.

  generic

    with function fA ( x : Vector ) return Vector;
    with function jfA ( x : Vector ) return Matrix;
    with function fB ( x : Vector ) return Vector;
    with function jfB ( x : Vector ) return Matrix;

  procedure Generic_Diagonal_Homotopy
              ( file : in file_type; nefA,nefB,n,a,b : in natural32;
                esols1,esols2 : in Solution_List; p1,p2 : in VecVec );

  -- DESCRIPTION :
  --   Generic version of the diagonal homotopy in intrinsic coordinates.

  -- REQUIRED : a >= b.

  -- ON ENTRY :
  --   file     to write diagnostics and results;
  --   nefA     dimension of vectors returned by fA;
  --   nefB     dimension of vectors returned by fB;
  --   n        ambient dimension, x'range = 1..n in fA and fB;
  --   a        dimension of 1st witness set defined by fA,jfA,esols1,p1;
  --   b        dimension of 2nd witness set defined by fB,jfB,esols2,p2;
  --   esols1   solutions in a-dimensional 1st witness set;
  --   esols2   solutions in b-dimensional 2nd witness set;
  --   p1       linear equations to cut out esols1, p1'range = 1..a;
  --   p2       linear equations to cut out esols2, p2'range = 1..b.

  generic

    with function fA ( x : Vector ) return Vector;
    with function jfA ( x : Vector ) return Matrix;

  procedure B_Call_Generic_Diagonal_Homotopy
              ( file : in file_type; nef1,n,a,b : in natural32;
                ep2 : in Poly_Sys; esols1,esols2 : in Solution_List;
                s1 : in VecVec );

  -- DESRIPTION :
  --   Given the functions for component A in complete intersection,
  --   this generic procedures prepares the functions for component B
  --   and instantiates the generic diagonal homotopy.

  procedure A_Call_Generic_Diagonal_Homotopy
              ( file : in file_type; ep1,ep2 : in Poly_Sys;
                esols1,esols2 : in Solution_List; a,b : in natural32 );

  -- DESCRIPTION :
  --   Prepares the functions for the first component A and calls
  --   then the procedure to prepare the functions for component B.

end Drivers_to_Intersect_Varieties;
