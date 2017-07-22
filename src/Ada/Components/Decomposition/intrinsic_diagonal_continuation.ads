with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Intrinsic_Diagonal_Continuation is

-- DESCRIPTION :
--   A diagonal homotopy allows to compute witness sets for all components
--   of the intersection of two positive dimensional solution sets.
--   The intrinsic version uses a basis for the linear spaces which cut out
--   the witness sets, while the extrinsic version uses the equations.

-- IMPORTANT REQUIREMENT :
--   The polynomial systems and system functions assume the right
--   number of equations, i.e.: complete intersections.
--   Randomize whenever necessary before applying diagonal homotopies.

  function Minimal_Intersection_Dimension
             ( n,a,b : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the minimal dimension of the intersection of 
  --   an a-dimensional with a b-dimensional component in n-space.

  procedure Intrinsic_Diagonal_Homotopy
              ( file : in file_type; report : in boolean; a,b : in natural32;
                p1,p2 : in Poly_Sys; sols1,sols2 : in Solution_List;
                plane1,plane2 : in Matrix;
                sols : out Solution_List; plane : out Matrix );

  -- DESCRIPTION :
  --   Defines the diagonal homotopy in intrinsic coordinates to
  --   intersect an a-dimensional component given by (p1,sols1,plane1)
  --   with a b-dimensional component given by (p2,sols2,plane2).
  --   Output is written to a file.

  -- REQUIRED : a >= b.

  -- ON ENTRY :
  --   file     file for output, if absent, then no output;
  --   report   indicates whether the path trackers will report diagnostics,
  --            if false, then path trackers will be silent;
  --   a        dimension of the 1st component;
  --   b        dimension of the 2nd component;
  --   p1       polynomials for the 1st witness set of dimension a;
  --   p2       polynomials for the 2nd witness set of dimension b;
  --   sols1    intrinsic coordinates for solutions in 1st witness set;
  --   sols2    intrinsic coordinates for solutions in 2nd witness set;
  --   plane1   basis for the affine space to define 1st witness set;
  --   plane2   basis for the affine space to define 2nd witness set.

  -- ON RETURN :
  --   sols     solutions at the last dimensional intersection component;
  --   plane    defines witness set for last dimensional component.

  generic

    with function fA ( x : Vector ) return Vector;
    with function jfA ( x : Vector ) return Matrix;
    with function fB ( x : Vector ) return Vector;
    with function jfB ( x : Vector ) return Matrix;

  procedure Generic_Diagonal_Homotopy
              ( file : in file_type; nefA,nefB,n,a,b : in integer32;
                sols1,sols2 : in Solution_List; plane1,plane2 : in Matrix;
                sols : out Solution_List; plane : out Matrix );

  -- DESCRIPTION :
  --   Generic version of the diagonal homotopy in intrinsic coordinates.

  -- REQUIRED : a >= b.

  -- ON ENTRY :
  --   file     to write diagnostics and results;
  --   nefA     dimension of vectors returned by fA;
  --   nefB     dimension of vectors returned by fB;
  --   n        ambient dimension, x'range = 1..n in fA and fB;
  --   a        dimension of 1st witness set defined by fA,jfA,sols1,plane1;
  --   b        dimension of 2nd witness set defined by fB,jfB,sols2,plane2;
  --   sols1    intrinsic coordinates for solutions in 1st witness set;
  --   sols2    intrinsic coordinates for solutions in 2nd witness set;
  --   plane1   basis for the affine space to define 1st witness set;
  --   plane2   basis for the affine space to define 2nd witness set.

  -- ON RETURN :
  --   sols     solutions at the last dimensional intersection component;
  --   plane    defines witness set for last dimensional component.

end Intrinsic_Diagonal_Continuation;
