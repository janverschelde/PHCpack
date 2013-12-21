with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Natural_Matrices;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;            use Standard_Complex_VecMats;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Solutions;          use Standard_Complex_Solutions;
with Localization_Posets;                 use Localization_Posets;

package Verification_with_Determinants is

-- DESCRIPTION :
--   This package contains facility to verify whether the curves computed by
--   the Pieri homotopies satisfy all interpolation-intersection conditions.

  function Determinant
              ( m1,m2 : Standard_Complex_Matrices.Matrix )
              return Complex_Number;

  -- DESCRIPTION :
  --   Returns the determinant of the matrix [m1|m2].

  -- REQUIRED :
  --   number of rows of m1 and m2 are both equal to m+p.

  -- ON ENTRY :
  --   m1      contains in its columns the generators of a p-plane;
  --   m2      contains in its columns the generators of an m-plane.

  -- ON RETURN :
  --   the determinant of the matrix with the columns of m1 and the
  --   columns m2 is zero when the planes intersect each other.

  procedure Determinant
              ( file : in file_type;
                m1,m2 : in Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   The matrix m1 contains the generators for a p-plane,
  --   the matrix m2 contains the generators for an m-plane.
  --   If the matrices both have m+p rows, then the determinant is
  --   computed and written to file.

  procedure Verify_Determinants
              ( file : in file_type; dim : integer32; nd : in Node;
                xpm : in Standard_Complex_Poly_Matrices.Matrix;
                locmap : in Standard_Natural_Matrices.Matrix;
                locsols : in Standard_Complex_VecVecs.VecVec;
                s : in Standard_Complex_Vectors.Vector; ip : in VecMat );
  procedure Verify_Determinants
              ( file : in file_type; nd : in Node;
                xpm : in Standard_Complex_Poly_Matrices.Matrix;
                locmap : in Standard_Natural_Matrices.Matrix;
                locsols : in Solution_List;
                s : in Standard_Complex_Vectors.Vector; ip : in VecMat );

  -- DESCRIPTION :
  --   Verifies the interpolation-intersection conditions by determinants.

  -- ON ENTRY :
  --   file     on this file all determinants will be written;
  --   dim      dimension of the solution vectors;
  --   nd       contains top and bottom pivots;
  --   xpm      curve in symbolic form, before the localization;
  --   locmap   localization pattern;
  --   locsols  solutions according to the localization pattern;
  --   s        values for the interpolation points;
  --   ip       m-planes that need to be met at the interpolation points.

end Verification_with_Determinants;
