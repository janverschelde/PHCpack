with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Vectors;
with Standard_Natural_Matrices;
with Standard_Floating_Matrices;
with Standard_Complex_Matrices;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;

package SAGBI_Homotopies is

-- DESCRIPTION :
--   Provides basic routines to set up the polynomials in the SAGBI homotopy.
--   Due to hexadecimal expansions, n and d are both limited to 16.
--   In practice, since the #roots grow so rapidly, this is no limitation.

  function Lifted_Localized_Laplace_Expansion ( n,d : natural32 ) return Poly;

  -- DESCRIPTION :
  --   Constructs the generic equation in the SAGBI homotopy.
  --   The coefficients are brackets in hexadecimal expansion.
  --   These brackets represents the selected rows for the maximal minors.
  --   The localization chooses the lower-right upper block of the d-plane
  --   as the identity matrix.

  function Lifted_Localized_Laplace_Expansion
             ( locmap : Standard_Natural_Matrices.Matrix ) return Poly;

  -- DESCRIPTION :
  --   The generic equation in the SAGBI homotopy is constructed using
  --   the localizaton map in locmap.  Zeros and ones indicate the
  --   position of the identity matrix, while free elements are twos.

  function Intersection_Coefficients
              ( m : Standard_Floating_Matrices.Matrix;
                c : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Vectors.Vector;
  function Intersection_Coefficients
              ( m : Standard_Complex_Matrices.Matrix;
                c : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Given a matrix m and the hexadecimal expansion of the coefficients
  --   in c, the vector of maximal minors of m is returned.

  function Intersection_Condition
              ( m : Standard_Floating_Matrices.Matrix; p : Poly ) return Poly;
  function Intersection_Condition
              ( m : Standard_Complex_Matrices.Matrix; p : Poly ) return Poly;

  -- DESCRIPTION :
  --   Generates the particular equations in the SAGBI homotopy, with
  --   as input a matrix m and the lifted localized Laplace expansion p.
  --   The matrix contains in its columns the generating points of the
  --   plane of intersection.

  -- REQUIRED : The dimensions of the matrix m are n times n-d.

end SAGBI_Homotopies;
