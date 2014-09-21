with Quad_Double_Matrices;
with QuadDobl_Complex_Matrices;

package QuadDobl_Matrix_Inversion is

-- DESCRIPTION :
--   The functions below return the inverse of a given square matrix,
--   using quad double arithmetic.
--   This type of operation is useful for change of basis.

  function Inverse ( m : Quad_Double_Matrices.Matrix )
                   return Quad_Double_Matrices.Matrix;

  function Inverse ( m : QuadDobl_Complex_Matrices.Matrix )
                   return QuadDobl_Complex_Matrices.Matrix;

end QuadDobl_Matrix_Inversion;
