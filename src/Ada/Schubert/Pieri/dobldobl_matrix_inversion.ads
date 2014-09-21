with Double_Double_Matrices;
with DoblDobl_Complex_Matrices;

package DoblDobl_Matrix_Inversion is

-- DESCRIPTION :
--   The functions below return the inverse of a given square matrix,
--   using double double arithmetic.
--   This type of operation is useful for change of basis.

  function Inverse ( m : Double_Double_Matrices.Matrix )
                   return Double_Double_Matrices.Matrix;

  function Inverse ( m : DoblDobl_Complex_Matrices.Matrix )
                   return DoblDobl_Complex_Matrices.Matrix;

end DoblDobl_Matrix_Inversion;
