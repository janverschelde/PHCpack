with Standard_Floating_Matrices;
with Standard_Complex_Matrices;

package Standard_Matrix_Inversion is

-- DESCRIPTION :
--   The functions below return the inverse of a given square matrix.
--   This type of operation is useful for change of basis.

  function Inverse ( m : Standard_Floating_Matrices.Matrix )
                   return Standard_Floating_Matrices.Matrix;

  function Inverse ( m : Standard_Complex_Matrices.Matrix )
                   return Standard_Complex_Matrices.Matrix;

end Standard_Matrix_Inversion;
