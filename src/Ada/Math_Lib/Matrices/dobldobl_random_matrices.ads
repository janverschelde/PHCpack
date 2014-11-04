with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Double_Double_Matrices;
with DoblDobl_Complex_Matrices;

package DoblDobl_Random_Matrices is

-- DESCRIPTION :
--   Offers routines to generate matrices of random double doubles.

  function Random_Matrix ( n,m : natural32 )
                         return Double_Double_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m with random double doubles.

  function Orthogonalize ( mat : Double_Double_Matrices.Matrix )
                         return Double_Double_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the orthogonal matrix with the same column span as mat.

  function Random_Orthogonal_Matrix
                ( n,m : natural32 ) return Double_Double_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m with random floating numbers,
  --   where the columns are orthogonal w.r.t. each other.

  function Random_Matrix ( n,m : natural32 )
                         return DoblDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix of range 1..n,1..m
  --   with random complex double double numbers.

  function Orthogonalize ( mat : DoblDobl_Complex_Matrices.Matrix )
                         return DoblDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the orthogonal matrix with the same column span as mat.

  function Random_Orthogonal_Matrix
             ( n,m : natural32 ) return DoblDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a random matrix where the columns form an orthonormal basis.

end DoblDobl_Random_Matrices;
