with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Matrices;          use DoblDobl_Complex_Matrices;

package DoblDobl_Complex_Matrix_Norms is

-- DESCRIPTION :
--   Defines the max norm for complex matrices of double double precision.

  function Max_Norm ( m : Matrix ) return double_double;

  -- DESCRIPTION :
  --   Returns the absolute value of the largest element in m.

end DoblDobl_Complex_Matrix_Norms;
