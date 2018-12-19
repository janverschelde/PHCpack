with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;

package Standard_Complex_Matrix_Norms is

-- DESCRIPTION :
--   Defines the max norm for complex matrices of standard double precision.

  function Max_Norm ( m : Matrix ) return double_float;

  -- DESCRIPTION :
  --   Returns the absolute value of the largest element in m.

end Standard_Complex_Matrix_Norms;
