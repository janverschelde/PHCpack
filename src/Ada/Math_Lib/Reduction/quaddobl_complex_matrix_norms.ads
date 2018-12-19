with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Matrices;          use QuadDobl_Complex_Matrices;

package QuadDobl_Complex_Matrix_Norms is

-- DESCRIPTION :
--   Defines the max norm for complex matrices of quad double precision.

  function Max_Norm ( m : Matrix ) return quad_double;

  -- DESCRIPTION :
  --   Returns the absolute value of the largest element in m.

end QuadDobl_Complex_Matrix_Norms;
