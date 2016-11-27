with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Complex_Numbers;         use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Dense_Vector_Series;
with Standard_Dense_Matrix_Series;

package Standard_Interpolating_Series is

-- DESCRIPTION :
--   Via interpolation at random points we can solve linear systems
--   with matrix series, even if several of the leading coefficients
--   of the matrix coefficients are singular matrices.

  function Eval ( v : Standard_Dense_Vector_Series.Vector;
                  t : Complex_Number )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector obtained by evaluating v at t.

  function Eval ( m : Standard_Dense_Matrix_Series.Matrix;
                  t : Complex_Number )
                return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the matrix obtained by evaluating m at t.

  function Full_Rank
             ( m : Standard_Dense_Matrix_Series.Matrix;
               d : integer32; verbose : boolean := true ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the matrix series up to degree d <= m.deg
  --   has full rank when evaluated at a random complex value.

  function Full_Rank
             ( m : Standard_Dense_Matrix_Series.Matrix;
               verbose : boolean := true ) return integer32;

  -- DESCRIPTION :
  --   Returns -1 if for all degrees d in 0..m.deg the matrix has no
  --   full rank when evaluated at a random complex value, otherwise
  --   returns the smallest d for which the evaluation of m at a
  --   random complex value yields a full rank matrix.

end Standard_Interpolating_Series;
