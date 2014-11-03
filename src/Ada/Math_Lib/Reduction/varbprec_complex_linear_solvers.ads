with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Matrices;

package Varbprec_Complex_Linear_Solvers is

-- DESCRIPTION :
--   The 'varbprec' is a contraction of variable precision.
--   Based on estimates for the condition number of the coefficient matrix
--   of the linear system, we can bound the expected loss of correct
--   decimal places in the computed solution.  If the expected loss is too
--   large to meet the wanted number of correct decimal places, then the
--   solving of the linear system should be done in a higher precision.
--   The condition number estimators are based on the LINPACK library.

  procedure Estimated_Loss_of_Decimal_Places
              ( mtx : in out Standard_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out double_float; loss : out integer32 );

  -- DESCRIPTION :
  --   Returns the estimated loss of decimal places that may occur
  --   when solving a linear system with coefficient matrix mat.
  --   As the estimate is computed in fixed double precision,
  --   the result is not reliable for matrices with condition numbers
  --   that are larger than 1.0E+15.

  -- REQUIRED : mat'range(1) = mat'range(2).

  -- ON ENTRY :
  --   mtx      matrix of floating-point complex numbers.

  -- ON RETURN :
  --   mtx      output of lufco, suitable for backsubstitution
  --            if nonsingular;
  --   piv      pivoting information computed by lufco;
  --   rco      estimate for the inverse of the condition number,
  --            as computed by lufco on the matrix;
  --   loss     logarithm of rco, indicates the loss of decimal
  --            places when solving a linear system with matrix mtx.

  function Estimated_Loss_of_Decimal_Places
              ( mtx : Standard_Complex_Matrices.Matrix )
              return integer32;

  -- DESCRIPTION :
  --   Returns the estimated loss of decimal places that may occur
  --   when solving a linear system with coefficient matrix mtx.
  --   As the estimate is computed in fixed double precision,
  --   the result is not reliable for matrices with condition numbers
  --   that are larger than 1.0E+15.

  -- REQUIRED : mtx'range(1) = mtx'range(2).

end Varbprec_Complex_Linear_Solvers;
