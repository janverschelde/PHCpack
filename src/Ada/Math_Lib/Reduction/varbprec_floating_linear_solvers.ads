with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Double_Double_Numbers;            use Double_Double_Numbers;
with Quad_Double_Numbers;              use Quad_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Matrices;
with Double_Double_Matrices;
with Quad_Double_Matrices;

package Varbprec_Floating_Linear_Solvers is

-- DESCRIPTION :
--   The 'varbprec' is a contraction of variable precision.
--   Based on estimates for the condition number of the coefficient matrix
--   of the linear system, we can bound the expected loss of correct
--   decimal places in the computed solution.  If the expected loss is too
--   large to meet the wanted number of correct decimal places, then the
--   solving of the linear system should be done in a higher precision.
--   The condition number estimators are based on the LINPACK library.

  procedure Estimated_Loss_of_Decimal_Places
              ( mtx : in out Standard_Floating_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out double_float; loss : out integer32 );
  procedure Estimated_Loss_of_Decimal_Places
              ( mtx : in out Double_Double_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out double_double; loss : out integer32 );
  procedure Estimated_Loss_of_Decimal_Places
              ( mtx : in out Quad_Double_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out quad_double; loss : out integer32 );

  -- DESCRIPTION :
  --   Returns the estimated loss of decimal places that may occur
  --   when solving a linear system with coefficient matrix mat.
  --   If the estimates are computed respectively in fixed double 
  --   double double, and quad double precision,
  --   the results are not reliable for matrices with condition numbers
  --   that are respectively larger than 1.0E+16, 1.0E+32, and 1.0E+64.

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
              ( mtx : Standard_Floating_Matrices.Matrix )
              return integer32;
  function Estimated_Loss_of_Decimal_Places
              ( mtx : Double_Double_Matrices.Matrix )
              return integer32;
  function Estimated_Loss_of_Decimal_Places
              ( mtx : Quad_Double_Matrices.Matrix )
              return integer32;

  -- DESCRIPTION :
  --   Returns the estimated loss of decimal places that may occur
  --   when solving a linear system with coefficient matrix mtx.
  --   This function encapsulates the above procedure with the same name.

  -- REQUIRED : mtx'range(1) = mtx'range(2).

end Varbprec_Floating_Linear_Solvers;
