with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Double_Double_Numbers;            use Double_Double_Numbers;
with Quad_Double_Numbers;              use Quad_Double_Numbers;
with Multprec_Floating_Numbers;        use Multprec_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Double_Double_Vectors;
with Quad_Double_Vectors;
with Multprec_Floating_Vectors;
with Standard_Floating_Matrices;
with Double_Double_Matrices;
with Quad_Double_Matrices;
with Multprec_Floating_Matrices;

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
  procedure Estimated_Loss_of_Decimal_Places
              ( mtx : in out Multprec_Floating_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out Floating_Number; loss : out integer32 );

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
  function Estimated_Loss_of_Decimal_Places
              ( mtx : Multprec_Floating_Matrices.Matrix )
              return integer32;

  -- DESCRIPTION :
  --   Returns the estimated loss of decimal places that may occur
  --   when solving a linear system with coefficient matrix mtx.
  --   This function encapsulates the above procedure with the same name.

  -- REQUIRED : mtx'range(1) = mtx'range(2).

  procedure Solve_to_Wanted_Decimal_Places
              ( mtx : in out Standard_Floating_Matrices.Matrix;
                rhs : in out Standard_Floating_Vectors.Vector;
                want : in integer32; fail : out boolean;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out double_float; loss : out integer32 );
  procedure Solve_to_Wanted_Decimal_Places
              ( mtx : in out Double_Double_Matrices.Matrix;
                rhs : in out Double_Double_Vectors.Vector;
                want : in integer32; fail : out boolean;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out double_double; loss : out integer32 );
  procedure Solve_to_Wanted_Decimal_Places
              ( mtx : in out Quad_Double_Matrices.Matrix;
                rhs : in out Quad_Double_Vectors.Vector;
                want : in integer32; fail : out boolean;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out quad_double; loss : out integer32 );

  -- DESCRIPTION :
  --   Estimates the loss of decimal places based on the estimated
  --   condition number, via the above procedures Estimated_Loss_.

  -- REQUIRED : mat'range(1) = mat'range(2).

  -- ON ENTRY :
  --   mtx      matrix of floating-point complex numbers;
  --   rhs      righthandside vector of a linear system;
  --   want     wanted number of correct decimal places.

  -- ON RETURN :
  --   mtx      output of lufco, suitable for backsubstitution
  --            if nonsingular;
  --   rhs      if not fail, then rhs contains the solution to
  --            the linear system with coefficient matrix mtx
  --            and as righthandside vector the given rhs,
  --            correct with as many decimal places as want;
  --   fail     if true, then the precision did not suffice
  --            to meet the wanted number of decimal places,
  --            otherwise, the solution in rhs is accurate enough;
  --   piv      pivoting information computed by lufco;
  --   rco      estimate for the inverse of the condition number,
  --            as computed by lufco on the matrix;
  --   loss     logarithm of rco, indicates the loss of decimal
  --            places when solving a linear system with matrix mtx.

  procedure Solve_to_Wanted_Decimal_Places
              ( mtx : in out Multprec_Floating_Matrices.Matrix;
                rhs : in out Multprec_Floating_Vectors.Vector;
                want : in integer32;
                piv : out Standard_Integer_Vectors.Vector;
                rco : out Floating_Number; loss : out integer32 );

  -- DESCRIPTION :
  --   Estimates the loss of decimal places based on the estimated
  --   condition number and raises the working precision to meet the
  --   wanted number of decimal places, taking into account the expected
  --   loss of decimal places, due to the conditioning of the matrix.
  --   Unlike the other solvers, we can always raise the precision to
  --   meet the desired accuracy, so failure is no option.

  -- REQUIRED : mat'range(1) = mat'range(2).

  -- ON ENTRY :
  --   mtx      matrix of floating-point complex numbers;
  --   rhs      righthandside vector of a linear system;
  --   want     wanted number of correct decimal places.

  -- ON RETURN :
  --   mtx      output of lufco, suitable for backsubstitution
  --            if nonsingular;
  --   rhs      if not fail, then rhs contains the solution to
  --            the linear system with coefficient matrix mtx
  --            and as righthandside vector the given rhs,
  --            correct with as many decimal places as want;
  --   piv      pivoting information computed by lufco;
  --   rco      estimate for the inverse of the condition number,
  --            as computed by lufco on the matrix;
  --   loss     logarithm of rco, indicates the loss of decimal
  --            places when solving a linear system with matrix mtx.

end Varbprec_Floating_Linear_Solvers;
