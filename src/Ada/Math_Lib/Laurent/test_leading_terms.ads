with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Boolean_Vectors;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Floating_Matrices;
with Standard_Complex_Matrices;

package Test_Leading_Terms is

-- DESCRIPTION :
--   Tests the computation of the leading terms of a power series
--   c(x)*[t^x] where both the * and ^ are componentwise operations:
--   (1) x is a vector of floating-point coefficients, and
--   (2) c(x) is a complex vector, all sufficiently generic,
--   as the solution of a linear system of series with real powers,
--   of the form C(A)*[t^A]*c(x)*[t^x] = C(B)*[t^B],
--   where A and B are matrices of floating-point numbers,
--   and C(A) and C(B) are complex matrices.

  procedure Random_Vector
              ( dim : in integer32;
                A,B : out Standard_Floating_Matrices.Matrix;
                x : out Standard_Floating_Vectors.Vector;
                cA,cB : out Standard_Complex_Matrices.Matrix;
                cx : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Generates a random matrix A of positive numbers in [0, 1],
  --   and leading powers in x, also in [0, 1], in one vector.
  --   The right hand side vectors in B are then computed.
  --   Random complex coefficients are generated in cA and cx,
  --   where cB is the result of a computation using cA and cx.

  procedure Random_Series
              ( dim,nbr : in integer32;
                X : out Standard_Floating_Matrices.Matrix;
                cX : out Standard_Complex_Matrices.Matrix );

  -- DESCRIPION :
  --   Generates dim many random series of nbr terms.

  -- ON ENTRY :
  --   dim      dimension of the vector of series;
  --   nbr      number of terms in each series;
  --   X        space for dim-by-nbr real matrix;
  --   cX       space for dim-by-nbr complex matrix.

  -- ON RETURN :
  --   X        matrix of exponents, range 1..dim and 1..nbr,
  --            where the rows are increasing;
  --   cX       a random dim-by-nbr matrix of coefficients.

  procedure Random_Series_System
              ( dim,nbr : in integer32;
                A,X,B : out Standard_Floating_Matrices.Matrix;
                cA,cX,cB : out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Generates a linear system of series with real powers 
  --   of dimension dim for a series with nbr terms.

  -- ON ENTRY :
  --   dim      dimension of the vector of series;
  --   nbr      number of terms in each series;
  --   A        space for dim-by-dim real matrix;
  --   cA       space for dim-by-dim complex matrix;
  --   X        space for dim-by-nbr real matrix;
  --   cX       space for dim-by-nbr complex matrix;
  --   B        space for dim-by-(nbr*dim) real matrix;
  --   cB       space for dim-by-(nbr*dim) complex matrix.

  -- ON RETURN :
  --   A        real powers of the leading terms of series,
  --            as the coefficients of a linear system;
  --   cA       coefficients of the leading terms, corresponding to A;
  --   X        matrix of exponents, range 1..dim and 1..nbr,
  --            where the rows are increasing;
  --   cX       a random dim-by-nbr matrix of coefficients.
  --   B        sorted exponents of the right hand sides,
  --            computed so that cX*[t^X] is a solution;
  --   cB       corresponding coefficients of the right hand sides.

  procedure Next_Coefficients_Check
              ( cy : in out Standard_Complex_Vectors.Vector;
                cA,cB : in Standard_Complex_Matrices.Matrix;
                idx1 : in Standard_Integer_Vectors.Vector;
                prev,next : in Boolean_Vectors.Vector;
                cx : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Computes the next coefficients of the leading power series,
  --   with a comparison to the correct values.

  -- ON ENTRY :
  --   cy       current values of the leading coefficients;
  --   cA       coefficients of the matrix of the linear system;
  --   cB       coefficients of the right hand side of the system;
  --   idx1     first set of indices of the tropical Cramer vector;
  --   prev     previously correct indices;
  --   next     current set of correct indices;
  --   cx       correct values of the leading coefficients,
  --            used for comparison.

  -- ON RETURN :
  --   cy       updated vector of coefficients.

  procedure Next_Series_Coefficients
              ( cy : in out Standard_Complex_Vectors.Vector;
                cA,cB : in Standard_Complex_Matrices.Matrix;
                idx1 : in Standard_Integer_Vectors.Vector;
                next : in Boolean_Vectors.Vector;
                correct : in out Standard_Integer_Vectors.Vector;
                cx : in Standard_Complex_Vectors.Vector;
                cZ : in out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Computes the next coefficients of the power series.

  -- ON ENTRY :
  --   cy       current values of the leading coefficients;
  --   cA       coefficients of the matrix of the linear system;
  --   cB       coefficients of the right hand side of the system;
  --   idx1     first set of indices of the tropical Cramer vector;
  --   next     current set of correct indices;
  --   correct  correct(i) indicates the number of correct values
  --            in the series cX*[t^X];
  --   cx       correct values of the leading coefficients,
  --            used for comparison;
  --   cZ       accumulates matrix of computed coefficients.

  -- ON RETURN :
  --   correct  updated indices of correct values;
  --   cy       updated vector of coefficients;
  --   cZ       updated matrix of computed coefficients.

  procedure Coefficient_Check
              ( tol : in double_float;
                cX,cY : in Standard_Complex_Vectors.Vector );
  procedure Coefficient_Check
              ( tol : in double_float;
                cX,cY : in Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Compares the computed coefficients cY to the original cX,
  --   using tol to decide whether a number is small enough.

  procedure Power_Check
              ( tol : in double_float;
                eX,eY : in Standard_Floating_Vectors.Vector );
  procedure Power_Check
              ( tol : in double_float;
                eX,eY : in Standard_Floating_Matrices.Matrix );

  -- DESCRIPTION :
  --   Compares the computed powers eY to the original eX,
  --   using tol to decide whether a number is small enough.

  procedure Test_Leading_Solver
              ( dim : in integer32; tol : in double_float;
                rA,rB : in Standard_Floating_Matrices.Matrix;
                cA,cB : in Standard_Complex_Matrices.Matrix;
                rx : in Standard_Floating_Vectors.Vector;
                cx : in Standard_Complex_Vectors.Vector;
                ry : out Standard_Floating_Vectors.Vector;
                cy : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Tests the solver to compute one vector of leading terms
  --   step-by-step, comparing with the generated solution.

  -- ON ENTRY :
  --   dim      dimension of the linear system;
  --   tol      tolerance to decide if zero or not;
  --   rA       leading exponents of the coefficient matrix;
  --   rB       exponents of the right-hand side;
  --   cA       leading coefficients of the coefficient matrix;
  --   cB       coefficients of the right-hand side;
  --   rx       exponents of the solution, for comparison;
  --   cx       coefficients of the solution, for comparison.

  -- ON RETURN :
  --   ry       computed exponents of the solution vector.
  --   cy       computed coefficients of the solution vector.

  procedure Test_Random_Vector ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates random input data of dimension dim,
  --   for one leading vector of a series.

  procedure Test_Random_Series ( dim,nbr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random linear system of dimension dim,
  --   for a series with nbr terms.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimension and then launches the test.

end Test_Leading_Terms;
