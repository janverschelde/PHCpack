with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Boolean_Vectors;
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

  procedure Series_Product
              ( A : in Standard_Floating_Matrices.Matrix;
                x : in Standard_Floating_Vectors.Vector;
                cA : in Standard_Complex_Matrices.Matrix;
                cx : in Standard_Complex_Vectors.Vector;
                B : out Standard_Floating_Matrices.Matrix;
                cB : out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Makes the product of the series in the power matrix A,
  --   coefficients in cA, with the series in the power vector x, 
  --   coefficients in cx, to return the power matrix B and cB
  --   as the corresponding coefficients.
  --   The columns in the power matrix B are sorted in increasing order.

  -- REQUIRED :
  --    All ranges of matrices and vectors are identical.
  --    Powers and coefficients are sufficiently generic.

  -- ON ENTRY:
  --   A        a power matrix in a linear system t^A;
  --   x        powers of the column vector t^x;
  --   cA       coefficients of t^A;
  --   cx       coefficients of t^x.

  -- ON RETURN :
  --   B        power matrix of [t^A]*[t^x], sorted columnwise;
  --   cB       corresponding coefficients of B.

  procedure Series_Product
              ( A : in Standard_Floating_Matrices.Matrix;
                x : in Standard_Floating_Vectors.Vector;
                cA : in Standard_Complex_Matrices.Matrix;
                cx : in Standard_Complex_Vectors.Vector;
                skip : in Boolean_Vectors.Vector;
                B : out Standard_Floating_Matrices.Matrix;
                cB : out Standard_Complex_Matrices.Matrix;
                nbrcols : out integer32 );

  -- DESCRIPTION :
  --   Similar as Series_Product(A,x,cA,cx,B,cB), columns k for which
  --   skip(k) is true are skipped.  The number of columns nbrcols in
  --   B and cB that count may be less than B'last(2).

  procedure Random_Input
              ( dim : in integer32;
                A,B : out Standard_Floating_Matrices.Matrix;
                x : out Standard_Floating_Vectors.Vector;
                cA,cB : out Standard_Complex_Matrices.Matrix;
                cx : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Generates a random matrix A of positive numbers in [0, 1],
  --   and leading powers in x, also in [0, 1].
  --   The right hand side vectors in B are then computed.
  --   Random complex coefficients are generated in cA and cx,
  --   where cB is the result of a computation using cA and cx.

  procedure Test_Random_Input ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates random input data of dimension dim.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimension and then launches the test.

end Test_Leading_Terms;
