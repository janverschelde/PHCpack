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

  procedure Sort ( B : in out Standard_Floating_Matrices.Matrix;
                   cB : in out Standard_Complex_Matrices.Matrix;
                   nbrcols : in integer32 );

  -- DESCRIPTION :
  --   Sorts the rows in B in increasing order,
  --   also swapping the corresponding coefficients in cB.
  --   The number of valid columns in B is nbrcols.

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

  procedure Series_Product
              ( A,X : in Standard_Floating_Matrices.Matrix;
                cA,cX : in Standard_Complex_Matrices.Matrix;
                B : out Standard_Floating_Matrices.Matrix;
                cB : out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Makes the product of the series cA*[t^A] with cX*[t^X]
  --   and stores the result in cB*[t^B].
 
  -- REQUIRED :
  --   If A is a dim-by-dim matrix and X is a dim-by-nbr matrix,
  --   then B is a dim-by-(dim*nbr) matrix.
  --   The coefficient matrices cA, cX, and cB have the same ranges
  --   as the exponent matrices A, X, and B.

  procedure Series_Product
              ( A,X : in Standard_Floating_Matrices.Matrix;
                cA,cX : in Standard_Complex_Matrices.Matrix;
                skip : in Standard_Integer_Vectors.Vector;
                B : out Standard_Floating_Matrices.Matrix;
                cB : out Standard_Complex_Matrices.Matrix;
                nbrcols : out integer32 );

  -- DESCRIPTION :
  --   Similar to Series_Product(A,X,cA,cX,B,cB),
  --   skipping the columns according to the indices in skip.
  --   The number of valid columns in B and cB is in nbrcols.

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

  procedure Next_Coefficients
              ( cy : in out Standard_Complex_Vectors.Vector;
                cA,cB : in Standard_Complex_Matrices.Matrix;
                idx1 : in Standard_Integer_Vectors.Vector;
                prev,next : in Boolean_Vectors.Vector;
                cx : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Computes the next coefficients of the leading power series.

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
                prev,next : in Boolean_Vectors.Vector;
                correct : in out Standard_Integer_Vectors.Vector;
                cx : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Computes the next coefficients of the power series.

  -- ON ENTRY :
  --   cy       current values of the leading coefficients;
  --   cA       coefficients of the matrix of the linear system;
  --   cB       coefficients of the right hand side of the system;
  --   idx1     first set of indices of the tropical Cramer vector;
  --   prev     previously correct indices;
  --   next     current set of correct indices;
  --   correct  correct(i) indicates the number of correct values
  --            in the series cX*[t^X];
  --   cx       correct values of the leading coefficients,
  --            used for comparison.

  -- ON RETURN :
  --   correct  updated indices of correct values;
  --   cy       updated vector of coefficients.

  procedure Coefficient_Check
              ( tol : in double_float;
                cx,cy : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Compares the computed coefficients cy to the original cx,
  --   using tol to decide whether a number is small enough.

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
