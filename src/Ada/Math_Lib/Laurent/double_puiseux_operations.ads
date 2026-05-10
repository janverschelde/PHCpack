with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Boolean_Vectors;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;

package Double_Puiseux_Operations is

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

  procedure Leading_Powers
              ( dim : in integer32; tol : in double_float;
                A : in Standard_Floating_Matrices.Matrix;
                b : in Standard_Floating_Vectors.Vector;
                wrkm : in Standard_Integer_VecVecs.VecVec;
                d : out Standard_Floating_Vectors.Vector;
                idxone,idxtwo : out Standard_Integer_Vectors.Vector;
                fail : out boolean; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Returns in d the leading powers via a tropical Cramer vector
  --   on A and b.

  -- ON ENTRY :
  --   dim      dimension, number of rows and columns of A;
  --   tol      tolerance to decide whether a number is zero, or equivalently,
  --            whether two numbers are different from each other;
  --   A        matrix of real powers;
  --   b        leading powers of the right hand side vector;
  --   wrkm     pointers to work space for the indices in the
  --            computation of the tropical Cramer vector,
  --            the range of wrkm is 0..dim;
  --   vrblvl   is the verbose level, silent when zero.

  -- ON RETURN :
  --   d        leading powers of the solution,
  --            computed via the tropical Cramer vector.
  --   idxone   indices where the minimum is first obtained;
  --   idxtwo   indices where the minimum is obtained the second time;
  --   fail     true if the minimum is not everywhere exactly obtained twice,
  --            false if the minimum is obtained exactly twice, everywhere.

  procedure Check_Correctness
              ( dim : in integer32; tol : in double_float;
                x,d : in Standard_Floating_Vectors.Vector;
                idx1,idx2 : in Standard_Integer_Vectors.Vector;
                correct : out Boolean_Vectors.Vector;
                cd : in out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Determines which of the entries in d are correct,
  --   using the indices computed by the tropical Cramer vector,
  --   updating the Boolean vector correct.

  -- ON ENTRY :
  --   dim      dimension, number of rows and columns of A;
  --   tol      tolerance to decide whether a number is zero, or equivalently,
  --            whether two numbers are different from each other;
  --   A        matrix with prescribed locations of the minima;
  --   x        original leading powers of the solution;
  --   d        leading powers of the solution,
  --            computed via the tropical Cramer vector;
  --   idx1     first set of indices in the tropical Cramer vector;
  --   idx2     second set of indices in the tropical Cramer vector;
  --   cd       current correct values;
  --   vrblvl   is the verbose level, silent if zero.

  -- ON RETURN :
  --   correct  updated vector of indices to correct values;
  --   cd       updated correct values.

  procedure Assign_Correctness
              ( dim : in integer32;
                d : in Standard_Floating_Vectors.Vector;
                idx1,idx2 : in Standard_Integer_Vectors.Vector;
                correct : out Boolean_Vectors.Vector;
                cd : in out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Determines which of the entries in d are correct,
  --   using the indices computed by the tropical Cramer vector,
  --   updating the Boolean vector correct.

  -- ON ENTRY :
  --   dim      dimension, number of rows and columns of A;
  --   A        matrix with prescribed locations of the minima;
  --   d        leading powers of the solution,
  --            computed via the tropical Cramer vector;
  --   idx1     first set of indices in the tropical Cramer vector;
  --   idx2     second set of indices in the tropical Cramer vector;
  --   cd       current correct values;
  --   vrblvl   is the verbose level, silent if zero.

  -- ON RETURN :
  --   correct  updated vector of indices to correct values;
  --   cd       updated correct values.

  procedure Next_Coefficients
              ( cy : in out Standard_Complex_Vectors.Vector;
                cA,cB : in Standard_Complex_Matrices.Matrix;
                idx1 : in Standard_Integer_Vectors.Vector;
                prev,next : in Boolean_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the next coefficients of the leading power series.

  -- ON ENTRY :
  --   cy       current values of the leading coefficients;
  --   cA       coefficients of the matrix of the linear system;
  --   cB       coefficients of the right hand side of the system;
  --   idx1     first set of indices of the tropical Cramer vector;
  --   prev     previously correct indices;
  --   next     current set of correct indices;
  --   vrblvl   is the verbose level, silent if zero.

  -- ON RETURN :
  --   cy       updated vector of coefficients.

  procedure Leading_Solver
              ( dim : in integer32; tol : in double_float;
                rA,rB : in Standard_Floating_Matrices.Matrix;
                cA,cB : in Standard_Complex_Matrices.Matrix;
                ry : out Standard_Floating_Vectors.Vector;
                cy : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Returns the leading powers and coefficients of the solution
  --   to a linear system with real-powered series coefficients.

  -- ON ENTRY :
  --   dim      dimension, number of rows and columns of A;
  --   tol      tolerance to decide if a number is zero or not;
  --   rA       leading powers of the coefficient matrix A;
  --   rB       powers of the right hand side of the linear system;
  --   cA       leading coefficients of the coefficient matrix A;
  --   cB       coefficients of the right hand side of the linear system;
  --   vrblvl   is the verbose level, silent if zero.

  -- ON RETURN :
  --   ry       leading powers of the solution vector;
  --   cy       leading coefficients of the solution vector.

end Double_Puiseux_Operations;
