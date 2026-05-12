with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
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
                cB : out Standard_Complex_Matrices.Matrix;
                tosort : in boolean := false );

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
  --   cx       coefficients of t^x;
  --   tosort   if true, then the right hand side B is sorted
  --            according in increasing order of the powers.

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
                cB : out Standard_Complex_Matrices.Matrix;
                tosort : in boolean := false );

  -- DESCRIPTION :
  --   Makes the product of the series cA*[t^A] with cX*[t^X]
  --   and stores the result in cB*[t^B].
 
  -- REQUIRED :
  --   If A is a dim-by-dim matrix and X is a dim-by-nbr matrix,
  --   then B is a dim-by-(dim*nbr) matrix.
  --   The coefficient matrices cA, cX, and cB have the same ranges
  --   as the exponent matrices A, X, and B.
  --   If tosort, then the right hand side (B, cB) is sorted
  --   in increasing order of the exponents.

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
                cA : in Standard_Complex_Matrices.Matrix;
                cb : in Standard_Complex_Vectors.Vector;
                idx1 : in Standard_Integer_Vectors.Vector;
                prev,next : in Boolean_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the next coefficients of the leading power series.

  -- ON ENTRY :
  --   cy       current values of the leading coefficients;
  --   cA       coefficients of the matrix of the linear system;
  --   cb       selected coefficients of the right hand side of the system,
  --            selected corresponding to the next minimal powers;
  --   idx1     first set of indices of the tropical Cramer vector;
  --   prev     previously correct indices;
  --   next     current set of correct indices;
  --   vrblvl   is the verbose level, silent if zero.

  -- ON RETURN :
  --   cy       updated vector of coefficients.

  procedure Next_Series_Coefficients
              ( cy : in out Standard_Complex_Vectors.Vector;
                cA,cB : in Standard_Complex_Matrices.Matrix;
                idx1,Bidx : in Standard_Integer_Vectors.Vector;
                next : in Boolean_Vectors.Vector;
                correct : in out Standard_Integer_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the next coefficients of the power series.

  -- ON ENTRY :
  --   cy       current values of the leading coefficients;
  --   cA       coefficients of the matrix of the linear system;
  --   cB       coefficients of the right hand side of the system;
  --   idx1     first set of indices of the tropical Cramer vector;
  --   Bidx     Bidx(i) is the column index in the right hand size
  --            for the i-th row;
  --   next     current set of correct indices;
  --   correct  correct(i) indicates the number of correct values
  --            in the series cX*[t^X];
  --   vrblvl   is the verbose level, silent if zero.

  -- ON RETURN :
  --   correct  updated indices of correct values;
  --   cy       updated vector of coefficients.

  function Leading_Right_Power
             ( rB : Standard_Floating_Matrices.Matrix; rowidx : integer32;
               skipcols : Boolean_Vectors.Vector;
               vrblvl : integer32 := 0 ) return double_float;

  -- DESCRIPTION :
  --   Returns the smallest number on the row of rB, with index rowidx,
  --   skipping the columns k for which skipcols(k) is true.

  procedure Leading_Right_Term
              ( rB : in Standard_Floating_Matrices.Matrix;
                cB : in Standard_Complex_Matrices.Matrix;
                rowidx : in integer32;
                skipcols : in Boolean_Vectors.Vector;
                minpow : out double_float; mincff : out Complex_Number;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Returns in minpow the smallest number on the row of rB,
  --   with index rowidx, skipping the columns k for which skipcols(k)
  --   is true, and returns in mincff the corresponding coefficients
  --   on the row of cB with index rowidx.

  function Leading_Right_Powers
             ( rB : Standard_Floating_Matrices.Matrix;
               skipcols : Boolean_Vectors.Vector; vrblvl : integer32 := 0 )
             return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector of smallest numbers on the rows of rB,
  --   skipping the columns in skipcols.

  procedure Leading_Right_Terms
              ( rB : in Standard_Floating_Matrices.Matrix;
                cB : in Standard_Complex_Matrices.Matrix;
                skipcols : in Boolean_Vectors.Vector;
                minpow : out Standard_Floating_Vectors.Vector;
                mincff : out Standard_Complex_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Returns in minpow the smallest number on the row of rB,
  --   skipping the columns k for which skipcols(k)is true, and
  --   returns in mincff the corresponding coefficients of cB.

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

  -- REQUIRED :
  --   The exponents of the right hand side vector must not be sorted.
  --   They must correspond to the order in which the coefficients of A
  --   are multiplied with the solution vector.

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

  procedure Right_Index_Terms
              ( rBidx : in out Standard_Integer_Vectors.Vector;
                rA : in Standard_Floating_Matrices.Matrix;
                rB : in out Standard_Floating_Matrices.Matrix;
                cB : in out Standard_Complex_Matrices.Matrix;
                vY : in Standard_Floating_Vectors.Vector;
                next : in Boolean_Vectors.Vector; tol : in double_float;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Updates the indices in the exponents of the right hand side,
  --   based on new computed exponents of the solution.
  --   Swaps the elements in rB and cB if the minima are not at
  --   the proper positions.

  function Locate_Right_Power
              ( B : Standard_Floating_Matrices.Matrix; x : double_float;
                row,startcol : integer32; tol : double_float )
              return integer32;
 
  -- DESCRIPTION :
  --   Returns the index of the element x on the given row of B,
  --   starting the search at the column startcol.

  procedure Series_Solver
              ( dim,nbr : in integer32; tol : in double_float;
                rA,rB : in Standard_Floating_Matrices.Matrix;
                cA,cB : in Standard_Complex_Matrices.Matrix;
                rY : out Standard_Floating_Matrices.Matrix;
                cY : out Standard_Complex_Matrices.Matrix;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Compute the first nbr terms of a dim-dimensional series vector,
  --   solving a linear system with real-powered series coefficients.

  -- REQUIRED : rB'range = 1..dim*nbr = cB'range,
  --   (rB, cB) stores the product of the coefficient matrix (rA, cA)
  --   with the solution series (rX, cX), and moreover: the rows
  --   of rB (with their corresponding cB) are in increasing order.

  -- ON ENTRY :
  --   dim      dimension of the linear system;
  --   nbr      number of terms in each component of the solution;
  --   tol      tolerance to decide if zero or not;
  --   rA       leading exponents of the coefficient matrix;
  --   rB       exponents of the right-hand side;
  --   cA       leading coefficients of the coefficient matrix;
  --   cB       coefficients of the right-hand side;
  --   vrblvl   is the verbose level, silent if zero.

  -- ON RETURN :
  --   rY       computed exponents of the solution vector.
  --   cY       computed coefficients of the solution vector.

end Double_Puiseux_Operations;
