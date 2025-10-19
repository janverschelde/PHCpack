with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;        use Standard_Floating_Matrices;

package Double_Weighted_Assignment is

-- DESCRIPTION :
--   Applies the Hungarian algorithm to solve the weighted assignment
--   problem in double arithmetic.

  procedure enumerate
              ( A : in Matrix; k : in integer32;
                selcols : in out Standard_Integer_Vectors.Vector;
                selsum : in double_float;
                mincols : out Standard_Integer_Vectors.Vector;
                minsum : in out double_float;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs a brute force method enumerating all choices of columns,
  --   running over all rows of A, with the current selection of columns
  --   in selcols and the current sum of the selection is selsum.

  -- ON ENTRY :
  --   A        matrix with weights;
  --   k        current column to consider, start at 1;
  --   selcols  work space for the selected columns;
  --   selsum   should be initialized to zero;
  --   minsum   should be initialized to a very large number,
  --            larger than any element in the matrix A;
  --   vrblvl   the verbose level, if positive, writes the columns
  --            and the corresponding minimum, each time a new minimum
  --            is reached.

  -- ON RETURN :
  --   selcols  updated work space of selected columns;
  --   mincols  selected columns with the minimum sum of weights;
  --   minsum   the minimum sum of weights.

  function Value_Selection
              ( A : Matrix; cols : Standard_Integer_Vectors.Vector ) 
              return double_float;

  -- DESCRIPTION :
  --   Given the weights in the matrix A and columns in cols,
  --   returns the value of the selected columns,
  --   which should match the minimum cost.

  procedure Hungarian
              ( A : in Matrix;
                u,v : out Standard_Floating_Vectors.Vector;
                p,way : out Standard_Integer_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs the Hungarian algorithm to solve the weighted assignment
  --   problem with weights given in the matrix A.

  -- ACKNOWLEDGEMENT :
  --   The code is based on the description at
  --   https://cp-algorithms.com/graph/hungarian-algorithm.html
  --   which attributes the short code for the Hungarian algorithm
  --   to Andrey Lopatin.

  -- ON ENTRY :
  --   A        matrix with weights;
  --   u        vector of range 0..A'last(1);
  --   v        vector of range 0..A'last(2);
  --   p        vector of range 0..A'last(2);
  --   way      vector of range 1..A'last(2);
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   u,v      defines the potential: u(i) + v(j) <= A(i,j),
  --            for all indices i and j;
  --   p        stores the matching;
  --   way      defines the augmenting path.

  function Row_Matching
              ( cols : Standard_Integer_Vectors.Vector;
                nbrows : integer32 )
              return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Given in cols the selected columns, as in p
  --   computed by running the Hungarian algorithm,
  --   returns an integer vector of range 1..nbrows,
  --   defining the matching row wise.

  procedure Cramer_Vector
              ( A : in Matrix;
                b : in Standard_Floating_Vectors.Vector;
                c : out Standard_Floating_Vectors.Vector;
                m : in Standard_Integer_VecVecs.VecVec;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies the Hungarian algorithm to compute the Cramer vector
  --   of the system with matrix A and right hand side b.

  -- ON ENTRY :
  --   A        matrix with weights;
  --   b        right hand side vector;
  --   c        vector of range 0..A'last(2);
  --   m        space of A'last(2)+1 vectors, available as
  --            a vector of range 0..A'last(2).

  -- ON RETURN :
  --   c        c(k) equals the minimum costs,
  --            where the k-th column of A is replaced by b,
  --            c(0) is the minimum cost of A;
  --   m        m(0) is the matching of A, and
  --            m(k) is the matching of replacing column k by b.

  function Second_Index
              ( m : Standard_Integer_VecVecs.VecVec )
              return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   The first index where the minimum is attained in [A|b]+c
  --   is in the vector m(0).  The second index is determined by
  --   the different indices in the other vectors in m,
  --   and all second indices are in the returned vector.

  function Abc_Matrix ( A : Matrix;
                        b,c : Standard_Floating_Vectors.Vector )
                      return Matrix;

  -- DESCRIPTION :
  --   The Abc matrix adds c to the columns of A extended with b,
  --   where c is the tropical Cramer vector.
  --   The matrix on return has one extra column than A.

  procedure Abc_Argmin
              ( abc : in Matrix; tol : in double_float;
                minvals : out Standard_Floating_Vectors.Vector;
                idxmins : out Standard_Integer_Vectors.Vector;
                fail : out boolean; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given in abc is the sum of the weights in A extended with b,
  --   augmented with the values of the Cramer vector,
  --   computes the minimum on each row and the indices
  --   where the minimum is attained.

  -- ON ENTRY :
  --   abc      [A|b]+c or the Abc matrix;
  --   tol      tolerance to decide whether two numbers are the same;
  --   minvals  space for as many doubles as the number of rows in abc;
  --   idxmins  space for twice as many doubles as the number of rows;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   minvals  values of the minimum on each row of abc;
  --   idxmins  pairs of indices where the minimum is attained twice,
  --            the first index appears on odd positions,
  --            followed by the second index on even positions;
  --   fail     true if the minimum is attained only once,
  --            or more than twice.

  procedure Abc_Indices
              ( abc : in Matrix; tol : in double_float;
                m : in Standard_Integer_VecVecs.VecVec;
                minvals : in Standard_Floating_Vectors.Vector;
                idx1,idx2 : out Standard_Integer_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the indices where the minimum in abc is attained,
  --   based on the indices of the Cramer vector,
  --   with a correction and a check on the right hand side index.

  -- ON ENTRY :
  --   abc      [A|b]+c or the Abc matrix;
  --   tol      tolerance to decide whether two numbers are the same;
  --   m        index vectors as output on the c vector computation;
  --   minvals  space for as many doubles as the number of rows in abc;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   idx1     first index where the minimum is attained;
  --   idx2     second index where the minimum is attained.

end Double_Weighted_Assignment;
