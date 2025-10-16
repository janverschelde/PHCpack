with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;        use Standard_Floating_Matrices;

package Double_Weighted_Assignment is

-- DESCRIPTION :
--   Applies the Hungarian algorithm to solve the weighted assignment
--   problem in double arithmetic.

  procedure enumerate ( A : in Matrix; k : in integer32;
                        selcols : out Standard_Integer_Vectors.Vector;
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
  --   selsum   should be initialized to zero;
  --   minsum   should be initialized to a very large number,
  --            larger than any element in the matrix A;
  --   vrblvl   the verbose level, if positive, writes the columns
  --            and the corresponding minimum, each time a new minimum
  --            is reached.

  -- ON RETURN :
  --   selcols  work space of selected columns;
  --   mincols  selected columns with the minimum sum of weights;
  --   minsum   the minimum sum of weights.

  function Value_Selection
             ( A : Matrix; cols : Standard_Integer_Vectors.Vector ) 
             return double_float;

  -- DESCRIPTION :
  --   Given the weights in the matrix A and columns in cols,
  --   returns the value of the selected columns,
  --   which should match the minimum cost.

  procedure Hungarian ( A : in Matrix;
                        u,v : out Standard_Floating_Vectors.Vector;
                        p,way : out Standard_Integer_Vectors.Vector;
                        vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Runs the Hungarian algorithm to solve the weighted assignment
  --   problem with weights given in the matrix A.

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

end Double_Weighted_Assignment;
