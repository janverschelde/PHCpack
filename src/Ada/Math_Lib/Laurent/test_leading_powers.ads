with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Matrices;        use Standard_Floating_Matrices;

package Test_Leading_Powers is

-- DESCRIPTION :
--   Tests the computation of the leading powers of a linear system
--   of real power series.

  function Row_Min_Plus
             ( A : Matrix; x : Standard_Floating_Vectors.Vector )
             return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the minimum of the sum of the elements on each row of A,
  --   augmented in each column with the corresponding element in x,
  --   representing the leading powers of t^A*t^x.

  function id ( n : integer32 ) return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the identity permutation, with numbers 1 to n.

  procedure Shuffle ( p : in out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Shuffles the indices in p, via p'last random swaps.

  function Random_Permutation
             ( n : integer32 ) return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a random permutation of n elements.

  procedure Random_Leading_Input 
              ( dim : in integer32;
                p : out Standard_Integer_Vectors.Vector;
                A : out Matrix;
                x,b : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Generates random input so the leading powers of a series
  --   with real powers can be recovered by one single Cramer vector.

  -- ON ENTRY :
  --   dim      dimension, number of rows and columns of A;
  --   p        permutation, p(i) is the column where the mininum
  --            will occur on row i of A;
  --   A        matrix with prescribed locations of the minima;
  --   x        sufficiently small powers for the minimum to
  --            occur exactly twice in the [A|b]+c matrix;
  --   b        leading powers of the right hand side vector.

  procedure Random_General_Input 
              ( dim : in integer32; A : out Matrix;
                x,b : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Generates a random matrix A of positive numbers in [0, 1],
  --   and leading powers in x, also in [0, 1].
  --   The right hand side vector is then computed.
  --   One single Cramer vector is unlikely to suffice to recover
  --   all leading powers for this general input.

  procedure Leading_Powers
              ( dim : in integer32; A : in Matrix;
                b : in Standard_Floating_Vectors.Vector;
                d : out Standard_Floating_Vectors.Vector;
                idxone,idxtwo : out Standard_Integer_Vectors.Vector;
                fail : out boolean );

  -- DESCRIPTION :
  --   Returns in d the leading powers via a tropical Cramer vector
  --   on A and b.

  -- ON ENTRY :
  --   dim      dimension, number of rows and columns of A;
  --   A        matrix with prescribed locations of the minima;
  --   b        leading powers of the right hand side vector.

  -- ON RETURN :
  --   d        leading powers of the solution,
  --            computed via the tropical Cramer vector.
  --   idxone   indices where the minimum is first obtained;
  --   idxtwo   indices where the minimum is obtained the second time;
  --   fail     true if the minimum is not everywhere exactly obtained twice,
  --            false if the minimum is obtained exactly twice, everywhere.

  procedure Check_Differences
              ( dim : in integer32; A : in Matrix;
                b,x,d : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Checks the differences between the original x
  --   and the computed leading degrees in d.

  -- ON ENTRY :
  --   dim      dimension, number of rows and columns of A;
  --   A        matrix with prescribed locations of the minima;
  --   b        leading powers of the right hand side vector;
  --   x        original leading powers of the solution;
  --   d        leading powers of the solution,
  --            computed via the tropical Cramer vector.

  procedure Test_Leading_Random ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates a random problem of dimension dim and then computes 
  --   the tropical Cramer vector to recover the leading powers.
  --   In this test case, one Cramer vector should suffice.

  procedure Test_General_Random ( dim : in integer32 );

  -- DESCRIPTION :
  --   Generates a random problem of dimension dim and then attempts
  --   to recover the leading powers.
  --   In this test case, one Cramer vector will not suffice.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a dimension and then launches the test.

end Test_Leading_Powers;
