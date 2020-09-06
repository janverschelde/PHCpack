with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Integer_Matrices;

package Test_Integer_Linear_Solvers is

-- DESCRIPTION :
--   Tests on integer linear algebra operations.

  procedure Random_Test_Standard_Solver
                 ( n,m,low,upp : in integer32 );

  -- DESCRIPTION :
  --   Generates a random n-by-m 32-bit integer matrix,
  --   with numbers between low and upp.
  --   Computes a vector in the kernel of the matrix
  --   using integer row reduction operations.

  procedure Random_Test_Standard64_Solver
                 ( n,m : in integer32; low,upp : in integer64 );

  -- DESCRIPTION :
  --   Generates a random n-by-m 64-bit integer matrix,
  --   with numbers between low and upp.
  --   Computes a vector in the kernel of the matrix
  --   using integer row reduction operations.

  procedure Random_Test_Multprec_Solver
              ( n,m : in integer32; sz : in natural32 );

  -- DESCRIPTION :
  --   Generates a random n-by-m multiprecision integer matrix,
  --   with numbers of size sz.
  --   Computes a vector in the kernel of the matrix
  --   using integer row reduction operations.

  procedure Random_Test_Standard_Solvers;

  -- DESCRIPTION :
  --   Prompts for the number of tests, the dimensions of the matrices
  --   and the lower and upper bounds on the 32-bit entries.
  --   Then runs a number of tests on 32-bit integer matrices.

  procedure Random_Test_Standard64_Solvers;

  -- DESCRIPTION :
  --   Prompts for the number of tests, the dimensions of the matrices
  --   and the lower and upper bounds on the 64-bit entries.
  --   Then runs a number of tests on 64-bit integer matrices.

  procedure Random_Test_Multprec_Solvers;

  -- DESCRIPTION :
  --   Prompts for the number of tests, the dimensions of the matrices
  --   and the size of the multiprecision random numbers.
  --   Then runs a number of tests on multiprecision integer matrices.

  procedure Interactive_Test_Standard_Solvers;

  -- DESCRIPTION :
  --   Prompts for the number of rows, the number of columns,
  --   and all entries of a 32-bit integer matrix
  --   to test the row reduction algorithms.

  procedure Interactive_Test_Standard64_Solvers;

  -- DESCRIPTION :
  --   Prompts for the number of rows, the number of columns,
  --   and all entries of a 64-bit integer matrix
  --   to test the row reduction algorithms.

  procedure Show_Differences
              ( A,B : in Multprec_Integer_Matrices.Matrix );

  -- DESCRIPTION :
  --   Shows the differences between the two matrices A and B.

  procedure Interactive_Test_Multprec_Solvers;

  -- DESCRIPTION :
  --   Prompts for the number of rows, the number of columns,
  --   and all entries of a 64-bit integer matrix
  --   to test the row reduction algorithms.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Integer_Linear_Solvers;
