with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_Matrices;

package Test_Multprec_Singular_Values is

-- DESCRIPTION :
--   Tests the Singular Value Decomposition in arbitrary multiprecision.

  procedure Set_Size ( a : in out Multprec_Complex_Matrices.Matrix;
                       size : in natural32 );

  -- DESCRIPTION :
  --   Converts the numbers in a to the precision given by size.

  function Read_Vector
             ( n : in integer32 ) return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Prompts the user to enter n complex numbers to define
  --   an n-vector.  Returns the vector of given numbers.

  function Read_Matrix
             ( n,m : in integer32 ) return Multprec_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Prompts the user to enter n*m complex numbers to define
  --   an n-by-m matrix.  Returns the matrix of given numbers.

  function Is_Identity ( a : Multprec_Complex_Matrices.Matrix;
                         tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the matrix a is the identity matrix,
  --   modulo the given tolerance tol.

  function Is_Orthogonal ( a : Multprec_Complex_Matrices.Matrix;
                           tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if a*a^H and a^H*a equal the identity matrix,
  --   modulo the given tolerance tol.

  function Is_SVD ( x,u,v : Multprec_Complex_Matrices.Matrix;
                    s : Multprec_Complex_Vectors.Vector;
                    tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if u^H*x*v is a diagonal matrix with the singular
  --   values in s on its diagonal, modulo the tolerance tol.

  procedure Test_SVD_Output
              ( x,u,v : in Multprec_Complex_Matrices.Matrix;
                s,e : in Multprec_Complex_Vectors.Vector;
                info : in integer32; output : in boolean := true );

  -- DESCRIPTION :
  --   Tests whether the matrices u and v are orthogonal and
  --   whether u^H*x*v yields a diagonal matrix with the singular
  --   values on the diagonal.  Diagnostics are written to screen.
  --   If output, then all vectors and matrices are written.

  procedure Test_SVD_Solver
              ( a,u,v : in Multprec_Complex_Matrices.Matrix;
                s,b : in Multprec_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Calls the solver to use SVD of a to solve the system a*x = b,
  --   and performs several diagnostic tests.

  procedure Test_SVD_on_Random_Matrix
              ( n,p : in integer32; size : in natural32 );

  -- DESCRIPTION :
  --   Uses SVD to solve a randomly generated linear system
  --   of n equations in p unknowns with number of the given size.

  procedure Test_SVD_on_Random_System
              ( n,p : in integer32; size : in natural32 );

  -- DESCRIPTION :
  --   Uses SVD to solve a randomly generated linear system
  --   of n equations in p unknowns with number of the given size.

  procedure Test_SVD_on_Given_System
             ( n,p : in integer32; size : in natural32 );

  -- DESCRIPTION :
  --   Uses SVD to solve a randomly generated linear system
  --   of n equations in p unknowns with number of the given size.

  procedure Test_SVD_on_Given_Matrix
              ( n,p : in integer32; size : in natural32 );

  -- DESCRIPTION :
  --   Uses SVD to solve a user given linear system
  --   of n equations in p unknowns with number of the given size.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Multprec_Singular_Values;
