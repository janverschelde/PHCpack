with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;

package Test_Standard_Singular_Values is

-- DESCRIPTION :
--   Tests the Singular Value Decomposition in standard double precision.

  function Read_Vector
             ( n : in integer32 ) return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Prompts the user to enter n complex numbers,
  --   returns the vector which contains the user given numbers.

  function Read_Matrix
             ( n,m : in integer32 ) return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Prompts the user to enter n*m complex numbers to define
  --   an n-by-m matrix.  Returns the matrix of given numbers.

  function Is_Identity ( a : Standard_Complex_Matrices.Matrix;
                         tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if a is the identity matrix modulo the tolerance.

  function Is_Orthogonal ( a : Standard_Complex_Matrices.Matrix;
                           tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if a*a^H and a^H*a equal the identity matrix
  --   modulo the given tolerance tol.

  function Is_SVD ( x,u,v : Standard_Complex_Matrices.Matrix;
                    s : Standard_Complex_Vectors.Vector;
                    tol : in double_float ) return boolean;

  -- DESCRIPTION :
  --   Forms the product u'*x*v and checks whether this product is a
  --   diagonal matrix with the singular values on the diagonal.

  procedure Test_SVD_Output
              ( x,u,v : in Standard_Complex_Matrices.Matrix;
                s,e : in Standard_Complex_Vectors.Vector;
                info : in integer32; output : in boolean := true );

  -- DESCRIPTION :
  --   Tests whether the matrices u and v are orthogonal and
  --   whether u'*x*v yields a diagonal matrix with the singular
  --   values on the diagonal.  Diagnostics are written to screen.
  --   If output, then all vectors and matrices are written.

  procedure Test_SVD_Solver
               ( a,u,v : in Standard_Complex_Matrices.Matrix;
                 s,b : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Calls the solver using the SVD of a to solve a*x = b,
  --   and prints the residual vector b - a*x.

  procedure Test_SVD_on_Given_Matrix ( n,p : in integer32 );

  -- DESCRIPTION :
  --   Computes the SVD of a user given n-by-p matrix.

  procedure Test_SVD_on_Given_System ( n,p : in integer32 );

  -- DESCRIPTION :
  --   Uses SVD to solve a given linear system of n equations in p unknowns.

  procedure Test_SVD_on_Random_Matrix ( n,p : in integer32 );

  -- DESCRIPTION :
  --   Computes a SVD on a random n-by-p matrix.

  procedure Test_SVD_on_Random_System ( n,p : in integer32 );

  -- DESCRIPTION :
  --   Uses SVD to solve a randomly generated linear system
  --   of n equations in p unknowns.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_Standard_Singular_Values;
