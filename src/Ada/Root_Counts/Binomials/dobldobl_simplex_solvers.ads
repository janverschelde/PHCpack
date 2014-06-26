with text_io;                          use text_io;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Double_Double_Numbers;            use Double_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer64_Matrices;
with Multprec_Integer_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Solutions;       use DoblDobl_Complex_Solutions;

package DoblDobl_Simplex_Solvers is

-- DESCRIPTION :
--   This package provides procedures to solve simplex systems
--   with standard complex coefficients.
--   There are 16 different solvers depending on 4 binary options:
--     1) no output or output to file: indicated by parameter file;
--     2) multiprecision Hermite normal form on exponent matrix or not;
--     3) lufac or lufco: respectively indicated by info or rcond;
--     4) extra output from binomial solver, many extra parameters.

  procedure Solve ( A : in Standard_Integer64_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; info : out integer32;
                    sols : out Solution_List );
  procedure Solve ( A : in Multprec_Integer_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; info : out integer32;
                    sols : out Solution_List );

  procedure Solve ( file : in file_type;
                    A : in Standard_Integer64_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; info : out integer32;
                    sols : out Solution_List );
  procedure Solve ( file : in file_type;
                    A : in Multprec_Integer_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; info : out integer32;
                    sols : out Solution_List );

  procedure Solve ( A : in Standard_Integer64_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; rcond : out double_double;
                    sols : out Solution_List );
  procedure Solve ( A : in Multprec_Integer_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; rcond : out double_double;
                    sols : out Solution_List );

  procedure Solve ( file : in file_type;
                    A : in Standard_Integer64_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; rcond : out double_double;
                    sols : out Solution_List );
  procedure Solve ( file : in file_type;
                    A : in Multprec_Integer_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; rcond : out double_double;
                    sols : out Solution_List );

  procedure Solve ( A : in Standard_Integer64_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; 
                    LUf : out DoblDobl_Complex_Matrices.Matrix;
                    pivots : out Standard_Integer_Vectors.Vector;
                    info : out integer32;
                    M,U : out Standard_Integer64_Matrices.Matrix;
                    Usols,Asols : out Solution_List );
  procedure Solve ( A : in Multprec_Integer_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; 
                    LUf : out DoblDobl_Complex_Matrices.Matrix;
                    pivots : out Standard_Integer_Vectors.Vector;
                    info : out integer32;
                    M,U : out Multprec_Integer_Matrices.Matrix;
                    Usols,Asols : out Solution_List );

  procedure Solve ( file : in file_type;
                    A : in Standard_Integer64_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; 
                    LUf : out DoblDobl_Complex_Matrices.Matrix;
                    pivots : out Standard_Integer_Vectors.Vector;
                    info : out integer32;
                    M,U : out Standard_Integer64_Matrices.Matrix;
                    Usols,Asols : out Solution_List );
  procedure Solve ( file : in file_type;
                    A : in Multprec_Integer_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; 
                    LUf : out DoblDobl_Complex_Matrices.Matrix;
                    pivots : out Standard_Integer_Vectors.Vector;
                    info : out integer32;
                    M,U : out Multprec_Integer_Matrices.Matrix;
                    Usols,Asols : out Solution_List );

  procedure Solve ( A : in Standard_Integer64_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; 
                    LUf : out DoblDobl_Complex_Matrices.Matrix;
                    pivots : out Standard_Integer_Vectors.Vector;
                    rcond : out double_double;
                    M,U : out Standard_Integer64_Matrices.Matrix;
                    Usols,Asols : out Solution_List );
  procedure Solve ( A : in Multprec_Integer_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; 
                    LUf : out DoblDobl_Complex_Matrices.Matrix;
                    pivots : out Standard_Integer_Vectors.Vector;
                    rcond : out double_double;
                    M,U : out Multprec_Integer_Matrices.Matrix;
                    Usols,Asols : out Solution_List );

  procedure Solve ( file : in file_type;
                    A : in Standard_Integer64_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; 
                    LUf : out DoblDobl_Complex_Matrices.Matrix;
                    pivots : out Standard_Integer_Vectors.Vector;
                    rcond : out double_double;
                    M,U : out Standard_Integer64_Matrices.Matrix;
                    Usols,Asols : out Solution_List );
  procedure Solve ( file : in file_type;
                    A : in Multprec_Integer_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    tol : in double_double; zero_y : out boolean;
                    r : out integer32; 
                    LUf : out DoblDobl_Complex_Matrices.Matrix;
                    pivots : out Standard_Integer_Vectors.Vector;
                    rcond : out double_double;
                    M,U : out Multprec_Integer_Matrices.Matrix;
                    Usols,Asols : out Solution_List );

  -- DESCRIPTION :
  --   Solves the simplex system C*x^A = b.

  -- ON ENTRY :
  --   file         for intermediate output and diagnostics;
  --   A            exponent vectors of the system;
  --   C            coefficient vector of the simplex system;
  --   b            constant righthand side vector;
  --   tol          tolerance to dismiss a complex floating point
  --                number as zero.

  -- ON RETURN :
  --   zero_y       true if the solution of C*y = b has zero components
  --                (smaller in modulus than tol) in which case there
  --                are no solutions with all components unequal zero,
  --                false otherwise;
  --   r            rank of the exponent matrix determines whether
  --                isolated solutions or not;
  --   LUf          LU factorization of the coefficient matrix;
  --   pivots       pivoting information of the LU factorization;
  --   info         returning info value from the LU factorization,
  --                if info = 0, then regular solutions,
  --                otherwise the coefficient matrix is singular;
  --   rcond        estimate for the inverse condition number of C;
  --   M            unimodular matrix to factor A;
  --   U            upper triangular matrix, M*U = A;
  --   Usols        solutions to x^U = y, with C*y = b;
  --   Asols        solutions to x^A = y, with C*y = b.

-- RESIDUAL CALCULATION :

  function Sum_Residuals
                  ( A : Standard_Integer64_Matrices.Matrix;
                    C : DoblDobl_Complex_Matrices.Matrix;
                    b : DoblDobl_Complex_Vectors.Vector;
                    sols : Solution_List ) return double_double;

  -- DESCRIPTION :
  --   Returns the sum of all residuals C*x^A - b,
  --   for x running through all solution vectors in sols.

  procedure Write_Residuals
                  ( file : in file_type;
                    A : in Standard_Integer64_Matrices.Matrix;
                    C : in DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    sols : in Solution_List );

  -- DESCRIPTION :
  --   Writes the residuals, max norm of C*x^A - b to the file,
  --   for all solutions in sols.

end DoblDobl_Simplex_Solvers;
