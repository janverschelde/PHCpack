with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers; 
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers; 
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Integer_Matrices;
with Standard_Integer64_Matrices;
with Lists_of_Integer_Vectors;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Standard_Binomial_Varieties is

-- DESCRIPTION :
--   Tools are offered to manipulate algebraic sets defined by binomial 
--   systems, using standard 32-bit or 64-bit arithmetic on the exponents.
--   The blackbox solvers are at the end.

-- STAGE I : computing and checking a cone of tropisms

  procedure Cone_of_Tropisms 
               ( A : in Standard_Integer_Matrices.Matrix;
                 rank : out integer32;
                 V : out Standard_Integer_Matrices.Link_to_Matrix );
  procedure Cone_of_Tropisms 
               ( A : in Standard_Integer64_Matrices.Matrix;
                 rank : out integer32;
                 V : out Standard_Integer64_Matrices.Link_to_Matrix );

  -- DESCRIPTION :
  --   Computes the cone of tropisms defined by the kernel of A,
  --   the exponent matrix of a binomial system.

  -- ON ENTRY :
  --   A         the columns of A store the exponent vectors of a binomial 
  --             system, i.e.: the number of rows equals the number of 
  --             variables and the number of columns equals the number
  --             of equations of the binomial system.

  -- ON RETURN :
  --   rank      the rank of the cone of tropisms;
  --   V         if r > 0, then columns of V contain generators
  --             for the cone of the tropisms.

  procedure Check_Inner_Products
               ( A : in Standard_Integer_Matrices.Matrix;
                 V : in Standard_Integer_Matrices.Matrix;
                 bug : out boolean );
  procedure Check_Inner_Products
               ( file : in file_type;
                 A : in Standard_Integer_Matrices.Matrix;
                 V : in Standard_Integer_Matrices.Matrix;
                 bug : out boolean );
  procedure Check_Inner_Products
               ( A : in Standard_Integer64_Matrices.Matrix;
                 V : in Standard_Integer64_Matrices.Matrix;
                 bug : out boolean );
  procedure Check_Inner_Products
               ( file : in file_type;
                 A : in Standard_Integer64_Matrices.Matrix;
                 V : in Standard_Integer64_Matrices.Matrix;
                 bug : out boolean );

  -- DESCRIPTION :
  --   Checks whether all inner products of the tropisms in V
  --   with all columns of A equal zero.

  -- ON ENTRY :
  --   file      for intermediate output;
  --   A         input matrix for Cone_of_Tropisms;
  --   V         output matrix of Cone_of_Tropisms;
  --   output    if intermediate output is wanted.

  -- ON RETURN :
  --   bug       true if some inner product was nonzero, false if otherwise.

  procedure Check_Rank
               ( V : in Standard_Integer_Matrices.Matrix;
                 rank : in integer32; bug : out boolean );
  procedure Check_Rank
               ( file : in file_type;
                 V : in Standard_Integer_Matrices.Matrix;
                 rank : in integer32; bug : out boolean );
  procedure Check_Rank
               ( V : in Standard_Integer64_Matrices.Matrix;
                 rank : in integer32; bug : out boolean );
  procedure Check_Rank
               ( file : in file_type;
                 V : in Standard_Integer64_Matrices.Matrix;
                 rank : in integer32; bug : out boolean );

  -- DESCRIPTION :
  --   Checks whether all the columns of V are linearly independent
  --   and if the dimension of V is as predicted by rank.

  -- ON ENTRY :
  --   file      for intermediate output;
  --   V         output matrix of Cone_of_Tropisms;
  --   rank      computed by Cone_of_Tropisms;
  --   output    if intermediate output is wanted.

  -- ON RETURN :
  --   bug       true if some inner product was nonzero, false if otherwise.

  procedure Check_Cone 
               ( A : in Standard_Integer_Matrices.Matrix;
                 V : in Standard_Integer_Matrices.Matrix;
                 rank : in integer32; bug : out boolean );
  procedure Check_Cone 
               ( file : in file_type;
                 A : in Standard_Integer_Matrices.Matrix;
                 V : in Standard_Integer_Matrices.Matrix;
                 rank : in integer32; bug : out boolean );
  procedure Check_Cone 
               ( A : in Standard_Integer64_Matrices.Matrix;
                 V : in Standard_Integer64_Matrices.Matrix;
                 rank : in integer32; bug : out boolean );
  procedure Check_Cone 
               ( file : in file_type;
                 A : in Standard_Integer64_Matrices.Matrix;
                 V : in Standard_Integer64_Matrices.Matrix;
                 rank : in integer32; bug : out boolean );

  -- DESCRIPTION :
  --   Given in V the generators of a cone of tropisms for the binomial
  --   system defined in A, this procedure checks if all tropisms lie in
  --   the kernel of A and if they are linearly independent.

  procedure Expected_Dimension
               ( A : in Standard_Integer_Matrices.Matrix;
                 rnk : in integer32; dim : out integer32 );
  procedure Expected_Dimension
               ( file : in file_type;
                 A,V : in Standard_Integer_Matrices.Matrix;
                 rnk : in integer32; dim : out integer32 );
  procedure Expected_Dimension
               ( A : in Standard_Integer64_Matrices.Matrix;
                 rnk : in integer32; dim : out integer32 );
  procedure Expected_Dimension
               ( file : in file_type;
                 A,V : in Standard_Integer64_Matrices.Matrix;
                 rnk : in integer32; dim : out integer32 );

  -- DESCRIPTION :
  --   Determines the expected dimension of the solution set of the
  --   binomial system with exponent matrix A and tropisms in V.

  -- ON ENTRY :
  --   file      for intermediate output;
  --   A         exponent matrix of a binomial system;
  --   V         generators of the cone of tropisms;
  --   rnk       rank of the exponent matrix A.

  -- ON RETURN :
  --   d         expected dimension of the solution set of the binomial
  --             system with exponent matrix in A and tropisms in V.

-- STAGE II : tropisms define a unimodular transformation,
--   this stage is refactored into standard_exponent_transformations

-- STAGE III : the coefficients of the representation

  procedure Upper_Transformed_Exponents
               ( A,M : in Standard_Integer_Matrices.Matrix;
                 dim : in integer32;
                 U : out Standard_Integer_Matrices.Matrix;
                 rank : out integer32;
                 pivots : out Standard_Integer_Vectors.Vector );
  procedure Upper_Transformed_Exponents
               ( file : in file_type;
                 A,M : in Standard_Integer_Matrices.Matrix;
                 dim : in integer32;
                 U : out Standard_Integer_Matrices.Matrix;
                 rank : out integer32;
                 pivots : out Standard_Integer_Vectors.Vector );
  procedure Upper_Transformed_Exponents
               ( A,M : in Standard_Integer64_Matrices.Matrix;
                 dim : in integer32;
                 U : out Standard_Integer64_Matrices.Matrix;
                 rank : out integer32;
                 pivots : out Standard_Integer_Vectors.Vector );
  procedure Upper_Transformed_Exponents
               ( file : in file_type;
                 A,M : in Standard_Integer64_Matrices.Matrix;
                 dim : in integer32;
                 U : out Standard_Integer64_Matrices.Matrix;
                 rank : out integer32;
                 pivots : out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Upper triangulates the product M*A returning rank and pivots
  --   to define the exponent matrix for the equations in the binomial
  --   system that will give the coefficients in the representation of
  --   the binomial variety.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   A         the exponent matrix of the binomial system;
  --   M         output of Unimodular_Coordinate_Transformation;
  --   dim       the expected dimension and the number of zero rows
  --             at the beginning of M*A.

  -- ON RETURN :
  --   U         an equivalent upper triangular matrix of A,
  --             the number of rows of U equals the number of rows of M
  --             minus the expected dimension,
  --             i.e.: U'range(1) = 1..M'last(1)-dim;
  --   rank      the rank of the matrix U;
  --   pivots    the pivots in the upper triangular matrix U.

  procedure Nonpivot_Selection
               ( A : in Standard_Integer_Matrices.Matrix;
                 p : in Standard_Integer_Vectors.Vector;
                 B : out Standard_Integer_Matrices.Matrix );
  procedure Nonpivot_Selection
               ( file : in file_type;
                 A : in Standard_Integer_Matrices.Matrix;
                 p : in Standard_Integer_Vectors.Vector;
                 B : out Standard_Integer_Matrices.Matrix );

  -- DESCRIPTION :
  --   Returns in B the same exponents as in A, except for those
  --   variables that occur in the pivot column information of p.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics;
  --   A         the exponent matrix of the binomial system;
  --   p         pivots of a rational coordinate transformation,
  --             computed by Pivots_of_Rational_Coordinate_Transformation.

  -- ON RETURN :
  --   B         the exponent matrix of the same binomial system x^A,
  --             but with variables indexed by p removed.

  function Product_of_Pivots
               ( U : Standard_Integer_Matrices.Matrix;
                 pivots : Standard_Integer_Vectors.Vector ) return integer32;
  function Product_of_Pivots
               ( U : Standard_Integer64_Matrices.Matrix;
                 pivots : Standard_Integer_Vectors.Vector ) return integer64;

  -- DESCRIPTION :
  --   Returns the product of the pivots of the matrix U.

  procedure Solve_Leading_Coefficient_System
               ( n,d : in integer32;
                 A : in Standard_Integer_Matrices.Matrix;
                 pivots : Standard_Integer_Vectors.Vector;
                 c : in Standard_Complex_Vectors.Vector;
                 sols : out Solution_List );
  procedure Solve_Leading_Coefficient_System
               ( file : in file_type; n,d : in integer32;
                 A : in Standard_Integer_Matrices.Matrix;
                 pivots : Standard_Integer_Vectors.Vector;
                 c : in Standard_Complex_Vectors.Vector;
                 sols : out Solution_List );

  -- DESCRIPTION :
  --   Solves the binomial system defined by the exponent matrix A
  --   with given pivots and the coefficients c.  This binomial system
  --   is obtained after an integer valued coordinate transformation.
  --   If file is provided, then extra output is written.

  procedure Solve_Projected_Coefficient_System 
               ( n : in integer32;
                 A : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector;
                 sols : out Solution_List );
  procedure Solve_Projected_Coefficient_System 
               ( file : in file_type; n : in integer32;
                 A : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector;
                 sols : out Solution_List );

  -- DESCRIPTION :
  --   Solves the binomial system defined by the expononent matrix A
  --   and coefficients c.  This binomial system is obtained after
  --   projecting onto the nonpivot variables for a general rational
  --   coordinate transformation.
  --   If file is provided, then extra output is written.

-- STAGE IV : testing and using the algebraic set

  function Evaluate_Tropisms
               ( T : Standard_Integer_Matrices.Matrix;
                 z : Standard_Complex_Vectors.Vector )
               return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns z^T, the algebraic set represented by tropisms in T
  --   at the values z for the parameters.

  -- REQUIRED : T'range(1) = 1..n and T'range(2) = z'range = 1..d,
  --   where n is the ambient dimension and d the dimension of the
  --   algebraic set defined by the binomial system.

  function Transform_Coefficients
               ( d : integer32; M : Standard_Integer_Matrices.Matrix;
                 c : Standard_Complex_Vectors.Vector )
               return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Applies the unimodular transformation M to the coefficients c,
  --   where the first d rows of M contain the tropisms.

  function Evaluate_Algebraic_Set
               ( M : Standard_Integer_Matrices.Matrix;
                 c,z : Standard_Complex_Vectors.Vector )
               return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the algebraic set at the values z for the parameters,
  --   using the formulas (c,z)^M, where the unimodular transformation M
  --   contains in its first d rows the tropisms.

  function Evaluate_Binomial_System
               ( A : Standard_Integer_Matrices.Matrix;
                 b,x : Standard_Complex_Vectors.Vector )
               return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns x^A - b, the binomial system defined by the exponent
  --   matrix in A and right hand side vector b, evaluated at x.
 
  -- REQUIRED : A'range(2) = b'range and A'range(1) = x'range.

  procedure Residual_Test
               ( A,M : in Standard_Integer_Matrices.Matrix;
                 b,c,z : in Standard_Complex_Vectors.Vector;
                 y : out Standard_Complex_Vectors.Vector;
                 r : out double_float );

  -- DESCRIPTION :
  --   Evaluates the solution curve defined by the tropisms in M
  --   the coefficients in c at the point x in the system x^A - b = y.
  --   If M and c are right, then the residual in r should be zero.

  procedure Random_Point_Filter
               ( A : in Standard_Integer_Matrices.Matrix;
                 b : in Standard_Complex_Vectors.Vector;
                 d : in integer32; M : in Standard_Integer_Matrices.Matrix;
                 c : in Solution_List; tol : in double_float;
                 filtered : out Solution_List );
  procedure Random_Point_Filter
               ( file : in file_type;
                 A : in Standard_Integer_Matrices.Matrix;
                 b : in Standard_Complex_Vectors.Vector;
                 d : in integer32; M : in Standard_Integer_Matrices.Matrix;
                 c : in Solution_List; tol : in double_float;
                 filtered : out Solution_List );

  -- DESCRIPTION :
  --   Generates random values for the d parameters and evaluates the 
  --   solution defined by the tropisms in M and the coefficients in c 
  --   at the binomial system x^A - b.  Coefficients for which the
  --   max norm of the residual is larger than the tolerance tol
  --   are filtered out and not in the filtered list on return.
 
-- THE SOLVERS :

  procedure Solve
               ( d : in integer32;
                 A,V : in Standard_Integer_Matrices.Matrix;
                 b : in Standard_Complex_Vectors.Vector;
                 tol : in double_float;
                 M : out Standard_Integer_Matrices.Matrix;
                 w : out Standard_Integer_Vectors.Vector;
                 c : out Solution_List );
  procedure Solve
               ( file : in file_type; d : in integer32;
                 A,V : in Standard_Integer_Matrices.Matrix;
                 b : in Standard_Complex_Vectors.Vector;
                 tol : in double_float;
                 M : out Standard_Integer_Matrices.Matrix;
                 w : out Standard_Integer_Vectors.Vector;
                 c : out Solution_List );

  -- DESCRIPTION :
  --   Given generators for a cone of tropisms in V,
  --   and with expected dimension d, computes solutions
  --   to the binomial system defined by x^A = c.
  
  -- ON ENTRY : 
  --   file      for output during tests;
  --   d         expected dimension of the solution set;
  --   A         exponent matrix of the binomial system x^A - b;
  --   V         columns contain the generators of cone of the tropisms;
  --   b         constants of the binomial system x^A - b;
  --   tol       tolerance to filter solutions for the coefficients;
 
  -- REQUIRED : d > 0.
  --
  -- ON RETURN :
  --   M         unimodular transformation containing the tropisms;
  --   w         denominators in the rational parameterization,
  --             if w(w'first) equals zero, then only M determines
  --             the parameterization which is then integer,
  --             otherwise the parameterization is determined by V,
  --             M and the denominators in w;
  --   c         coefficients of the solutions.

  procedure Black_Box_Solver
              ( p : in Laur_Sys;
                fail : out boolean; d : out integer32;
                M : out Standard_Integer_Matrices.Link_to_Matrix;
                c : out Solution_List );
  procedure Black_Box_Solver
              ( file : in file_type; p : in Laur_Sys;
                fail : out boolean; d : out integer32;
                M : out Standard_Integer_Matrices.Link_to_Matrix;
                c : out Solution_List );

  -- DESCRIPTION :
  --   Black box solver for a binomial system.

  -- ON ENTRY :
  --   file     for extra checks and intermediate output;
  --   p        a Laurent polynomial system.

  -- ON RETURN :
  --   fail     true if p is not a binomial system;
  --   d        the expected dimension of the solution set;
  --   M        unimodular transformation if it exists;
  --   c        coefficients of the solutions.

-- COMPUTING THE DEGREE OF A SOLUTION :

  function Support ( T : Standard_Integer_Matrices.Matrix )
                   return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION :
  --   Stores the tropisms in the rows of T into a list,
  --   adding the origin.  This is the support of a generic linear
  --   system where the ambient variables have been replaced by the
  --   representation defined by the tropisms in T.

  function Volume ( A : Lists_of_Integer_Vectors.List ) return natural32;

  -- DESCRIPTION :
  --   Applies the dynamic lifting algorithm to compute the volume
  --   of the polytope spanned by the points in A.

  function Degree ( T : Standard_Integer_Matrices.Matrix ) return natural32;

  -- DESCRIPTION :
  --   Returns the degree of the solution set defined by the tropisms
  --   in the matrix T.

end Standard_Binomial_Varieties;
