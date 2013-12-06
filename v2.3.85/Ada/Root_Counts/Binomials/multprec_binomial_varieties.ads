with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Integer_Numbers;
with Standard_Integer_Vectors;
with Multprec_Integer_Matrices;

package Multprec_Binomial_Varieties is

-- DESCRIPTION :
--   Tools are offered to manipulate algebraic sets defined by binomial 
--   systems, using multiprecision integer arithmetic on the exponents.

-- STAGE I : computing and checking a cone of tropisms

  procedure Cone_of_Tropisms 
               ( A : in Multprec_Integer_Matrices.Matrix;
                 rank : out integer32;
                 V : out Multprec_Integer_Matrices.Link_to_Matrix );

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
               ( A : in Multprec_Integer_Matrices.Matrix;
                 V : in Multprec_Integer_Matrices.Matrix;
                 output : in boolean; bug : out boolean );

  -- DESCRIPTION :
  --   Checks whether all inner products of the tropisms in V
  --   with all columns of A equal zero.

  -- ON ENTRY :
  --   A         input matrix for Cone_of_Tropisms;
  --   V         output matrix of Cone_of_Tropisms;
  --   output    if intermediate output is wanted.

  -- ON RETURN :
  --   bug       true if some inner product was nonzero, false if otherwise.

  procedure Check_Rank
               ( V : in Multprec_Integer_Matrices.Matrix;
                 rank : in integer32; output : in boolean;
                 bug : out boolean );

  -- DESCRIPTION :
  --   Checks whether all the columns of V are linearly independent
  --   and if the dimension of V is as predicted by rank.

  -- ON ENTRY :
  --   V         output matrix of Cone_of_Tropisms;
  --   rank      computed by Cone_of_Tropisms;
  --   output    if intermediate output is wanted.

  -- ON RETURN :
  --   bug       true if some inner product was nonzero, false if otherwise.

  procedure Check_Cone 
               ( A : in Multprec_Integer_Matrices.Matrix;
                 V : in Multprec_Integer_Matrices.Matrix;
                 rank : in integer32; output : in boolean;
                 bug : out boolean );

  -- DESCRIPTION :
  --   Given in V the generators of a cone of tropisms for the binomial
  --   system defined in A, this procedure checks if all tropisms lie in
  --   the kernel of A and if they are linearly independent.

  procedure Expected_Dimension
               ( A,V : in Multprec_Integer_Matrices.Matrix;
                 rnk : in integer32; output : in boolean;
                 dim : out integer32 );

  -- DESCRIPTION :
  --   Determines the expected dimension of the solution set of the
  --   binomial system with exponent matrix A and tropisms in V.

  -- ON ENTRY :
  --   A         exponent matrix of a binomial system;
  --   V         generators of the cone of tropisms;
  --   rnk       rank of the exponent matrix A;
  --   output    if true, then additional output is written.

  -- ON RETURN :
  --   d         expected dimension of the solution set of the binomial
  --             system with exponent matrix in A and tropisms in V.

-- STAGE II : tropisms define a unimodular transformation

  procedure UV_Coordinate_Transformation
              ( U,V : in Multprec_Integer_Matrices.Matrix;
                M : out Multprec_Integer_Matrices.Matrix );

  -- DESCRIPTION :
  --   Computes the unimodular coordinate transformation M by
  --   multiplication of an extended inverse of U with the inverse of V.

  -- REQUIRED : U and V are unimodular matrices.

  function Diagonal_Product 
               ( A : Multprec_Integer_Matrices.Matrix )
               return Multprec_Integer_Numbers.Integer_Number;

  -- DESCRIPTION :
  --   Returns the product of the elements on the diagonal of A.

  procedure Unimodular_Coordinate_Transformation
               ( C : in Multprec_Integer_Matrices.Matrix;
                 U,T,V,M : out Multprec_Integer_Matrices.Matrix;
                 p : out Multprec_Integer_Numbers.Integer_Number;
                 uv,fail : out boolean );

  -- DESCRIPTION :
  --   Computes a unimodular transformation M based on the tropisms in C.

  -- ON ENTRY :
  --   C         a cone of tropisms, computed as V in Cone_of_Tropisms.

  -- ON RETURN :
  --   U         U in the Smith form of the transpose of C;
  --   T         the Smith form of the transpose of C;
  --   V         V in the Smith form of the transpose of C;
  --   M         a unimodular coordinate transformation if not fail;
  --   p         product of element on the diagonal of T;
  --   uv        true if the type of unimodular transformation used
  --             and extension of U inverse an the inverse of V,
  --             false if otherwise;
  --   fail      true if no rational parameterization is possible,
  --             false if M has meaning.

  procedure Unimodular_Coordinate_Transformation
               ( C : in Multprec_Integer_Matrices.Matrix;
                 output : in boolean; fail : out boolean;
                 M : out Multprec_Integer_Matrices.Matrix );

  -- DESCRIPTION :
  --   Wrapper of the previous call, with extra output when output is true.

  function Is_Zero_Row
               ( A : Multprec_Integer_Matrices.Matrix; i : integer32 )
               return boolean;

  -- DESCRIPTION :
  --   Returns true if all elements on the i-th row of A are zero,
  --   returns false if otherwise.

  procedure Test_Unimodular_Coordinate_Transformation
               ( A,M : in Multprec_Integer_Matrices.Matrix;
                 dim : in integer32; output : in boolean;
                 nzr : out natural32; fail : out boolean );

  -- DESCRIPTION :
  --   The product M*A should have at least as many zero rows as
  --   the expected dimension of the binomial variety.

  -- ON ENTRY :
  --   A         the exponent matrix of the binomial system;
  --   M         the unimodular coordinate transformation to bring
  --             the binomial system in normal form and to compute
  --             the leading coefficients of the parameterization
  --             of the solution set of the binomial system;
  --   dim       the expected dimension of the binomial variety;
  --   output    if true, then M*A will be displayed on screen.

  -- ON RETURN :
  --   nzr       the number of nonzero rows in M*A;
  --   fail      true if the first dim rows in M*A are not all zero,
  --             false if otherwise.

-- STAGE III : the coefficients of the representation

  procedure Upper_Transformed_Exponents
               ( A,M : in Multprec_Integer_Matrices.Matrix;
                 dim : in integer32; output : in boolean;
                 U : out Multprec_Integer_Matrices.Matrix;
                 rank : out integer32;
                 pivots : out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Upper triangulates the product M*A returning rank and pivots
  --   to define the exponent matrix for the equations in the binomial
  --   system that will give the coefficients in the representation of
  --   the binomial variety.

  -- ON ENTRY :
  --   A         the exponent matrix of the binomial system;
  --   M         output of Unimodular_Coordinate_Transformation;
  --   dim       the expected dimension and the number of zero rows
  --             at the beginning of M*A;
  --   output    if true, then extra output will be displayed.

  -- ON RETURN :
  --   U         an equivalent upper triangular matrix of A,
  --             the number of rows of U equals the number of rows of M
  --             minus the expected dimension,
  --             i.e.: U'range(1) = 1..M'last(1)-dim;
  --   rank      the rank of the matrix U;
  --   pivots    the pivots in the upper triangular matrix U.

  function Product_of_Pivots
               ( U : Multprec_Integer_Matrices.Matrix;
                 pivots : Standard_Integer_Vectors.Vector ) 
               return Multprec_Integer_Numbers.Integer_Number;

  -- DESCRIPTION :
  --   Returns the product of the pivots of the matrix U.

end Multprec_Binomial_Varieties;
