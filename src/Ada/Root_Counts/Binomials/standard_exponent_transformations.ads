with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Standard_Integer64_Matrices;

package Standard_Exponent_Transformations is

-- DESCRIPTION :
--   The procedures in this package construct unimodular matrices defined
--   by cones of tropisms, with 32-bit and 64-bit integer arithmetic.
--   This is the second stage in the solver of standard_binomial_varieties.
--   The operations in this package generalize the power transformations.

  procedure UV_Coordinate_Transformation
              ( U,V : in Standard_Integer_Matrices.Matrix;
                M : out Standard_Integer_Matrices.Matrix );
  procedure UV_Coordinate_Transformation
              ( U,V : in Standard_Integer64_Matrices.Matrix;
                M : out Standard_Integer64_Matrices.Matrix );

  -- DESCRIPTION :
  --   Computes the unimodular coordinate transformation M by
  --   multiplication of an extended inverse of U with the inverse of V.

  -- REQUIRED : U and V are unimodular matrices.

  function Diagonal_Product 
               ( A : Standard_Integer_Matrices.Matrix ) return integer32;
  function Diagonal_Product 
               ( A : Standard_Integer64_Matrices.Matrix ) return integer64;

  -- DESCRIPTION :
  --   Returns the product of the elements on the diagonal of A.

  procedure Unimodular_Coordinate_Transformation
               ( C : in Standard_Integer_Matrices.Matrix;
                 U,T,V,M : out Standard_Integer_Matrices.Matrix;
                 p : out integer32; uv,fail : out boolean );
  procedure Unimodular_Coordinate_Transformation
               ( C : in Standard_Integer64_Matrices.Matrix;
                 U,T,V,M : out Standard_Integer64_Matrices.Matrix;
                 p : out integer64; uv,fail : out boolean );

  -- DESCRIPTION :
  --   Computes a unimodular transformation M based on the tropisms in C.

  -- ON ENTRY :
  --   C         a cone of tropisms, generators are in the columns of C.

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
               ( C : in Standard_Integer_Matrices.Matrix;
                 fail : out boolean;
                 M : out Standard_Integer_Matrices.Matrix );
  procedure Unimodular_Coordinate_Transformation
               ( file : in file_type;
                 C : in Standard_Integer_Matrices.Matrix;
                 fail : out boolean;
                 M : out Standard_Integer_Matrices.Matrix );
  procedure Unimodular_Coordinate_Transformation
               ( C : in Standard_Integer64_Matrices.Matrix;
                 fail : out boolean;
                 M : out Standard_Integer64_Matrices.Matrix );
  procedure Unimodular_Coordinate_Transformation
               ( file : in file_type;
                 C : in Standard_Integer64_Matrices.Matrix;
                 fail : out boolean;
                 M : out Standard_Integer64_Matrices.Matrix );

  -- DESCRIPTION :
  --   Wrapper of the previous call, with extra output if file is provided.

  procedure Pivots_of_Rational_Coordinate_Transformation
               ( C : in Standard_Integer_Matrices.Matrix;
                 H : out Standard_Integer_Matrices.Matrix;
                 p : out Standard_Integer_Vectors.Vector );
  procedure Pivots_of_Rational_Coordinate_Transformation
               ( file : in file_type;
                 C : in Standard_Integer_Matrices.Matrix;
                 H : out Standard_Integer_Matrices.Matrix;
                 p : out Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   A rational coordinate transformation for a k-dimensional set
  --   is defined by k vectors that span the cone of tropisms
  --   multiplied by the inverse of the diagonal matrix determined
  --   by the pivots of the Hermite normal form of the matrix of tropisms.
  --   This product is padded with standard basis vectors to form a
  --   square matrix.

  -- ON ENTRY :
  --   file      for intermediate output and diagnostics
  --   C         a cone of tropisms, computed as V in Cone_of_Tropisms.

  -- ON RETURN :
  --   H         Hermite normal form of the transpose of C;
  --   p         pivot columns for the first nonzero element
  --             on the rows of H.

  function Is_Zero_Row
               ( A : Standard_Integer_Matrices.Matrix; i : integer32 )
               return boolean;
  function Is_Zero_Row
               ( A : Standard_Integer64_Matrices.Matrix; i : integer32 )
               return boolean;

  -- DESCRIPTION :
  --   Returns true if all elements on the i-th row of A are zero,
  --   returns false if otherwise.

  procedure Test_Unimodular_Coordinate_Transformation
               ( A,M : in Standard_Integer_Matrices.Matrix;
                 dim : in integer32; nzr : out natural32; fail : out boolean );
  procedure Test_Unimodular_Coordinate_Transformation
               ( file : in file_type;
                 A,M : in Standard_Integer_Matrices.Matrix;
                 dim : in integer32; nzr : out natural32; fail : out boolean );
  procedure Test_Unimodular_Coordinate_Transformation
               ( A,M : in Standard_Integer64_Matrices.Matrix;
                 dim : in integer32; nzr : out natural32; fail : out boolean );
  procedure Test_Unimodular_Coordinate_Transformation
               ( file : in file_type;
                 A,M : in Standard_Integer64_Matrices.Matrix;
                 dim : in integer32; nzr : out natural32; fail : out boolean );

  -- DESCRIPTION :
  --   The product M*A should have at least as many zero rows as
  --   the expected dimension of the binomial variety.

  -- ON ENTRY :
  --   file      for intermediate output;
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

  function Rational_Coordinate_Transformation
               ( V : in Standard_Integer_Matrices.Matrix;
                 pivots : in Standard_Integer_Vectors.Vector )
               return Standard_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the nominators of a rational coordinate transformation,
  --   extending the matrix V with tropisms in its columns with rows
  --   of standard basis vectors, taking into account the pivots of
  --   the Hermite normal form of the transpose of V.

  procedure Test_Rational_Coordinate_Transformation
              ( file : in file_type;
                M : in Standard_Integer_Matrices.Matrix;
                w : in Standard_Integer_Vectors.Vector;
                tol : in double_float; fail : out boolean );

  -- DESCRIPTION :
  --   Tests whether the determinant of the rational coordinate transformation
  --   defined by M and w is within the given tolerance equal to +1  or -1.

  -- ON ENTRY :
  --   file     for intermediate results and writing of the determinant;
  --   M        nominators of the rational coordinate transformation;
  --   w        denominators for the rows of M that contain the tropisms;
  --   tol      tolerance for use in the determinant test.

  -- ON RETURN :
  --   fail     true if the determinant differs more than tol
  --            from +1 or -1, false otherwise.

end Standard_Exponent_Transformations;
