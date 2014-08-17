with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;         use QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Polynomials;      use QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;     use QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;        use QuadDobl_Complex_Solutions;

package QuadDobl_Linear_Span is

-- DESCRIPTION :
--   This package allows to compute the span of a solution component,
--   using quad double complex arithmetic.

  function Create ( sols : Solution_List ) return Matrix;

  -- DESCRIPTION :
  --   The matrix on return has in its rows the differences of all
  --   vectors in sols with the first solution vector.

  -- REQUIRED : Length_Of(sols) >= 2.

  procedure Rank ( v : in out Matrix;
                   tol : in double_float; r : out natural32 );

  -- DESCRIPTION :
  --   Triangulates the matrix v to determine the rank.

  -- ON ENTRY :
  --   v        configuration of vectors in rows of matrix;
  --   tol      tolerance to decide if number is zero or not.

  -- ON RETURN :
  --   v        triangulated matrix;
  --   r        rank of the vector configuration.

  function Pivots ( v : Matrix; tol : double_float; r : natural32 )
                  return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the column indices of the nonzero elements in the upper
  --   triangular matrix used to determine the rank.

  -- ON ENTRY :
  --   v        output of procedure Rank from above;
  --   tol      tolerance to decide whether number is zero;
  --   r        rank of the matrix v, obtained from above.

  -- ON RETURN :
  --   a vector of range 1..r with pivots to compute kernel of v.

  function Kernel ( v : Matrix; tol : double_float; r : natural32;
                    pivots : Standard_Integer_Vectors.Vector;
                    point : QuadDobl_Complex_Vectors.Vector ) return Matrix;

  -- DESCRIPTION :
  --   Returns n-r equations (n = ambient dimension) which describe the
  --   linear span of the vectors in the upper triangular matrix v and
  --   the given point.

  -- ON ENTRY :
  --   v        output of the procedure Rank from above;
  --   tol      tolerance to decide whether a number is zero;
  --   r        rank of the matrix v, obtained from above;
  --   pivots   free variables, obtained from function Pivots.

  -- ON RETURN :
  --   matrix of row range 1..n-r and column range 0..n with the
  --   coefficients of the hyperplanes defining the linear span.

  function Equations ( m : Matrix ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the hyperplane equations in m in polynomial format.

  function Eliminators ( kernel_eqs : Matrix;
                         pivots : Standard_Integer_Vectors.Vector )
                       return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the polynomials in function of the pivot variables
  --   to be used to eliminate the nonpivot variables.

  function Eliminate_non_Pivots
             ( p : Poly_Sys; pivots : Standard_Integer_Vectors.Vector;
               elm : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   Eliminates the non-pivot variables from the polynomial system p,
  --   using the eliminator equations in elm.

  function Filter ( p : Poly; tol : double_float ) return Poly;

  -- DESCRIPTION :
  --   Only terms with coefficient whose absolute value is higher
  --   than tol are copied to the output.

  function Filter ( p : Poly_Sys; tol : double_float ) return Poly_Sys;

  -- DESCRIPTION :
  --   The system on return has no terms whose coefficient is less
  --   than tol in absolute value.  Also zero polynomials are eliminated.

  function Eliminate_non_Pivots
             ( v : QuadDobl_Complex_Vectors.Vector;
               pivots : Standard_Integer_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;
  function Eliminate_non_Pivots
             ( s : Solution; pivots : Standard_Integer_Vectors.Vector )
             return Solution;
  function Eliminate_non_Pivots
             ( sols : Solution_List;
               pivots : Standard_Integer_Vectors.Vector )
             return Solution_List;

  -- DESCRIPTION :
  --   Eliminates all non pivot entries from the solution vectors.
 
end QuadDobl_Linear_Span;
