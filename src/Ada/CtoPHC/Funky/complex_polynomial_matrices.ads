with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Floating_Matrices;
with Standard_Complex_Matrices;

package Complex_Polynomial_Matrices is

-- DESCRIPTION :
--   This package defines matrices of polynomials in one variable
--   with complex floating-point coefficients.

  type Polynomial_Matrix is
    array ( integer32 range <>, integer32 range <> ) 
    of Standard_Complex_Vectors.Link_to_Vector;

  type Link_to_Polynomial_Matrix is access Polynomial_Matrix;

  type Array_of_Polynomial_Matrices is
    array ( integer32 range <> ) of Link_to_Polynomial_Matrix;

  -- REPRESENTATION :
  --   A polynomial in one variable with complex coefficients of degree d
  --   is represented by its coefficient vector: the i-th entry in this 
  --   vector is the coefficient of x^i.

  function Degrees ( pm : Polynomial_Matrix )
                   return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 1..pm'length(1)*pm'length(2),
  --   with the degrees of all polynomials.  If some entry is null,
  --   then -1 is the degree of the corresponding polynomial.
  --   The total number of coefficients in the matrix pm is the sum of
  --   all the integers in Degrees(pm) plus the length of Degrees(pm).

  function Degrees ( apm : Array_of_Polynomial_Matrices )
                   return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector with all degree vectors of the matrices in apm,
  --   stacked one after the other.

  -- REQUIRED :
  --   All matrices in the array must have the same dimensions.

  function Coefficients ( k : integer32; pm : Polynomial_Matrix )
                        return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of size 1..k with all the coefficients of pm.
  --   The k can be computed via Degrees(pm).

  function Coefficients ( k : integer32; apm : Array_of_Polynomial_Matrices )
                        return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   The vector on return has range 1..k and contains all coefficients
  --   of all polynomial matrix in the array apm.

  function Create ( n,m : integer32;
                    d : Standard_Integer_Vectors.Vector;
                    c : Standard_Complex_Vectors.Vector )
                  return Polynomial_Matrix;

  -- DESCRIPTION :
  --   Returns an n-by-m polynomial matrix with degrees in d and
  --   coefficients in c.

  -- REQUIRED :
  --   d'range = 1..n*m and c'range = sum(d) + n*m.

  function Left_Multiply
             ( a : Standard_Complex_Matrices.Matrix; b : Polynomial_Matrix )
             return Polynomial_Matrix;
  function Left_Multiply
             ( a : Standard_Complex_Matrices.Matrix;
               b : Array_of_Polynomial_Matrices )
             return Array_of_Polynomial_Matrices;
  function Left_Multiply
             ( a : Standard_Floating_Matrices.Matrix; b : Polynomial_Matrix )
             return Polynomial_Matrix;
  function Left_Multiply
             ( a : Standard_Floating_Matrices.Matrix;
               b : Array_of_Polynomial_Matrices )
             return Array_of_Polynomial_Matrices;

  -- DESCRIPTION :
  --   returns a*b, multiplying b from the left with a.

  -- REQUIRED : a'range(2) = b'range(1), for b a Polynomial_Matrix.

  function Eval ( pm : Polynomial_Matrix; x : Complex_Number )
                return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the function evaluation of pm at the point x.

  procedure Clear ( pm : in out Polynomial_Matrix );
  procedure Clear ( lpm : in out Link_to_Polynomial_Matrix );

  -- DESCRIPTION :
  --   Release of all allocated memory.

end Complex_Polynomial_Matrices;
