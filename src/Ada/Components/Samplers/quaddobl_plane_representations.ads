with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with QuadDobl_Complex_Vectors;          use QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;          use QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;         use QuadDobl_Complex_Matrices;

package QuadDobl_Plane_Representations is

-- DESCRIPTION :
--   An affine plane in n-space of co-dimension k can be represented
--   (1) extrinsically (or implicitly) by k linear equations, or
--   (2) intrinsically (or explicitly) by a parametric representation.
--   In particular, we use the following conventions:
--   (1) A hyperplane is defined by a linear equation in x, as
--         c(0) + c(1)*x(1) + c(2)*x2 + .. + c(n)*xn = 0,
--       where c is the coefficient vector of range 0..n.
--   (2) The parametric representation consists of an offset vector
--       (also called basis) and n-k direction vectors:
--         x = b + t(1)*v(1) + t(2)*v(2) + .. + t(n-k)*v(n-k),
--       where x, b, and v(i) are n-vectors of range 1..n.
--   In projective space, we work with homogeneous coordinates:
--   (1) c(0)*x(0) + c(1)*x(1) + c(2)*x2 + .. + c(n)*xn = 0; and
--   (2) x = t(0)*b + t(1)*v(1) + t(2)*v(2) + .. + t(n-k)*v(n-k).
--   The representation above lends itself naturally to vectors of vectors.
--   For efficiency of linear algebra, we may prefer to work with matrices.
--   (1) The k linear equations are then represented by a k-by-(n+1) matrix:
--         c(i,0) + c(i,1)*x(1) + c(i,2)*x(2) + .. + c(i,n)*x(n) = 0,
--       for i ranging from 1 to k.
--   (2) The parametric representation is then the orthogonal complement:
--         x(i) = v(i,0) + t(1)*v(i,1) + t(2)*v(i,2) + .. + t(n-k)*v(i,n-k),
--       where v is an n-by-(n-k+1) matrix.  The columns 1 to n-k are
--       perpendicular to the coefficient matrix c.
--   This package offers operations to convert (1) into (2), and vice versa,
--   using quad double complex arithmetic.

-- DATA STRUCTURE CONVERSIONS :

  function Equations_to_Matrix ( c : VecVec ) return Matrix;
  function Equations_to_VecVec ( c : Matrix ) return VecVec;

  -- DESCRIPTION :
  --   The coefficients for the linear equations are converted
  --   from vector of vectors into matrix format, or vice versa.

  function Equations_to_Matrix ( c : VecVec; n : integer32 ) return Matrix;

  -- DESCRIPTION :
  --   Puts the first coefficients (in the range 0..n) of c
  --   in the rows of a matrix; useful for embedded hyperplane.

  function Generators_to_Matrix ( b : Vector; v : VecVec ) return Matrix;

  -- DESCRIPTION :
  --   Returns a matrix with in its zero-th column the vector b,
  --   and in the columns 1 to v'length the vectors in v.
  -- REQUIRED : v'first = 1.

  procedure Generators_to_VecVec
              ( g : in Matrix; b : out Vector; v : out VecVec );

  -- DESCRIPTION :
  --   Extracts the offset vector b from the first column of g,
  --   and puts the other columns of g in v.
  -- REQUIRED : b'range = g'range(1), v'first = 1 and v'last = g'last(2).

-- FROM EQUATIONS TO GENERATORS :

  procedure Generators1 ( hyp : in Vector;
                          basis : out Vector; directions : out VecVec );

  -- DESCRIPTION :
  --   Computes the basis and directions needed to compute a representation
  --   of the plane by its generators, i.e.: as the basis plus a linear
  --   combination of the directions.

  -- REQUIRED : co-dimension k = 1.

  procedure Generators ( n,k : in integer32; hyp : in VecVec;
                         basis : out Vector; directions : out VecVec );

  -- DESCRIPTION :
  --   Computes a basis point and n-k direction vectors as a parametric
  --   representation of the plane.

  -- REQUIRED : 
  --   The hyperplanes form a complete intersection and moreover,
  --   the first k-by-k minor in the coefficient matrix is of full rank.

  function Generators ( hyp : Matrix ) return Matrix;

  -- DESCRIPTION :
  --   Returns a matrix which has as 0-th column an offset vector,
  --   and in the next columns the direction vectors of the plane
  --   defined by the equations in hyp.

  function Orthogonalize ( v : Matrix ) return Matrix;
  function Orthogonalize ( v : VecVec ) return VecVec;

  -- DESCRIPTION :
  --   Returns an orthonormal basis for the space spanned by v.

  procedure Orthogonalize ( v : in out Matrix );

  -- DESCRIPTION :
  --   Makes the directions in v orthogonal.

-- FROM GENERATORS TO EQUATIONS :

  function Equations1 ( basis,direction : Vector ) return VecVec;

  -- DESCRIPTION :
  --   Given a parametric representation of a line,
  --   n-1 linear equations defining the line are returned.

  -- NOTE : suffix 1 is dimension, not co-dimension as in Generators1.

  function Equations ( basis : Vector; directions : VecVec ) return VecVec;

  -- DESCRIPTION :
  --   Returns the equations of a k-dimensional affine plane, whose points
  --   can all be written as the sum of the basis point and a linear 
  --   combination of the direction vectors.

  -- REQUIRED : the direction vectors are linearly independent.

  function Equations ( g : Matrix ) return Matrix;

  -- DESCRIPTION :
  --   Returns coefficients of hyperplanes defining the plane
  --   generated by the columns in g.

end QuadDobl_Plane_Representations;
