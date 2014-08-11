with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;          use DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;          use DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;         use DoblDobl_Complex_Matrices;

package DoblDobl_Point_Coordinates is

-- DESCRIPTION :
--   A point in a k-plane has a coordinate vector of length k,
--   which is obtained by projection.  The coordinate vector in
--   the ambient n-space is obtained by expansion of the short
--   coordinate k-vector using the generators of the plane.
--   All floating-point operations are done in double double precision.

  function Affine_Coordinates ( x : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns a vector y of range 1..x'last where y(i) = x(i)/x(0).
  -- REQUIRED : x'first = 0 and x(0) /= 0.

  function Projective_Coordinates ( x : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns a vector y with y(0) = 1 and y(x'range) = x.
  -- REQUIRED : x'first = 1.

  procedure Max_Norm ( x : in Vector;
                       ind : out integer32; nrm : out double_double );

  -- DESCRIPTION :
  --   On return, AbsVal(x(ind)) = nrm is the maximum norm of x.

  procedure Scale ( x : in out Vector; ind : in integer32 );

  -- DESCRIPTION :
  --   Divides x by x(ind).

  function Affine_Expand ( c : Complex_Number; b,v : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns b + c*v, the point on the line with offset b and direction v.

  function Affine_Expand ( c,b : Vector; v : VecVec ) return Vector;

  -- DESCRIPTION :
  --   Given a linear space spanned by the directions in v and offset b,
  --   the coordinates c are expanded to the ambient space.

  function Affine_Expand ( c : Vector; p : Matrix ) return Vector;

  -- DESCRIPTION :
  --   Expands the coordinates c of the point in the plane with offset
  --   the 0-th column of p and the directions in the next columns.
  --   The range of the vector on return is p'range(1).
  -- REQUIRED : c'range = 1..p'last(2).

  function Projective_Expand ( c : Vector; p : Matrix ) return Vector;

  -- DESCRIPTION :
  --   Expands the homogeneous coordinates c of the point in the plane p.
  --   The range of the vector on return is p'range(1).
  -- REQUIRED : c'range = p'range(2).

  function Inner_Product ( u,v : Vector ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the complex inner product of the two vectors u and v.

  function Project ( x,b : Vector; v : VecVec ) return Vector;

  -- DESCRIPTION :
  --   Projects the point x onto the linear space spanned by b and v.
  --   For all x = Expand(c,b,v), we have c = Project(x,b,v).
  -- REQUIRED : vectors in v form an orthonormal basis.

  function Project ( x : Vector; p : Matrix ) return Vector;

  -- DESCRIPTION :
  --   Projects the point x onto the plane represented by the matrix.
  --   For all x = Expand(c,p), we have c = Project(x,p).
  -- REQUIRED : the vectors in columns 1 to p'last(2) are orthonormal.

end DoblDobl_Point_Coordinates;
