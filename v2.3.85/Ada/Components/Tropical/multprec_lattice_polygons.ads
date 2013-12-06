--with Multprec_Integer64_Numbers;         use Multprec_Integer64_Numbers;
--with Multprec_Integer64_Vectors;         use Multprec_Integer64_Vectors;
--with Multprec_Integer64_Matrices;        use Multprec_Integer64_Matrices;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Vectors;           use Multprec_Integer_Vectors;
with Multprec_Integer_Matrices;          use Multprec_Integer_Matrices;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;

package Multprec_Lattice_Polygons is

-- DESCRIPTION :
--   A lattice polygon is spanned by points with integer coordinates,
--   as obtained as the Newton polygon of a Laurent polynomial.
--   The V-representation of a polygon is counterclockwise ordered list
--   of vertex points.  Two consecutive points in the V-representation
--   span an edge of the polygon.  The H-representation consists of list
--   of inner normals, perpendicular to the edges of the polygon.
--   At each edge the inner normal makes its minimal inner product and
--   this minimal value is the constant of the half-plane representation.

-- I. convert set into lexicographically decreasing matrix

  procedure Convert ( s : in List; A : out Matrix );

  -- DESCRIPTION :
  --   Converts the list s into a 2-by-m matrix A, m = Length_Of(s).

  procedure Convert ( A : in Matrix; s : out List );

  -- DESCRIPTION :
  --   Converts the matrix A into a list s.  The elements of s on return
  --   are the columns of the matrix A.

  procedure Lexicographic_Decreasing_Sort ( A : in out Matrix );

  -- DESCRIPTION :
  --   Sorts the columns in A lexicographically, in decreasing order,
  --   using insertion sort of the largest element.

-- II. compute vertex set and inner normals to the edges

  function Convex_Hull_2D ( A : Matrix ) return Matrix;

  -- DESCRIPTION :
  --   Returns the convex hull of the points with coordinates
  --   in the columns of the matrix A.

  -- REQUIRED :
  --   Points in A are sorted lexicographically in decreasing order.

  function Inner_Normals ( V : Matrix ) return Matrix;

  -- DESCRIPTION :
  --   Given in the columns of V an ordered vertex list of the convex hull,
  --   the i-th column of the matrix on return contains the inner normal
  --   to the edge spanned by the i-th and (i+1)-th column of V,
  --   for i < V'last(2).  The last column of the matrix on return contains
  --   the inner normal to the edge spanned by the first and last vertex.
  --   Components of the inner normals are relatively prime.

  function Inner_Normals ( s : List ) return List;

  -- DESCRIPTION :
  --   Wrapper function: Convert(Inner_Normals(Convex_Hull_2D(Convert(s)))).

-- III. hyperplane representations

  function Rank ( A : Matrix; i : integer32; v : Vector ) return Integer_Number;

  -- DESCRIPTION :
  --   The rank of the point with coordinates in the i-th column of A
  --   is the value of the inner product with the vector v.

  function Rank ( A : Matrix; v : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns a vector of range equal to A'range(2) with the i-th entry
  --   the value of the inner product of v with the i-th column of A.

  function Minimum ( A : Matrix; v : Vector ) return Integer_Number;

  -- DESCRIPTION :
  --   Returns the minimum inner product of the vector v
  --   with the columns in the matrix A.

  function Minima ( A,N : Matrix ) return Vector;
   
  -- DESCRIPTION :
  --   Returns a vector of range equal to N'range(2), with the minimal inner
  --   product for each vector in the columns of N with the columns of A.

  function Number_of_Minima ( v : Vector ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of minima in the vector v.

  function Rank ( A,N : Matrix ) return Matrix;

  -- DESCRIPTION :
  --   The matrix on return has as many columns as A and as many rows as N.
  --   The (i,j)-th entry equals the value of the inner product of the
  --   j-th column of A with the i-th column of N.

-- IV. sanity check

  procedure Check ( A,V,N : in Matrix; output : in boolean;
                    bug : out boolean );

  -- DESCRIPTION :
  --   Performs a check on the set of vertices V and inner normals N,
  --   compute from the points with coordinates in A.

  -- ON ENTRY :
  --   A        matrix with in columns coordinates of points in the plane;
  --   V        vertex set, output of Convex_Hull_2D(A), 
  --            where A is lexicographically in descreasing order;
  --   N        set of inner normals to the edges of the polygon.
  --   output   if output is desired, otherwise check is done silently.

  -- ON RETURN :
  --   bug      true if test with hyperplane representations failed,
  --            due either to a bug in the code or overflow because
  --            of limitations of standard machine arithmetic.

end Multprec_Lattice_Polygons;
