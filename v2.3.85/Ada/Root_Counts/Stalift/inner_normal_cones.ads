with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Integer_Faces_of_Polytope;          use Integer_Faces_of_Polytope;

package Inner_Normal_Cones is

-- DESCRIPTION :
--   This package contains facilities for dealing with inner normal cones,
--   that are cones normal to some faces of a given polytope, with apex
--   at the origin.  Especially the normal cones to vertices are of great
--   interest to determine which points contribute to the mixed volume.
--   The primal representation of a normal cone consists in the generators
--   of the cone.  A matrix of inequalities defines the dual representation.

-- CONSTRUCTORS FOR PRIMAL REPRESENTATION :

  function Generators ( L : List; facets : Faces; x : Vector ) return List;

  -- DESCRIPTION :
  --   Returns the generators of the polyhedral cone, which lies in the
  --   dual space of the polytope spanned by l and is enclosed by the facets,
  --   which all contain the vector x.
  --   The list on return contains all inner normals to the facets.

  function Generators ( L : List; x : Vector ) return List;

  -- DESCRIPTION :
  --   Returns the generators of the polyhedral cone to conv(l),
  --   normal at the vector x.

-- CONSTRUCTOR FOR DUAL REPRESENTATION :

  function Inner_Normal_Cone ( L : List; x : Vector ) return Matrix;

  -- DESCRIPTION :
  --   Returns the inequalities that define the normal cone to the vector x,
  --   normal to the polytope spanned by the points in l.
  --   The range of the matrix on return is (x'range,1..Length_Of(l)-1),
  --   which implies that the inequalities are stored columnwise.
  --   The type of the inequalities is >= 0.

  function Included_Normal_Cone ( gx : List; x : Vector ) return Matrix;

  -- DESCRIPTION :
  --   Returns the inequalities that define the inner normal cone to the
  --   vector x, such that these normals are pointing in the interior of
  --   the polytope spanned by l, when their basis is the point x.
  --   The walls of this cone are the facets, the list gx should
  --   contain all inner normal to the facets at the vertex x.

  -- RANGE FORMATS :
  --   The matrix on return has ranges (x'first-1..x'last,1..Length_Of(gx)).
  --   The first row of the matrix contains <gx,x>,
  --   to be used in the right-hand side of the inequalities.

-- PRIMITIVE SELECTORS :

  function Evaluate ( m : Matrix; i : integer32; v : Vector ) return integer32;

  -- DESCRIPTION :
  --   Evaluates the vector in the ith hyperplane of the inequality matrix.

  function Satisfies ( m : Matrix; i : integer32; v : Vector ) return boolean;
  function Strictly_Satisfies ( m : Matrix; i : integer32; v : Vector )
                              return boolean;

  -- DESCRIPTION :
  --   Returns true if the vector v satisfies the ith inequality of the matrix.
  --   With Strictly_ true is only returned if the evaluation yields a
  --   positive result.

  function Satisfies ( m : Matrix; v : Vector ) return boolean;
  function Strictly_Satisfies ( m : Matrix; v : Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the vector v satisfies all inequalities defined by m.

  function Satisfies ( m : Matrix; v : List ) return boolean;
  function Strictly_Satisfies ( m : Matrix; v : List ) return boolean;

  -- DESCRIPTION :
  --   Returns true if all vectors in v satisfy all inequalities defined by m.

-- SECONDARY SELECTORS :

  function Satisfy ( m : Matrix; L : List ) return List;
  function Strictly_Satisfy ( m : Matrix; L : List ) return List;

  -- DESCRIPTION :
  --   Returns a list of points from l that (strictly) satisfy the 
  --   inequalities in the matrix.

  function Contained_in_Cone ( L : List; x : Vector; v : List ) return boolean;
  function Strictly_Contained_in_Cone ( L : List; x : Vector; v : List )
                                      return boolean;

  -- DESCRIPTION :
  --   Returns true when all vectors in v are contained in the normal cone
  --   to the points in l, at the vector x.

  function Contained_in_Cone ( L,v : List ) return boolean;
  function Strictly_Contained_in_Cone ( L,v : List ) return boolean;

  -- DESCRIPTION :
  --   Returns true when there exists a normal cone to the points in l
  --   that entirely contains the points in v.

-- CONCERNING THE UNION OF CONES :

  function In_Union ( v1,v2 : Vector; k1,k2 : Matrix ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the convex combination of v1 and v2 is contained
  --   in the union of the convex cones k1 and k2, given by their
  --   inequality matrix.

  -- REQUIRED : Satisfies(k1,v1) and Satisfies(k2,v2).

  function In_Union ( v1,v2 : List; k1,k2 : Matrix ) return boolean;

  -- DESCRIPTION :
  --   Returns true if every convex combination of a vector from v1
  --   and a vector from v2 is contained in the union of the convex
  --   cones k1 and k2, given by their inequality matrix.

  -- REQUIRED : Satisfies(k1,v1) and Satisfies(k2,v2).

end Inner_Normal_Cones;
