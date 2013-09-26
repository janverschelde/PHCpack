with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Integer64_Vectors;
with Standard_Integer64_Matrices;        use Standard_Integer64_Matrices;
with Lists_of_Integer_Vectors;
with Standard_Lattice_3d_Facets;
with Generic_Lists;

package Standard_Lattice_Facets is

-- DESCRIPTION :
--   A lattice polytope is spanned by points with integer coordinates.
--   This package provides an implementation of the gift wrapping method
--   for any dimension for points with coordinates of modest size,
--   because we work with standard integer 64-bit arithmetic.
--   We assume here that the dimension is at least 4.

  type Facet;
  type Link_to_Facet is access Facet;
  type Array_of_Facets is array ( integer32 range <> ) of Link_to_Facet;

  type Facet ( d,n,m : integer32 ) is record
    label : integer32;
    normal : Standard_Integer64_Vectors.Vector(1..d);
    points : Standard_Integer_Vectors.Vector(1..n);
    neighbors : Array_of_Facets(1..m);
    case d is
      when 4 => ridges : Standard_Lattice_3d_Facets.Array_of_3d_Facets(1..m);
      when others
             => faces : Array_of_Facets(1..m);
    end case;
  end record;

  package Lists_of_Facets is new Generic_Lists(Link_to_Facet);
  type Facet_List is new Lists_of_Facets.List;

-- CONSTRUCTORS :

  procedure Convex_Hull_of_Ridge
              ( A : in Matrix; v : in Standard_Integer64_Vectors.Vector;
                s : in Standard_Integer_Vectors.Vector;
                p : out integer32; M : out Matrix; F : out Facet_List );

  -- DESCRIPTION :
  --   Returns the list of facets of the ridge of A supported by v.

  -- ON ENTRY :
  --   A        matrix of more than 4 rows;
  --   v        normal to a facet of convex hull spanned by A;
  --   s        indices to points of A supported by v,
  --            i.e.: s = Standard_Lattice_Supports.Support(A,v).

  -- ON RETURN :
  --   p        pivot of v, p = Standard_Power_Transformations.Pivot(v);
  --   M        matrix to turn facet perpendicular to coordinate axes;
  --   F        list of 3d faces of convex hull of ridge.

  function Filter_non_Vertices
              ( f : Facet; v : Standard_Integer_Vectors.Vector ) return Facet;

  -- DESCRIPTION :
  --   Given in v the vertex points for the facet f, if f.points /= v,
  --   then the facet on return has its points correctly adjusted.

  function Faces_of_Facet
              ( A : Matrix; v : Standard_Integer64_Vectors.Vector )
              return Facet;

  -- DESCRIPTION :
  --   Computes the ridges of a facet of the polytope spanned by the points
  --   with coordinates in the columns of A, supported in the direction v,
  --   in the general case where A'last(1) > 4.

  function Ridges_of_Facet
              ( A : Matrix; v : Standard_Integer64_Vectors.Vector )
              return Facet;

  -- DESCRIPTION :
  --   Computes the ridges of a facet of the polytope spanned by the points
  --   with coordinates in the columns of A, supported in the direction v,
  --   in all cases where A'last(1) >= 4.

  procedure Neighbors ( A : in Matrix; f : in Link_to_Facet;
                        idcnt : in out integer32 );

  -- DESCRIPTION :
  --   Performs giftwrapping for every ridge of the facet f
  --   to compute facets of conv(A) that are neighbors to f.
  --   The current highest identification number in idcnt gets updated.

  function Convex_Hull
               ( A : Matrix; v : Standard_Integer64_Vectors.Vector )
               return Facet_List;

  -- DESCRIPTION :
  --   Returns a list of facets supporting the convex hull of A.
  --   The vector v is the first facet normal.

  function Convex_Hull ( A : Matrix ) return Facet_List;

  -- DESCRIPTION :
  --   Returns a list of facets supporting the convex hull of A.

-- SELECTORS :

  function Is_Connected ( f : Link_to_Facet ) return boolean;

  -- DESCRIPTION :
  --   Returns true if all neighbors of f are defined.

  function Pop ( f : Facet_List ) return Link_to_Facet;

  -- DESCRIPTION :
  --   Returns the first facet in the list for which Is_Connected()
  --   returns false.  If all facets in f are connected, then the
  --   facets in the list define the convex hull of a point set.

  function Vertices ( n : integer32; f : Facet_List )
                    return Standard_Integer_Vectors.Vector;
  function Vertices ( n : integer32; f : Array_of_Facets )
                    return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   The vector on return contains labels to the vertices of a polytope
  --   spanned by n points with facet in f.

  function Edges ( f : Facet_List )
                 return Lists_of_Integer_Vectors.List;
  function Edges ( f : Array_of_Facets )
                 return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION :
  --   Returns the list of points that span the edges of the faces in f.

  function Ridges ( f : Facet_List )
                  return Lists_of_Integer_Vectors.List;
  function Ridges ( f : Array_of_Facets )
                  return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION :
  --   Returns the list of points that span the ridges of the faces in f.

  function Faces ( f : Facet_List; d : integer32 )
                 return Lists_of_Integer_Vectors.List;
  function Faces ( f : Array_of_Facets; d : integer32 )
                 return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION :
  --   Returns the list of points that span the d-faces of the facets in f.

  function fvector ( n,d : integer32; f : Facet_List )
                   return Standard_Natural_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the f-vector of a d-dimensional polytope, spanned by n points.

-- DESTRUCTORS :

  procedure Clear ( f : in out Facet );
  procedure Clear ( f : in out Link_to_Facet );
  procedure Clear ( f : in out Array_of_Facets );

  -- DESCRIPTION :
  --   Deallocation of all memory occupied by the facets.

end Standard_Lattice_Facets;
