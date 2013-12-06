with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer64_Vectors;
with Standard_Integer64_Matrices;        use Standard_Integer64_Matrices;
with Lists_of_Integer_Vectors;
with Standard_Lattice_3d_Facets;
with Generic_Lists;

package Standard_Lattice_4d_Facets is

-- DESCRIPTION :
--   Lattice polytopes in 4-dimensional space are a jumping board to get
--   a grip on data structures for polytopes in any dimension.
--   Ridges of a 4D polytope are facets of a 3D polytope.

  type Facet_in_4d;
  type Link_to_4d_Facet is access Facet_in_4d;
  type Array_of_4d_Facets is array ( integer32 range <> ) of Link_to_4d_Facet;

  type Facet_in_4d ( d,n,m : integer32 ) is record
    label : integer32;
    normal : Standard_Integer64_Vectors.Vector(1..d);
    points : Standard_Integer_Vectors.Vector(1..n);
    ridges : Standard_Lattice_3d_Facets.Array_of_3d_Facets(1..m);
    neighbors : Array_of_4d_Facets(1..m);
  end record;

  package Lists_of_4d_Facets is new Generic_Lists(Link_to_4d_Facet);
  type Facet_4d_List is new Lists_of_4d_Facets.List;

-- CONSTRUCTORS :

  function Insert_Zero
              ( v : Standard_Integer64_Vectors.Vector; p : integer32 )
              return Standard_Integer64_Vectors.Vector;

  -- DESCRIPTION :
  --   Inserts a zero at position p of the inner normal v.

  procedure Multiply_with_Transpose
              ( v : in out Standard_Integer64_Vectors.Vector;
                M : in Matrix );

  -- DESCRIPTION :
  --   Multiplies v with the transpose of M.

  function Relabel ( p,s : Standard_Integer_Vectors.Vector )
                   return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the relabeling of p (points supported on a ridge),
  --   using the indices of the supported points on the facet.

  procedure Convex_Hull_of_Ridge
              ( A : in Matrix; v : in Standard_Integer64_Vectors.Vector;
                s : in Standard_Integer_Vectors.Vector;
                p : out integer32; M : out Matrix;
                F : out Standard_Lattice_3d_Facets.Facet_3d_List );

  -- DESCRIPTION :
  --   Returns the list of facets of the ridge of A supported by v.

  -- ON ENTRY :
  --   A        matrix of 4 rows;
  --   v        normal to a facet of convex hull spanned by A;
  --   s        indices to points of A supported by v,
  --            i.e.: s = Standard_Lattice_Supports.Support(A,v).

  -- ON RETURN :
  --   p        pivot of v, p = Standard_Power_Transformations.Pivot(v);
  --   M        matrix to turn facet perpendicular to coordinate axes;
  --   F        list of 3d faces of convex hull of ridge.

  function Filter_non_Vertices
              ( f : Facet_in_4d; v : Standard_Integer_Vectors.Vector )
              return Facet_in_4d;

  -- DESCRIPTION :
  --   Given in v the vertex points for the facet f, if f.points /= v,
  --   then the facet on return has its points correctly adjusted.

  function Ridges_of_Facet
              ( A : Matrix; v : Standard_Integer64_Vectors.Vector )
              return Facet_in_4d;

  -- DESCRIPTION :
  --   Computes the ridges of a facet of the polytope spanned by the points
  --   with coordinates in the columns of A, supported in the direction v.

  function Extreme_Angle
              ( A : Matrix; f : Standard_Integer_Vectors.Vector;
                v,w : Standard_Integer64_Vectors.Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns the index to a column in A, not in f, which makes the
  --   most extreme angle to the facet supported by v pivoting along
  --   the ridge normal w.

  function Inner_Normal
              ( A : Matrix; p : Standard_Integer_Vectors.Vector;
                k : integer32 ) return Standard_Integer64_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the inner normal to the points of A,
  --   with column indices in p and k.

  function Is_Same ( a,b : Standard_Integer_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   The index sets a and b are the same if they have the same length
  --   and every element of a occurs as well in b.

  procedure Neighbors ( A : in Matrix; f : in Link_to_4d_Facet;
                        idcnt : in out integer32 );

  -- DESCRIPTION :
  --   Performs giftwrapping for every ridge of the facet f
  --   to compute facets of conv(A) that are neighbors to f.
  --   The current highest identification number in idcnt gets updated.

  function Convex_Hull_4D 
               ( A : Matrix; v : Standard_Integer64_Vectors.Vector )
               return Facet_4d_List;

  -- DESCRIPTION :
  --   Returns a list of facets supporting the convex hull of A.
  --   The vector v is the first facet normal.

-- SELECTORS :

  function Is_Connected ( f : Link_to_4d_Facet ) return boolean;

  -- DESCRIPTION :
  --   Returns true if all neighbors of f are defined.

  function Pop ( f : Facet_4d_List ) return Link_to_4d_Facet;

  -- DESCRIPTION :
  --   Returns the first facet in the list for which Is_Connected()
  --   returns false.  If all facets in f are connected, then the
  --   facets in the list define the convex hull of a point set.

-- for plain check of Euler characteristic :

  function Vertices ( n : integer32; f : Facet_4d_List )
                    return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   The vector on return contains labels to the vertices of a polytope
  --   spanned by n points with facets in f.

  function Belongs_To ( R : Lists_of_Integer_Vectors.List;
                        p : Standard_Integer_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if p is the same as one of the vectors in R.

  function Edges ( f : Facet_4d_List ) return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION :
  --   The list on return contains indices of points that span edges of
  --   the polytope with facets in f.

  function Ridges ( f : Facet_4d_List ) return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION :
  --   The list on return contains indices of points that span ridges of
  --   the polytope with facets in f.

-- DESTRUCTORS :

  procedure Clear ( f : in out Facet_in_4d );
  procedure Clear ( f : in out Link_to_4d_Facet );
  procedure Clear ( f : in out Array_of_4d_Facets );

  -- DESCRIPTION :
  --   Deallocation of all memory occupied by the facets.

end Standard_Lattice_4d_Facets;
