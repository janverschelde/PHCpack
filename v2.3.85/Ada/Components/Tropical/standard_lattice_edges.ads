with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer64_Matrices;        use Standard_Integer64_Matrices;
with Standard_Lattice_3d_Facets;         use Standard_Lattice_3d_Facets;
with Generic_Lists;

package Standard_Lattice_Edges is

-- DESCRIPTION :
--   Every edge of a three dimensional polytope is the intersection
--   of exactly two facets.

-- DATA STRUCTURES :

  type Edge;
  type Link_to_Edge is access Edge;

  type Edge is record       -- an edge has a label,
    label,a,b : integer32;  -- is spanned by vertices a and b, and is
    f,g : Link_to_3d_Facet; -- in the intersection of two facets f and g
  end record;

  package Lists_of_Edges is new Generic_Lists(Link_to_Edge);
  type Edge_List is new Lists_of_Edges.List;

-- CONSTRUCTORS :

  function Edges_of_3D_Hull ( A : Matrix ) return Edge_List;

  -- DESCRIPTION :
  --   Returns the list of edges of the convex hull spanned by the
  --   points with coordinates in the columns of A.

  function Edges_of_3D_Hull
              ( m : integer32; f : Facet_3d_List ) return Edge_List;

  -- DESCRIPTION :
  --   Returns the list of edges of the polytope defined by the facets in f,
  --   of a polytope spanned by m points.

-- DESTRUCTORS :

  procedure Clear ( e : in out Link_to_Edge );
  procedure Clear ( e : in out Edge_List );

  -- DESCRIPTION :
  --   Deallocation of all memory occupied by the edges.

end Standard_Lattice_Edges;
