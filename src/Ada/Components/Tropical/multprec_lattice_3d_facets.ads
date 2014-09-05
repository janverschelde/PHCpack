with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Vectors;   
with Multprec_Integer_Matrices;          use Multprec_Integer_Matrices;
with Lists_of_Integer_Vectors;
with Lists_of_Integer64_Vectors;
with Generic_Lists;

package Multprec_Lattice_3d_Facets is

-- DESCRIPTION :
--   A lattice polytope is spanned by points with integer coordinates.
--   The operations in this package compute the entire face lattice for
--   a three dimensional lattice polytope, using multprecision integers.

-- DATA STRUCTURES :

  type Facet_in_3d;
  type Link_to_3d_Facet is access Facet_in_3d;
  type Array_of_3d_Facets is array ( integer32 range <> ) of Link_to_3d_Facet;

  type Facet_in_3d ( d,n : integer32 ) is record
    label : integer32;                        -- identification number
    normal : Multprec_Integer_Vectors.Vector(1..d); 
      -- inner normal to the facet is a vector of multiprecision integers
    points : Standard_Integer_Vectors.Vector(1..n);
      -- indices to vertices of facet is a vector of standard integers
    neighbors : Array_of_3d_Facets(1..n);   -- neighbors share an edge
  end record;

  package Lists_of_3d_Facets is new Generic_Lists(Link_to_3d_Facet);
  type Facet_3d_List is new Lists_of_3d_Facets.List;

  type Array_of_Facet_3d_Lists is
    array ( integer32 range <> ) of Facet_3d_List;

  type Boolean_Array is array ( integer32 range <> ) of boolean;
  type Boolean_Matrix is 
     array ( integer32 range <> , integer32 range <> ) of boolean;

-- INITIAL FACET : compute the first facet

  function Lower ( A : Matrix; i,j : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the i-th column of A has coordinates that are 
  --   lexicographically lower than those in the j-th column.

  function Lowest ( A : Matrix ) return integer32;

  -- DESCRIPTION :
  --   Returns the index of the column in A that has lexicographically 
  --   the lowest coordinates.

  function Second_Lowest ( A : Matrix; k : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the index of the column different from k with the
  --   smallest first coordinate.

  function Largest_Angle ( A : Matrix; k : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the index j of the column in A, j /= k, for which the
  --   vector A(:,j) - A(:,k) makes the largest angle with (1,0,0).

  function Initial_Edge ( A : Matrix; k : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the index of the point in A which will jointly with k
  --   span the initial edge of the polytope spanned by the columns of A.

  function Edge_Normal ( A : Matrix; i,j : integer32 )
                       return Multprec_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns an inner normal for points that span the initial edge,
  --   with coordinates in i-th and j-th column of A.

  function Independent ( A : Matrix; i,j,k : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the columns i, j and k contain points that are
  --   not collinear and are linearly independent.

  function Third_Point
              ( A : Matrix; i,j : integer32; m : Integer_Number;
                v : Multprec_Integer_Vectors.Vector ) return integer32;

  -- DESCRIPTION :
  --   Checks whether there is a facet with inner normal v.
 
  -- ON ENTRY :
  --   A        coordinates for points in space;
  --   i        index to lexicographically lowest point of A;
  --   j        index to second point in A that spans an edge;
  --   m        equals Inner_Product(A,i,v) = Inner_Product(A,j,v);
  --   v        inner normal of edge spanned by points i and j.

  -- ON RETURN :
  --   0 or an index k of a column in A, k /= i, k /= j for which
  --   Inner_Product(A,k,v) = m and Independent(A,i,j,k).

  procedure Normalize ( v : in out Multprec_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Divides the components of the vector v by the greatest common divisor.

  function Normal ( A : Matrix; i,j : integer32;
                    v : Multprec_Integer_Vectors.Vector )
                  return Multprec_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the normal vector perpendicular to the vector along the
  --   edge spanned by points with coordinates in columns i and j of A,
  --   and perpendicular to the vector v.

  -- REQUIRED : v(v'last) = 0.

  function Shift ( A : Matrix; i,j : integer32 )
                 return Multprec_Integer_Vectors.Vector;

  -- DESCRIPTION : returns A(:,j) - A(:,i).

  function Largest_Angle
              ( A : Matrix; i,j : integer32;
                v,w : Multprec_Integer_Vectors.Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns the index to the point with coordinates in A which makes
  --   the largest angle with the vector v.

  -- ON ENTRY :
  --   A        coordinates for a set of points in 3-space;
  --   i        index to lexicographically lowest point of A;
  --   j        index to the second point of the initial edge;
  --   v        normal to the initial edge, spanned by columns i and j;
  --   w        vector perpendicular to v and the initial edge.

  function Extreme_Angle
              ( A : Matrix; i,j : integer32;
                f : Standard_Integer_Vectors.Vector;
                v,w : Multprec_Integer_Vectors.Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns the index to the point with coordinates in A which makes
  --   the largest angle in absolute value with the vector v.

  -- REQUIRED : v is perpendicular to w and moreover, v is inner normal
  --   relative to the other points that span the facet.

  -- ON ENTRY :
  --   A        coordinates for a set of points in 3-space;
  --   i        index to the first vertex points of an edge; 
  --   j        index to end point of an edge;
  --   f        index of points that span the facet;
  --   v        normal to an facet edge, with a vertex in column i;
  --   w        normal to the face spanned by the points in f.

  function Normal ( u,v : Multprec_Integer_Vectors.Vector )
                  return Multprec_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector perpendicular to u and v.

  function Normal ( A : Matrix; i,j,k : integer32 )
                  return Multprec_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector perpendicular to the points with coordinates
  --   in columns i, j, and k.

  function Inner_Normal ( A : Matrix; i,j,k : integer32 )
                        return Multprec_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the inner normal to the facet spanned by the points with
  --   coordinates in columns i, j, and k of A.

  procedure Initial_Facet_Normal
              ( A : in Matrix; i,j,k : out integer32;
                v : out Multprec_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Returns three points of A which span the initial facet of A,
  --   with coordinaetes in columns i, j, and k of A, where v is the
  --   inner normal of the facet.

  function Drop ( A : Matrix; p : integer32 ) return Matrix;

  -- DESCRIPTION :
  --   The matrix on return has the same columns as A,
  --   but with the p-th coordinate dropped.

  function Match ( A,B : Matrix; p,i,j : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if column i of A matches column j of B,
  --   where the p-th coordinate was dropped from B.

  function Match_Vertices
              ( A,V : Matrix; p : integer32 )
              return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   The vector on return has range V'range(2) and contains the
  --   column indices of A for the corresponding vertex.

  -- REQUIRED :
  --   A is the transformed matrix of 3d points
  --   and p is the dropped coordinate from V.

  function Edges_of_Facet
              ( A : Matrix; v : Multprec_Integer_Vectors.Vector )
              return Facet_in_3d;

  -- DESCRIPTION :
  --   Computes the edges of the facet with inner normal v.

  function Initial_Facet ( A : Matrix ) return Facet_in_3d;

  -- DESCRIPTION :
  --   Computes the initial facet of the 3d lattice polygon
  --   spanned by the points with coordinates in A.

-- ROTATE facet normal :

  procedure Invert ( p : in out Standard_Integer_Vectors.Vector;
                     k : in integer32 );

  -- DESCRIPTION :
  --   Inverts the contents of p, reading p backwards starting at k.

  procedure Connect ( f : in Link_to_3d_Facet; p,i,j : in integer32 );

  -- DESCRIPTION :
  --   The p-th neighbor of f shares the edge starting at i and ending 
  --   at the vertex j.  Before connecting the p-th neighbor to f,
  --   the orientation of the vertices may be reversed.

  -- CONNECTIVITY RULE :
  --   Facet f has neighbor g at p, if the edge starting at i = f.points(p)
  --   and ending at j occurs in g as an edge starting at j and ending at i.
  --   Likewise, the q-th neighbor of g will be f, where j = g.points(q).

  procedure Previous_and_Next_Edge
              ( f : in Link_to_3d_Facet; i : in integer32;
                p,q : out integer32 );

  -- DESCRIPTION :
  --   Sets p and q to the index of previous and next edge of f,
  --   both sharing i respectively as end and begin point of an edge.

  procedure Connect ( f,g : in Link_to_3d_Facet; i,j : in integer32;
                      connected : out boolean );

  -- DESCRIPTION :
  --   For facets f and g sharing a common vertex respectively
  --   as their i-th and j-th vertex, connected on return is true
  --   if they share an edge.  If connected is true on return,
  --   then they are connected as neighbors.

  procedure Connect ( f,g : in Link_to_3d_Facet );

  -- DESCRIPTION :
  --   Checks whether the facets f and g share an edge.
  --   If f and g share an edge, then they are connected.

  function Is_Connected ( f : Link_to_3d_Facet ) return boolean;

  -- DESCRIPTION :
  --   Returns true if all neighbors of f are defined,
  --   returns false otherwise.

  procedure Neighbors ( A : in Matrix; f : in Link_to_3d_Facet;
                        idcnt : in out integer32 );

  -- DESCRIPTION :
  --   Computes the neighbors to the facet f of the polytope spanned by A.

  -- ON ENTRY :
  --   A        coordinates of points spanning the polytope;
  --   f        a facet of the polytope;
  --   idcnt    current highest identification number.

  -- ON RETURN :
  --   f        new neighbors of f have label > idcnt on entry;
  --   idcnt    updated highest identification number.

-- MAIN ALGORITHM :

  function Pop ( f : Facet_3d_List ) return Link_to_3d_Facet;
 
  -- DESCRIPTION :
  --   Returns the link to the first facet in f that is not connected.

  procedure Connect ( f : in Facet_3d_List; lft : in Link_to_3d_Facet );

  -- DESCRIPTION :
  --   Connects the facet lft to the facets in the list f.

  function Convex_Hull_3D ( A : Matrix ) return Facet_3d_List;

  -- DESCRIPTION :
  --   Returns a list of facets supporting the convex hull of A.

-- SELECTORS :

  function Support_Value_of_Facet
             ( A : Matrix; f : Facet_in_3d ) return Integer_Number;

  -- DESCRIPTION :
  --   Returns the value of the supporting function defined by the
  --   facet normal and the points spanned by the facet.

  function Facet_Normals
             ( f : Facet_3d_List ) return Lists_of_Integer64_Vectors.List;

  -- DESCRIPTION :
  --   Returns the list of inner normals to the facets in f.

  function Is_Facet_Normal
              ( f : Facet_3d_List; v : Multprec_Integer_Vectors.Vector )
              return boolean;

  -- DESCRIPTION :
  --   Returns true if the vector v occurs among the facet normals in f.

  function Select_Facet_Normals
              ( f : Facet_3d_List; v : Lists_of_Integer64_Vectors.List )
              return Lists_of_Integer64_Vectors.List;

  -- DESCRIPTION :
  --   Returns the list of vectors in v that occur as facet normals in f.

  function Get_Facet
              ( f : Facet_3d_List; v : Multprec_Integer_Vectors.Vector )
              return Link_to_3d_Facet;

  -- DESCRIPTION :
  --   Returns the facet in the list f with inner normal equal to v.

  function Edges ( f : Facet_3d_List ) return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION :
  --   Returns a list of all edges in the list of facets f.

  function Vertices ( n : integer32; f : Facet_3d_List )
                    return Standard_Integer_Vectors.Vector;
  function Vertices ( n : integer32; f : Array_of_3d_Facets )
                    return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector of indices to vertex points of facets in f.

  -- REQUIRED : n = total number of points in the support set.

  function First_Incident_Vertex ( f : Facet_3d_List ) return integer32;

  -- DESCRIPTION :
  --   Returns the first vertex incident to the first facet of f,
  --   as a default vertex to start the walk in the enumerators below.

  procedure Check_Euler_Characteristic
               ( m : in integer32; f : in Facet_3d_List );

  -- DESCRIPTION :
  --   The Euler characteristic of a 3D polytope equals
  --   #facets - #edges + #vertices = 2.

  -- ON ENTRY :
  --   m        total number of points in the support;
  --   f        list of facets spanned by the polytope.

-- ENUMERATORS : walk to enumerate vertices, crawl for edges

  generic

    with procedure Report ( w : in integer32; continue : out boolean );

    -- DESCRIPTION :
    --   For every vertex w connected by an edge path to the original v,
    --   Report is called.  No neighbors to w are visited is continue
    --   is set to false by Report.  Every vertex is visited only once.

  procedure Walk ( f : in Facet_3d_List; v : in integer32;
                   b : in out Boolean_Array );

  -- DESCRIPTION :
  --   Starting at the vertex v, walks to neighboring vertices,
  --   calling Report at every newly visited vertex.

  -- REQUIRED : b'range = A'range(2) and v in A'range(2),
  --   where A is the support set of the polytope with facets in f.

  -- ON ENTRY :
  --   f        list of facets of a polytope;
  --   v        index to a column in v;
  --   b        if b(i) then the i-th column of A is a visited vertex,
  --            otherwise i is not yet visited or not a vertex;
  --            initially all entries in b should be false.

  -- ON RETURN :
  --   b        if b(i), then the i-th point is a vertex, else not.

  generic

    with procedure Report ( v,w : in integer32; continue : out boolean );

    -- DESCRIPTION :
    --   The pair (v,w) spans an edge connected by a path from the original
    --   vertex v given on entry of the Crawl procedure.  Report is called
    --   with every new edge.  No neighbors to w are visited is continue
    --   is set to false by Report.  Every edge is visited only once.

  procedure Crawl ( f : in Facet_3d_List; v : in integer32;
                    b : in out Boolean_Matrix );

  -- DESCRIPTION :
  --   Starting at the vertex v, walks to neighboring vertices,
  --   calling Report at every newly visited edge.

  -- REQUIRED : b'range = A'range(2) and v in A'range(2),
  --   where A is the support set of the polytope with facets in f.

  -- ON ENTRY :
  --   f        list of facets of a polytope;
  --   v        index to a column in v;
  --   b        if b(i,j) then vertices i and j span a visited edge,
  --            otherwise the edge spanned by i and j is not yet visited
  --            or there is no edge spanned by i and j;
  --            initially all entries in b should be false.

  -- ON RETURN :
  --   b        if b(i,j), then there is an edge spanned by i and j,
  --            else there is no edge spanned by i and j.

-- DESTRUCTORS :

  procedure Clear ( f : in out Link_to_3d_Facet );
  procedure Clear ( f : in out Array_of_3d_Facets );
  procedure Clear ( f : in out Facet_3d_List );

  -- DESCRIPTION :
  --   Deallocation of all memory occupied by the facets.

end Multprec_Lattice_3d_Facets;
