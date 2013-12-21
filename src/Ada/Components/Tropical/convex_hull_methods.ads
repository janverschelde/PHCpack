with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer64_Vectors;
with Standard_Integer64_Matrices;
with Multprec_Integer_Vectors;
with Multprec_Integer_Matrices;
with Standard_Lattice_Edges;
with Standard_Lattice_3d_Facets;
with Standard_Lattice_4d_Facets;
with Multprec_Lattice_Edges;
with Multprec_Lattice_3d_Facets;
with Multprec_Lattice_4d_Facets;

package Convex_Hull_Methods is

-- DESCRIPTION :
--   This package contains drivers to call several gift wrapping algorithms
--   to construct the convex hull of a set of lattice points, with as well
--   procedures to test the output.

  function Is_In ( A : Standard_Integer64_Matrices.Matrix;
                   k : integer32 ) return boolean;
  function Is_In ( A : Multprec_Integer_Matrices.Matrix;
                   k : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the kth column in A already occurs
  --   in a column j of A, for j < k.

  function Filter_Duplicates
             ( A : Standard_Integer64_Matrices.Matrix )
             return Standard_Integer64_Matrices.Matrix;
  function Filter_Duplicates
             ( A : Multprec_Integer_Matrices.Matrix )
             return Multprec_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   The matrix on return contains no duplicate columns.

  function User_Input ( n,m : integer32 )
                      return Standard_Integer64_Matrices.Matrix;
  function User_Input ( n,m : integer32 )
                      return Multprec_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns an n-by-m integer matrix, entered by the user.

  function Convert ( A : Standard_Integer64_Matrices.Matrix )
                   return Multprec_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Converts a standard integer matrix into multiprecision format.

  function Random_Data ( n,m : integer32 )
                       return Standard_Integer64_Matrices.Matrix;
  function Random_Data ( n,m : integer32 )
                       return Multprec_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Prompts the user for lower and upper bounds on the entries
  --   of a random n-by-m integer matrix.

  function Cyclic_Polytope
              ( d : integer32; t : Standard_Integer64_Vectors.Vector )
              return Standard_Integer64_Matrices.Matrix;
  function Cyclic_Polytope
              ( d : integer32; t : Multprec_Integer_Vectors.Vector )
              return Multprec_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a d-by-n matrix where n = t'length with in the columns
  --   the points that span a cyclic d-polytope.

  function Cyclic_Polytope 
              ( d,n : integer32 )
              return Standard_Integer64_Matrices.Matrix;
  function Cyclic_Polytope 
              ( d,n : integer32 )
              return Multprec_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Prompts the user for n points to interpolate the moment curve,
  --   or allows for random generation of these points, and then
  --   returns the points that span a cyclic polytope in d-space.

  procedure Standard_Planar_Hull
              ( A : in Standard_Integer64_Matrices.Matrix );

  -- DESCRIPTION :
  --   Computes the convex hull for a set of points in the plane,
  --   using standard machine 64-bit arithmetic.

  procedure Multprec_Planar_Hull
              ( A : in Multprec_Integer_Matrices.Matrix );

  -- DESCRIPTION :
  --   Computes the convex hull for a set of points in the plane,
  --   using exact multiprecision integer arithmetic.

  procedure Standard_Show_Initial_Facet
               ( A : in Standard_Integer64_Matrices.Matrix );

  -- DESCRIPTION :
  --   Shows the computation of the initial facet of the three dimensional
  --   polytope spanned by the points in the columns of A.

  procedure Multprec_Show_Initial_Facet
               ( A : in Multprec_Integer_Matrices.Matrix );

  -- DESCRIPTION :
  --   Shows the computation of the initial facet of the three dimensional
  --   polytope spanned by the points in the columns of A.

  procedure Standard_Show_Neighbors
               ( A : in Standard_Integer64_Matrices.Matrix;
                 f : in Standard_Lattice_3d_Facets.Link_to_3d_Facet );

  -- DESCRIPTION :
  --   Computes the neighbors to the facet f of the polytope spanned by A.

  procedure Multprec_Show_Neighbors
               ( A : in Multprec_Integer_Matrices.Matrix;
                 f : in Multprec_Lattice_3d_Facets.Link_to_3d_Facet );

  -- DESCRIPTION :
  --   Computes the neighbors to the facet f of the polytope spanned by A.

  procedure Standard_Start_3D_Giftwrapping
              ( A : in Standard_Integer64_Matrices.Matrix );

  -- DESCRIPTION :
  --   Shows initial stages of the giftwrapping method to compute
  --   the convex hull for a set of points in 3-space,
  --   using standard machine 64-bit arithmetic.

  procedure Multprec_Start_3D_Giftwrapping
              ( A : in Multprec_Integer_Matrices.Matrix );

  -- DESCRIPTION :
  --   Shows initial stages of the giftwrapping method to compute
  --   the convex hull for a set of points in 3-space,
  --   using exact multiprecision integer arithmetic.

  procedure Walk_Vertices
               ( m : in integer32;
                 f : in Standard_Lattice_3d_Facets.Facet_3d_List );

  -- DESCRIPTION :
  --   Applies the Walk procedure to enumerate all vertices
  --   of the polytope spanned by m points and facets in f.

  procedure Walk_Edges
               ( m : in integer32;
                 f : in Standard_Lattice_3d_Facets.Facet_3d_List );

  -- DESCRIPTION :
  --   Applies the Crawl procedure to enumerate all edges
  --   of the polytope spanned by m points and facets in f.

  procedure Check_Edge
              ( e : in Standard_Lattice_Edges.Link_to_Edge;
                bug : out boolean );
  procedure Check_Edge
              ( e : in Multprec_Lattice_Edges.Link_to_Edge;
                bug : out boolean );

  -- DESCRIPTION :
  --   Verifies that the points that span e are indeed in the 
  --   intersecting facets.

  procedure Check_Edges
              ( A : in Standard_Integer64_Matrices.Matrix;
                f : in Standard_Lattice_3d_Facets.Facet_3d_List );
  procedure Check_Edges
              ( A : in Multprec_Integer_Matrices.Matrix;
                f : in Multprec_Lattice_3d_Facets.Facet_3d_List );

  -- DESCRIPTION :
  --   Checks the property that every edge is the intersection
  --   of exactly two facets.  

  procedure Standard_3D_Giftwrapping
                ( A : in Standard_Integer64_Matrices.Matrix );

  -- DESCRIPTION :
  --   Performs giftwrapping on the points stored in the columns of A,
  --   using standard 64-bit integer arithmetic.

  -- REQUIRED :
  --   The matrix A has three rows and
  --   all duplicate columns have been removed from A.

  procedure Multprec_3D_Giftwrapping
                ( A : in Multprec_Integer_Matrices.Matrix );

  -- DESCRIPTION :
  --   Performs giftwrapping on the points stored in the columns of A,
  --   using multiprecision integer arithmetic.

  -- REQUIRED :
  --   The matrix A has three rows and
  --   all duplicate columns have been removed from A.

  procedure Standard_3D_Hull
              ( A : in Standard_Integer64_Matrices.Matrix );

  -- DESCRIPTION :
  --   Computes the convex hull for a set of points in 3D,
  --   using standard machine 64-bit arithmetic.

  procedure Multprec_3D_Hull
              ( A : in Multprec_Integer_Matrices.Matrix );

  -- DESCRIPTION :
  --   Computes the convex hull for a set of points in 3D 
  --   using multiprecision integer arithmetic.

  procedure Standard_Start_4D_hull
              ( A : in Standard_Integer64_Matrices.Matrix;
                v : in Standard_Integer64_Vectors.Vector );

  -- DESCRIPTION :
  --   Starts the construction of the convex hull for a point configuration A
  --   in 4-space, with v given as the first facet normal.

  procedure Multprec_Start_4D_hull
              ( A : in Multprec_Integer_Matrices.Matrix;
                v : in Multprec_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Starts the construction of the convex hull for a point configuration A
  --   in 4-space, with v given as the first facet normal.

  procedure Standard_4D_Euler_Characteristic
              ( n : in integer32;
                f : in Standard_Lattice_4d_Facets.Facet_4d_List );
  procedure Multprec_4D_Euler_Characteristic
              ( n : in integer32;
                f : in Multprec_Lattice_4d_Facets.Facet_4d_List );

  -- DESCRIPTION :
  --   Plain check on the Euler characteristic of a 4-dimensional
  --   polytope spanned by n points and with facets in f.

  procedure Standard_4D_Hull
              ( A : in Standard_Integer64_Matrices.Matrix );

  -- DESCRIPTION :
  --   Computes the convex hull in four dimensions,
  --   using standard 64-bit integer arithmetic.

  -- REQUIRED : A has exactly four rows.

  procedure Multprec_4D_Hull
              ( A : in Multprec_Integer_Matrices.Matrix );

  -- DESCRIPTION :
  --   Computes the convex hull in four dimensions,
  --   using multiprecision integer arithmetic.

  -- REQUIRED : A has exactly four rows.

  procedure Standard_Start_hull
              ( A : in Standard_Integer64_Matrices.Matrix;
                v : in Standard_Integer64_Vectors.Vector );

  -- DESCRIPTION :
  --   Starts the construction of the convex hull for a point configuration A
  --   in d-space, d = A'last(1) with v given as the first facet normal.

  procedure Standard_General_Hull
              ( A : in Standard_Integer64_Matrices.Matrix );

  -- DESCRIPTION :
  --   Computes the convex hull in general dimensions,
  --   using standard 64-bit integer arithmetic.

  procedure Multprec_General_Hull
              ( A : in Multprec_Integer_Matrices.Matrix );

  -- DESCRIPTION :
  --   Computes the convex hull in general dimensions,
  --   using multiprecision integer arithmetic.

end Convex_Hull_Methods;
