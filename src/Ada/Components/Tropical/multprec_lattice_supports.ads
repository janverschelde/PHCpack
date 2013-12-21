with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
--with Multprec_Integer64_Numbers;         use Multprec_Integer64_Numbers;
--with Multprec_Integer64_Vectors;
--with Multprec_Integer64_Matrices;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Vectors;
with Multprec_Integer_Matrices;

package Multprec_Lattice_Supports is

-- DESCRIPTION :
--   This package collects utility functions to deal with supports
--   of lattice polytopes, spanned by points with exact multiprecision
--   integer coordinates.

  function Inner_Product
              ( u,v : Multprec_Integer_Vectors.Vector )
              return Integer_Number;

  -- DESCRIPTION :
  --   Returns the inner product of the vectors u and v.

  function Inner_Product
              ( A : Multprec_Integer_Matrices.Matrix; k : integer32;
                v : Multprec_Integer_Vectors.Vector ) return Integer_Number;

  -- DESCRIPTION :
  --   Returns the inner product of the point with coordinates
  --   in the k-th column of A with the vector v.

  function Inner_Product
              ( A : Multprec_Integer_Matrices.Matrix; i,j : integer32 )
              return Integer_Number;

  -- DESCRIPTION :
  --   Returns the inner product of the columns i and j of A.

  function Inner_Products
              ( A : Multprec_Integer_Matrices.Matrix;
                v : Multprec_Integer_Vectors.Vector )
              return Multprec_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector of all inner products of the points with
  --   coordinates in the columns of A with the vector v.

  function Minimum ( A : Multprec_Integer_Matrices.Matrix;
                     v : Multprec_Integer_Vectors.Vector )
                   return Integer_Number;

  -- DESCRIPTION :
  --   Returns the minimal values of the inner product of the points in
  --   all columns with the vector v.

  function Support ( A : Multprec_Integer_Matrices.Matrix;
                     v : Multprec_Integer_Vectors.Vector )
                   return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector of indices for those points of A supported by v.

  function Support ( A,B : Multprec_Integer_Matrices.Matrix;
                     v : Multprec_Integer_Vectors.Vector )
                   return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the vector of indices to points of A and B supported by v.

  procedure Inner ( A : in Multprec_Integer_Matrices.Matrix;
                    i,j : in integer32;
                    v : in out Multprec_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Changes the orientation of the normal v to the edge spanned by
  --   points with coordinates in columns i and j of A.

  procedure Inner ( A : in Multprec_Integer_Matrices.Matrix;
                    i,j : in integer32;
                    f : in Standard_Integer_Vectors.Vector;
                    v : in out Multprec_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Changes the orientation of the normal v to the edge spanned by
  --   points with coordinates in columns i and j of A, relative to
  --   the other points in the facet with indices in the vector f.

  procedure Inner ( A : in Multprec_Integer_Matrices.Matrix;
                    i,j,k : in integer32;
                    v : in out Multprec_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Changes the orientation of the normal vector v, perpendicular to the
  --   vectors defined by the points in columns i, j, and k of A, so the
  --   vector v is an inner normal to the facet spanned by i, j, and k.

  function Support_Points
              ( A : Multprec_Integer_Matrices.Matrix;
                s : Standard_Integer_Vectors.Vector )
              return Multprec_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a matrix with points selected from the columns of A,
  --   defined by indices in s.

  function Equal ( A,B : Multprec_Integer_Matrices.Matrix; i,j : integer32 )
                 return boolean;

  -- DESCRIPTION :
  --   Returns true if column i of A equals column j of B.

  -- REQUIRED : i in A'range(2), j in B'range(2) and A'range(1) = B'range(1).

  function Indices ( A,V : Multprec_Integer_Matrices.Matrix )
                   return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the column indices for A of the vertex points in V.

end Multprec_Lattice_Supports;
