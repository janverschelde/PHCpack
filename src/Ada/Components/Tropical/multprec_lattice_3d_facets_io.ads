with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Integer_Vectors;
with Multprec_Integer_Matrices;          use Multprec_Integer_Matrices;
with Multprec_Lattice_3d_Facets;         use Multprec_Lattice_3d_Facets;

package Multprec_Lattice_3d_Facets_io is

-- DESCRIPTION :
--   Provides basic output routines for facets.

  procedure Write_Coordinates ( A : in Matrix; k : in integer32 );

  -- DESCRIPTION :
  --   Writes the coordinates of the point stored in the k-th column of A.

  procedure Write_Supported_Points
               ( A : in Matrix; v : in Multprec_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the vector v and the indices to the points in A
  --   on the supporting hyperplane defined by v.

  procedure Write_Facet ( A : in Matrix; f : in Facet_in_3d );

  -- DESCRIPTION :
  --   Writes a facet of conv(A) to multprec output.

  procedure Write_Facets ( A : in Matrix; f : in Facet_3d_List );

  -- DESCRIPTION :
  --   Writes the list of facets of conv(A) to multprec output.

end Multprec_Lattice_3d_Facets_io;
