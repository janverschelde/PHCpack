with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer64_Vectors;
with Standard_Integer64_Matrices;        use Standard_Integer64_Matrices;
with Standard_Lattice_3d_Facets;         use Standard_Lattice_3d_Facets;

package Standard_Lattice_3d_Facets_io is

-- DESCRIPTION :
--   Provides basic output routines for facets.

  procedure Write_Coordinates ( A : in Matrix; k : in integer32 );

  -- DESCRIPTION :
  --   Writes the coordinates of the point stored in the k-th column of A.

  procedure Write_Supported_Points
               ( A : in Matrix; v : in Standard_Integer64_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the vector v and the indices to the points in A 
  --   on the supporting hyperplane defined by v.

  procedure Write_Facet ( A : in Matrix; f : in Facet_in_3d );

  -- DESCRIPTION :
  --   Writes a facet of conv(A) to standard output.

  procedure Write_Facets ( A : in Matrix; f : in Facet_3d_List );

  -- DESCRIPTION :
  --   Writes the list of facets of conv(A) to standard output.

end Standard_Lattice_3d_Facets_io;
