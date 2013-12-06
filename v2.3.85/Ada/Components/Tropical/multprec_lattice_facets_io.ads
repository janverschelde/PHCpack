--with Multprec_Integer64_Matrices;        use Multprec_Integer64_Matrices;
with Multprec_Integer_Matrices;          use Multprec_Integer_Matrices;
with Multprec_Lattice_Facets;            use Multprec_Lattice_Facets;

package Multprec_Lattice_Facets_io is

-- DESCRIPTION :
--   Provides basic output routines for facets.

  procedure Write_Coordinates ( A : in Matrix; k : in integer );

  -- DESCRIPTION :
  --   Writes the coordinates of the point stored in the k-th column of A.

  procedure Write_Facet ( A : in Matrix; f : in Facet );

  -- DESCRIPTION :
  --   Writes a facet of conv(A) to multprec output.

  procedure Write_Facets ( A : in Matrix; f : in Facet_List );

  -- DESCRIPTION :
  --   Writes the list of facets of conv(A) to multprec output.

end Multprec_Lattice_Facets_io;
