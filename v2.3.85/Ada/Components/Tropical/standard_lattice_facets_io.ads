with Standard_Integer64_Matrices;
with Standard_Lattice_Facets;

package Standard_Lattice_Facets_io is

-- DESCRIPTION :
--   Provides basic output and check for visual inspection of the output
--   of d dimensional convex hull computations, for d >= 4.

  procedure Write_Facet
              ( A : in Standard_Integer64_Matrices.Matrix;
                f : in Standard_Lattice_Facets.Link_to_Facet );
 
  -- DESCRIPTION :
  --   Writes information about a facet of the convex hull of A. 

  procedure Write_Facets
              ( A : in Standard_Integer64_Matrices.Matrix;
                f : in Standard_Lattice_Facets.Facet_List );

  -- DESCRIPTION :
  --   Writes all facets of the convex hull of A.

end Standard_Lattice_Facets_io;
