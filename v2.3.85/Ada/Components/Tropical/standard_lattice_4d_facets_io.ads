with Standard_Integer_Vectors;
with Standard_Integer64_Vectors;
with Standard_Integer64_Matrices;
with Standard_Lattice_4d_Facets;

package Standard_Lattice_4d_Facets_io is

-- DESCRIPTION :
--   Provides basic output and check for visual inspection of the output
--   of 4 dimensional convex hull computations.

  procedure Check_Inner_Products
              ( A : in Standard_Integer64_Matrices.Matrix;
                v : in Standard_Integer64_Vectors.Vector;
                s : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Checks whether all inner products of the columns in A with v
  --   are all equal for indices in s.

  procedure Write_4D_Facet
              ( A : in Standard_Integer64_Matrices.Matrix;
                f : in Standard_Lattice_4D_Facets.Facet_in_4d );
  procedure Write_4D_Facet
              ( A : in Standard_Integer64_Matrices.Matrix;
                f : in Standard_Lattice_4D_Facets.Link_to_4d_Facet );
 
  -- DESCRIPTION :
  --   Writes information about a 4d facet of the convex hull of A. 

  procedure Write_4D_Facets
              ( A : in Standard_Integer64_Matrices.Matrix;
                f : in Standard_Lattice_4D_Facets.Facet_4D_List );

  -- DESCRIPTION :
  --   Writes all facets of the 4D convex hull of A.

end Standard_Lattice_4d_Facets_io;
