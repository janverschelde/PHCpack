with Standard_Integer_Vectors;
with Multprec_Integer_Vectors;
with Multprec_Integer_Matrices;
with Multprec_Lattice_4d_Facets;

package Multprec_Lattice_4d_Facets_io is

-- DESCRIPTION :
--   Provides basic output and check for visual inspection of the output
--   of 4 dimensional convex hull computations.

  procedure Check_Inner_Products
              ( A : in Multprec_Integer_Matrices.Matrix;
                v : in Multprec_Integer_Vectors.Vector;
                s : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Checks whether all inner products of the columns in A with v
  --   are all equal for indices in s.

  procedure Write_4D_Facet
              ( A : in Multprec_Integer_Matrices.Matrix;
                f : in Multprec_Lattice_4D_Facets.Link_to_4d_Facet );
 
  -- DESCRIPTION :
  --   Writes information about a 4d facet of the convex hull of A. 

  procedure Write_4D_Facets
              ( A : in Multprec_Integer_Matrices.Matrix;
                f : in Multprec_Lattice_4D_Facets.Facet_4D_List );

  -- DESCRIPTION :
  --   Writes all facets of the 4D convex hull of A.

end Multprec_Lattice_4d_Facets_io;
