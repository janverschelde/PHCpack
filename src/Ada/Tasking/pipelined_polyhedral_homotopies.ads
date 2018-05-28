with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_VecVecs;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;

package Pipelined_Polyhedral_Homotopies is

-- DESCRIPTION :
--   A 2-stage pipeline feeds the mixed cells computed by DEMiCs
--   directly to the multitasked path trackers to solve a random
--   coefficient start system.

  procedure Pipeline_Cells_to_Paths
              ( dim,nt : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lif : in Standard_Floating_VecVecs.Link_to_VecVec;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Computes the lifted supports and runs the pipeline with nt tasks.

  -- ON ENTRY :
  --   dim      dimension of the points before the lifting;
  --   nt       number of tasks, must be at least two;
  --   mix      type of mixture;
  --   sup      supports of the system, orded along occurrence;
  --   lif      random lifting values for the points in the supports;
  --   verbose  if additional output is needed;

  -- ON RETURN :
  --   q        random coefficient start system with same supports as sup;
  --   qsols    solutions of the start system q.

end Pipelined_Polyhedral_Homotopies;
