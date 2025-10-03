with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_VecVecs;
with Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;

package Test_Multitasked_Pipelining is

-- DESCRIPTION :
--   Development of multitasked pipelined processing of the cell indices
--   produced by DEMiCs.

  procedure Write_Cell_Indices;

  -- DESCRIPTION :
  --   Writes the cell indices stored in DEMiCs_Output_Data.

  procedure Write_DEMiCs_Output
              ( dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stable : in boolean; stlb : in double_float;
                verbose : in boolean := true );

  -- DESCRIPTIN :
  --   Writes the cells computed by DEMiCs.

  -- ON ENTRY :
  --   dim      dimension of the points in each support;
  --   mix      type of mixture;
  --   sup      supports of a polynomial system;
  --   verbose  flag for more information.

  procedure Compute_Mixed_Volume ( p : in Laur_Sys );

  -- DESCRIPTION :
  --   Computes the mixed volume of the Newton polytopes
  --   spanned by the supports of p.

  procedure Test_Pipeline
              ( dim,nbtasks : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Computes the lifted supports and runs the pipeline as a test.

  procedure Random_Coefficient_System
              ( dim,nbtasks : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Generates random lifting values, prompts the user for a file name,
  --   and then runs a 2-stage pipeline with number of tasks nbtasks
  --   to solve a random coefficient system with polyhedral homotopies.

  procedure Construct_Mixed_Cells
              ( p : in Laur_Sys; randstart : in boolean;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Constructs the mixed cells in a regular subdivision of 
  --   the Newton polytopes spanned by the supports of p.
  --   If randstart, then a random coefficient system will be
  --   constructed and solved by pipelined polyhedral homotopies.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a polynomial system
  --   and then prepares the input for DEMiCs.

end Test_Multitasked_Pipelining;
