with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Embeddings_and_Cascades is

-- DESCRIPTION :
--   The procedures in this package embed a given system,
--   apply the blackbox solver to the top dimensional system
--   and launch the cascades of homotopies to compute all
--   lower dimensional solution sets.
--   The are two types of polynomial systems accepted on input :
--   (1) ordinary polynomial systems, with nonnegative exponents;
--   (2) Laurent polynomial systems, with exponents which may be negative.
--   Three different levels of precision are supported :
--   (1) standard double precision, as provided by the hardware;
--   (2) double double precision, the double of hardware precision;
--   (3) quad double precision quadruples the hardware precision.
--   Parallelism is indicated by the number of tasks :
--   (1) if the number of tasks is zero, then no multitasking is applied;
--   (2) a positive number of tasks defines the level of multitasking.
--   The output comes in three different forms :
--   (1) all output is written to screen;
--   (2) all output is written to file;
--   (3) all output is passed to a callback procedure.
--   So in total, there are 2x3x2x3 = 36 procedures to run the
--   homotopy cascades to compute a numerical irreducible decomposition.

  function Lower_Dimension ( nq,nv : in natural32 ) return natural32;

  -- DESCRIPTION :
  --   Given in nq are the number of polynomials
  --   and in nv are the number of variables of a polynomial system,
  --   the lower bound on the dimension of the solution set
  --   is nv - nq, if nq >= nv, otherwise, the lower bound is 0.

  procedure Prompt_for_Top_Dimension
              ( nq,nv : in natural32; topdim,lowdim : out natural32 );

  -- DESCRIPTION :
  --   Displays the number of equations and variables,
  --   prompts the user to change the default top dimension
  --   and returns the top dimension and the lowest dimension.
  --   Checks whether what gets entered is not larger than nv-1
  --   and that it is at least larger than the lowest dimension
  --   for underdetermined inputs.

  -- ON ENTRY :
  --   nq       number of equations in the system;
  --   nv       number of variables in the system.

  -- ON RETURN :
  --   topdim   top dimension of the solution set;
  --   lowdim   computed lowest dimension of the solution set.

-- ALL OUTPUT WRITTEN TO SCREEN :

  procedure Standard_Embed_and_Cascade
              ( nt : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure Standard_Embed_and_Cascade
              ( nt : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure DoblDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure DoblDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure QuadDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure QuadDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts the user for the top dimension of the solution set
  --   of the system p and then runs the algorithms to compute
  --   generic points on all components of the solution set,
  --   in standard double, double double, or quad double precision.
  --   If requested, homotopy membership tests are applied and
  --   the filtered generic points are classified according to the
  --   irreducible factors in a numerical irreducible decomposition.
  --   All output is written to screen.

  -- ON ENTRY :
  --   nt       number of tasks, zero for no multitasking;
  --   p        a (Laurent) polynomial system;
  --   filter   flag to indicate that the homotopy membership tests
  --            will remove the junk points from the output of cascades,
  --            if false, the output will be superwitness sets;
  --   factor   if filter and factor, then numerical representations for
  --            the irreducible factors will be computed,
  --            otherwise, the output sets may still be reducible;
  --   verbose  is the verbose level.

-- ALL OUTPUT WRITTEN TO FILE :

  procedure Standard_Embed_and_Cascade
              ( file : in file_type; name : in string;
                nt : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure Standard_Embed_and_Cascade
              ( file : in file_type; name : in string;
                nt : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure DoblDobl_Embed_and_Cascade
              ( file : in file_type; name : in string;
                nt : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure DoblDobl_Embed_and_Cascade
              ( file : in file_type; name : in string;
                nt : in natural32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure QuadDobl_Embed_and_Cascade
              ( file : in file_type; name : in string;
                nt : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure QuadDobl_Embed_and_Cascade
              ( file : in file_type; name : in string;
                nt : in natural32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Prompts the user for the top dimension of the solution set
  --   of the system p and then runs the algorithms to compute
  --   generic points on all components of the solution set,
  --   in standard double, double double, or quad double precision.
  --   If requested, homotopy membership tests are applied and
  --   the filtered generic points are classified according to the
  --   irreducible factors in a numerical irreducible decomposition.
  --   Output is written to files.

  -- ON ENTRY :
  --   file     for output, to write diagnostics and results;
  --   name     witness set output will be written to files
  --            which have the string name as prefix;
  --   nt       number of tasks, zero for no multitasking;
  --   p        a (Laurent) polynomial system;
  --   filter   flag to indicate that the homotopy membership tests
  --            will remove the junk points from the output of cascades,
  --            if false, the output will be superwitness sets;
  --   factor   if filter and factor, then numerical representations for
  --            the irreducible factors will be computed,
  --            otherwise, the output sets may still be reducible;
  --   verbose  is the verbose level.

-- ALL OUTPUT TO CALLBACK PROCEDURE :

  procedure Standard_Solve_with_Callback
              ( nt,topdim,lowdim : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Standard_Solve_with_Callback
              ( nt,topdim,lowdim : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure DoblDobl_Solve_with_Callback
              ( nt,topdim,lowdim : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure DoblDobl_Solve_with_Callback
              ( nt,topdim,lowdim : in natural32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure QuadDobl_Solve_with_Callback
              ( nt,topdim,lowdim : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure QuadDobl_Solve_with_Callback
              ( nt,topdim,lowdim : in natural32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );

  -- DESCRIPTION :
  --   Given the top dimension of the solution set of the system p,
  --   computes generic points on all solution components, 
  --   in standard double, double double, or quad double precision.
  --   If requested, homotopy membership tests are applied and
  --   the filtered generic points are classified according to the
  --   irreducible factors in a numerical irreducible decomposition.
  --   No output is written to screen or file.

  -- ON ENTRY :
  --   nt       number of tasks, zero for no multitasking;
  --   topdim   the top dimension of the solution set;
  --   lowdim   lower bound on the dimension to stop the cascade;
  --   p        a (Laurent) polynomial system;
  --   filter   flag to indicate that the homotopy membership tests
  --            will remove the junk points from the output of cascades,
  --            if false, the output will be superwitness sets;
  --   factor   if filter and factor, then numerical representations for
  --            the irreducible factors will be computed,
  --            otherwise, the output sets may still be reducible;
  --   Report_Witness_Set is a callback procedure, called each time
  --            a new witness set is computed.

  -- ON RETURN :
  --   idxfac   indices to the irreducible factors in the decomposition,
  --            if filter and factor are both set to true on input.

end Embeddings_and_Cascades;
