with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;

package Embeddings_and_Cascades is

-- DESCRIPTION :
--   The procedures in this package embed a given system,
--   apply the blackbox solver to the top dimensional system
--   and launch the cascades of homotopies to compute all
--   lower dimensional solution sets.

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

  procedure Standard_Embed_and_Cascade
              ( nt : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean );
  procedure Standard_Embed_and_Cascade
              ( nt : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean );
  procedure DoblDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean );
  procedure DoblDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean );
  procedure QuadDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean );
  procedure QuadDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean );

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
  --            otherwise, the output sets may still be reducible.

  procedure Standard_Embed_and_Cascade
              ( file : in file_type; name : in string;
                nt : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean );
  procedure Standard_Embed_and_Cascade
              ( file : in file_type; name : in string;
                nt : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean );
  procedure DoblDobl_Embed_and_Cascade
              ( file : in file_type; name : in string;
                nt : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean );
  procedure DoblDobl_Embed_and_Cascade
              ( file : in file_type; name : in string;
                nt : in natural32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean );
  procedure QuadDobl_Embed_and_Cascade
              ( file : in file_type; name : in string;
                nt : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean );
  procedure QuadDobl_Embed_and_Cascade
              ( file : in file_type; name : in string;
                nt : in natural32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean );

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
  --            otherwise, the output sets may still be reducible.

end Embeddings_and_Cascades;
