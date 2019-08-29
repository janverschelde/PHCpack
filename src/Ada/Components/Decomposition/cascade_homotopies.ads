with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
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

package Cascade_Homotopies is

-- DESCRIPTION :
--   A cascade homotopy defines a sequence of homotopies to compute
--   to compute generic points on all components of the solution set.
--   Versions of the cascade homotopies are provided
--   1) for ordinary and for Laurent polynomial systems, and
--   2) in double, double double, and quad double precision.
--   The six versions are combined with various levels of output:
--   1) the output is written to one single file, or
--   2) each superwitness set is written to a separate file.

  procedure Witness_Generate
               ( outfile : in file_type; nt : in natural32;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );
  procedure Witness_Generate
               ( outfile : in file_type; nt : in natural32;
                 ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );
  procedure Witness_Generate
               ( outfile : in file_type; nt : in natural32;
                 ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );
  procedure Witness_Generate
               ( outfile : in file_type; nt : in natural32;
                 ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );
  procedure Witness_Generate
               ( outfile : in file_type; nt : in natural32;
                 ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );
  procedure Witness_Generate
               ( outfile : in file_type; nt : in natural32;
                 ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );

  -- DESCRIPTION :
  --   Calculates candidate witness points on every component,
  --   starting at the component of dimension k,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   outfile   file for intermediate results and diagnostics;
  --   nt        number of tasks, set to zero for no tasking;
  --   ep        embedded polynomial system;
  --   sols      solutions to the system ep (unfiltered);
  --   topdim    number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   lowdim    lower bound on the dimension to stop the cascade;
  --   zerotol   tolerance to decide whether a number is zero or not;
  --   tolsing   tolerance on rco to decide whether a solution is singular.

  procedure Witness_Generate
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );
  procedure Witness_Generate
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );
  procedure Witness_Generate
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );
  procedure Witness_Generate
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );
  procedure Witness_Generate
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );
  procedure Witness_Generate
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );

  -- DESCRIPTION :
  --   Calculates candidate witness points on every component,
  --   starting at the component of dimension k,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   outfile   file for intermediate results and diagnostics;
  --   resfile   file for final results:
  --               1) system at each level, for k down to 0; and
  --               2) solutions with zero slack variables;
  --   nt        number of tasks, set to zero for no tasking;
  --   ep        embedded polynomial system;
  --   sols      solutions to the system ep (unfiltered);
  --   topdim    number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   lowdim    lower bound on the dimension to stop the cascade;
  --   zerotol   tolerance to decide whether a number is zero or not;
  --   tolsing   tolerance on rco to decide whether a solution is singular.

  procedure Witness_Generate
              ( name : in string; outfile : in file_type;
                nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );

  -- DESCRIPTION :
  --   This witness generate writes the witness supersets to files,
  --   and returns the superwitness sets,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   name      file name for the top embedded system;
  --   outfile   file for all intermediate and final results;  
  --   nt        number of tasks for multitasking, set to zero for no tasking;
  --   ep        embedded polynomial system;
  --   sols      solutions to the system ep (unfiltered);
  --   topdim    number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   lowdim    lower bound on the dimension to stop the cascade;
  --   zerotol   tolerance to decide whether a number is zero or not;
  --   tolsing   tolerance on rco to decide whether a solution is singular.

  -- ON RETURN :
  --   embsys    sequence of embedded polynomial systems;
  --   esols0    candidate witness points at each dimension;
  --   pathcnts  table with path counts during the cascade homotopies;
  --   times     CPU time at each stage in the cascade homotopy;
  --   alltime   the total elapsed CPU time.

  procedure Witness_Generate
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );
  procedure Witness_Generate
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );
  procedure Witness_Generate
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );
  procedure Witness_Generate
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );
  procedure Witness_Generate
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );
  procedure Witness_Generate
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32;
                 zerotol,tolsing : in double_float );

  -- DESCRIPTION :
  --   This witness generate writes the witness supersets to files,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   name      file name for the top embedded system;
  --   outfile   file for all intermediate and final results;  
  --   nt        number of tasks for multitasking, set to zero for no tasking;
  --   ep        embedded polynomial system;
  --   sols      solutions to the system ep (unfiltered);
  --   topdim    number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   lowdim    lower bound on the dimension to stop the cascade;
  --   zerotol   tolerance to decide whether a number is zero or not;
  --   tolsing   tolerance on rco to decide whether a solution is singular.

-- SILENT VERSIONS :

  procedure Witness_Generate
              ( nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Witness_Generate
              ( nt : in natural32;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Witness_Generate
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Witness_Generate
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Witness_Generate
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Witness_Generate
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );

  -- DESCRIPTION :
  --   Silent version of the witness generate to compute witness supersets,
  --   in double, double double, and quad double precision.

  -- ON ENTRY :
  --   nt        number of tasks for multitasking, set to zero for no tasking;
  --   ep        embedded polynomial system;
  --   sols      solutions to the system ep (unfiltered);
  --   topdim    number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   lowdim    lower bound on the dimension to stop the cascade;
  --   zerotol   tolerance to decide whether a number is zero or not;
  --   tolsing   tolerance on rco to decide whether a solution is singular.

  -- ON RETURN :
  --   embsys    sequence of embedded polynomial systems;
  --   esols0    candidate witness points at each dimension;
  --   pathcnts  table with path counts during the cascade homotopies;
  --   times     CPU time at each stage in the cascade homotopy;
  --   alltime   the total elapsed CPU time.

-- OUTPUT WITH CALLBACK PROCEDURE :

  procedure Witness_Generate_Callback
              ( nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Witness_Generate_Callback
              ( nt : in natural32;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Witness_Generate_Callback
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Witness_Generate_Callback
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Witness_Generate_Callback
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Witness_Generate_Callback
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );

  -- DESCRIPTION :
  --   Silent version of the witness generate to compute witness supersets,
  --   in double, double double, and quad double precision,
  --   with a callback procedure, which is called on each new witness set.

  -- ON ENTRY :
  --   nt        number of tasks for multitasking, set to zero for no tasking;
  --   ep        embedded polynomial system;
  --   sols      solutions to the system ep (unfiltered);
  --   topdim    number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   lowdim    lower bound on the dimension to stop the cascade;
  --   zerotol   tolerance to decide whether a number is zero or not;
  --   tolsing   tolerance on rco to decide whether a solution is singular.

  -- ON RETURN :
  --   embsys    sequence of embedded polynomial systems;
  --   esols0    candidate witness points at each dimension;
  --   pathcnts  table with path counts during the cascade homotopies;
  --   times     CPU time at each stage in the cascade homotopy;
  --   alltime   the total elapsed CPU time.

  procedure Witness_Generate
              ( file : in file_type; nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Witness_Generate
              ( file : in file_type; nt : in natural32;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Witness_Generate
              ( file : in file_type; nt : in natural32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Witness_Generate
              ( file : in file_type; nt : in natural32;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Witness_Generate
              ( file : in file_type; nt : in natural32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );
  procedure Witness_Generate
              ( file : in file_type; nt : in natural32;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration );

  -- DESCRIPTION :
  --   Reporting version of the witness generate to compute witness supersets,
  --   in double, double double, and quad double precision.

  -- ON ENTRY :
  --   file      must be opened for output;
  --   nt        number of tasks for multitasking, set to zero for no tasking;
  --   ep        embedded polynomial system;
  --   sols      solutions to the system ep (unfiltered);
  --   topdim    number of slack variables and random hyperplanes,
  --             equals the top dimension of the solution sets;
  --   lowdim    lower bound on the dimension to stop the cascade;
  --   zerotol   tolerance to decide whether a number is zero or not;
  --   tolsing   tolerance on rco to decide whether a solution is singular.

  -- ON RETURN :
  --   embsys    sequence of embedded polynomial systems;
  --   esols0    candidate witness points at each dimension;
  --   pathcnts  table with path counts during the cascade homotopies;
  --   times     CPU time at each stage in the cascade homotopy;
  --   alltime   the total elapsed CPU time.

end Cascade_Homotopies;
