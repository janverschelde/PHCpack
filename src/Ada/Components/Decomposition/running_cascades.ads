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

package Running_Cascades is

-- DESCRIPTION :
--   A cascade is a sequence of homotopies to compute a numerical
--   irreducible decomposition in a top down manner.
--   The procedures in this package take on input the solutions of
--   an embedded system and then compute witness sets.
--   Both ordinary and Laurent polynomial systems are supported.
--   the three different levels of precision are double,
--   double double, and quad double.

  procedure Standard_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure Standard_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure DoblDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure DoblDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure QuadDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure QuadDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given an embedding of the top dimensional solution set,
  --   runs a cascade of homotopies, in standard double, double double,
  --   or quad double precision.  All output is written to screen.

  -- ON ENTRY :
  --   nt       number of tasks for multitasking,
  --            if zero, then no multitasking will be used;
  --   topdim   the top dimension of the solution set;
  --   lowdim   lower bound on the dimension to stop the cascade;
  --   embsys   an embedded system for the top dimension;
  --   sols     solutions of the system embsys;
  --   filter   if true, then junk points will be removed,
  --            otherwise, the output will be superwitness sets.
  --   factor   if true and filter, then the filtered witness sets will be
  --            factored into irreducible components,
  --            otherwise, the output sets may still be reducible;
  --   verbose  is the verbose level.

  procedure Standard_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure Standard_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure DoblDobl_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure DoblDobl_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure QuadDobl_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 );
  procedure QuadDobl_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given an embedding of the top dimensional solution set,
  --   runs a cascade of homotopies, in standard double, double double,
  --   or quad double precision.  All output is written to file.

  -- ON ENTRY :
  --   file     file opened for output;
  --   name     file name for the top embedded system;
  --   nt       number of tasks for multitasking,
  --            if zero, then no multitasking will be used;
  --   topdim   the top dimension of the solution set;
  --   lowdim   lower bound on the dimension to stop the cascade;
  --   embsys   an embedded system for the top dimension;
  --   sols     solutions of the system embsys;
  --   filter   if true, then junk points will be removed,
  --            otherwise, the output will be superwitness sets.
  --   factor   if true and filter, then the filtered witness sets will be
  --            factored into irreducible components,
  --            otherwise, the output sets may still be reducible;
  --   verbose  is the verbose level.

  procedure Standard_Cascade_Callback
              ( nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure Standard_Cascade_Callback
              ( nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure DoblDobl_Cascade_Callback
              ( nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure DoblDobl_Cascade_Callback
              ( nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure QuadDobl_Cascade_Callback
              ( nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );
  procedure QuadDobl_Cascade_Callback
              ( nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) );

  -- DESCRIPTION :
  --   Given an embedding of the top dimensional solution set,
  --   runs a cascade of homotopies, in standard double, double double,
  --   or quad double precision.  No output is written to screen or file,
  --   the callback procedure is call for each new witness set.

  -- ON ENTRY :
  --   nt       number of tasks for multitasking,
  --            if zero, then no multitasking will be used;
  --   topdim   the top dimension of the solution set;
  --   lowdim   lower bound on the dimension to stop the cascade;
  --   embsys   an embedded system for the top dimension;
  --   sols     solutions of the system embsys;
  --   filter   if true, then junk points will be removed,
  --            otherwise, the output will be superwitness sets.
  --   factor   if true and filter, then the filtered witness sets will be
  --            factored into irreducible components,
  --            otherwise, the output sets may still be reducible.

  -- ON RETURN :
  --   pathcnt  number of paths tracked at each level in the cascade;
  --   filtcnt  if filter, then filtcnt contains the counts of each stage
  --            of removing the junk points from the witness supersets;
  --   idxfac   indices to the irreducible factors in the decomposition,
  --            if filter and factor are both set to true on input.

end Running_Cascades;
