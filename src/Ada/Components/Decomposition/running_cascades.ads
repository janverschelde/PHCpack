with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
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
                filter,factor : in boolean );
  procedure Standard_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean );
  procedure DoblDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean );
  procedure DoblDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean );
  procedure QuadDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean );
  procedure QuadDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean );

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
  --            otherwise, the output sets may still be reducible.

  procedure Standard_Run_Cascade
              ( file : in file_type; nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean );
  procedure Standard_Run_Cascade
              ( file : in file_type; nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean );
  procedure DoblDobl_Run_Cascade
              ( file : in file_type; nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean );
  procedure DoblDobl_Run_Cascade
              ( file : in file_type; nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean );
  procedure QuadDobl_Run_Cascade
              ( file : in file_type; nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean );
  procedure QuadDobl_Run_Cascade
              ( file : in file_type; nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean );

  -- DESCRIPTION :
  --   Given an embedding of the top dimensional solution set,
  --   runs a cascade of homotopies, in standard double, double double,
  --   or quad double precision.  All output is written to file.

  -- ON ENTRY :
  --   file     file opened for output;
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

end Running_Cascades;
