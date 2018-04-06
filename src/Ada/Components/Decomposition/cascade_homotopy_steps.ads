with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Cascade_Homotopy_Steps is

-- DESCRIPTION :
--   One step in the cascade homotopy removes the last hyperplane of
--   the embedded system.  The embedded system itself is used as the start
--   system in the homotopy to the system with the last hyperplane removed.
--   There are twelve (3x2x2) procedures in this package:
--   1) for double, double double, and quad double precision;
--   2) for ordinary and Laurent polynomial systems;
--   3) with and without output to file.
--   Every procedure allows for multitasking.

  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                sols0,sols1 : out Standard_Complex_Solutions.Solution_List;
                time : out duration );
  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                sols0,sols1 : out Standard_Complex_Solutions.Solution_List;
                time : out duration );
  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out DoblDobl_Complex_Solutions.Solution_List;
                time : out duration );
  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out DoblDobl_Complex_Solutions.Solution_List;
                time : out duration );
  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out QuadDobl_Complex_Solutions.Solution_List;
                time : out duration );
  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out QuadDobl_Complex_Solutions.Solution_List;
                time : out duration );

  -- DESCRIPTION :
  --   Runs a continuation to remove the slice from the embedded
  --   (Laurent) polynomial system without output,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   nt       the number of tasks, for no tasking set nt = 0;
  --   embsys   the embedded (Laurent) polynomial system;
  --   level    equals the current dimension in the cascade;
  --   zerotol  is the tolerance whether the slack variable is zero or not;
  --   tolsing  tolerance on rco to remove singular solutions with
  --            nonzero slack variable value;
  --   sols     solutions with nonzero slack variables of embsys.

  -- ON RETURN :
  --   sols     solutions at the end of the homotopy;
  --   sols0    solutions with zero slack variable;
  --   sols1    solutions with nonzero slack variable;
  --   time     the elapsed CPU time.

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                sols0,sols1 : out Standard_Complex_Solutions.Solution_List;
                time : out duration );
  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                sols0,sols1 : out Standard_Complex_Solutions.Solution_List;
                time : out duration );
  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out DoblDobl_Complex_Solutions.Solution_List;
                time : out duration );
  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out DoblDobl_Complex_Solutions.Solution_List;
                time : out duration );
  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out QuadDobl_Complex_Solutions.Solution_List;
                time : out duration );
  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out QuadDobl_Complex_Solutions.Solution_List;
                time : out duration );

  -- DESCRIPTION :
  --   Runs a continuation to remove the slice from the embedded (Laurent)
  --   polynomial system with output to file,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nt       the number of tasks, for no tasking set nt = 0;
  --   embsys   the embedded Laurent polynomial system;
  --   level    equals the current dimension in the cascade;
  --   zerotol  is the tolerance whether the slack variable is zero or not;
  --   tolsing  tolerance on rco to remove singular solutions with
  --            nonzero slack variable value;
  --   sols     solutions with nonzero slack variables of embsys.

  -- ON RETURN :
  --   sols     solutions at the end of the homotopy;
  --   sols0    solutions with zero slack variable;
  --   sols1    solutions with nonzero slack variable;
  --   time     the elapsed CPU time.

end Cascade_Homotopy_Steps;
