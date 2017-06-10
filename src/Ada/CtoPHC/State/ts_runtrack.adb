with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with PHCpack_Operations;
with Crude_Path_Trackers;

procedure ts_runtrack is

-- DESCRIPTION :
--   Development of path trackers which write the results immediately
--   in the solution container.

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts the user for a target and a start system
  --   with start solution and then launches the tracking,
  --   in standard double precision.

    lp,lq : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading a target system ...");
    get(lp);
    PHCpack_Operations.Store_Target_System(lp.all);
    Standard_System_and_Solutions_io.get(lq,sols);
    PHCpack_Operations.Store_Start_System(lq.all);
    PHCpack_Operations.Store_Start_Solutions(sols);
    new_line;
    put("Read ");
    put(Standard_Complex_Solutions.Length_Of(sols),1);
    put_line(" start solutions.");
    new_line;
    PHCpack_Operations.Create_Standard_Homotopy;
    Crude_Path_Trackers.Standard_Track_Paths(true);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a target and a start system
  --   with start solution and then launches the tracking,
  --   in double double precision.

    lp,lq : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading a target system ...");
    get(lp);
    PHCpack_Operations.Store_Target_System(lp.all);
    DoblDobl_System_and_Solutions_io.get(lq,sols);
    PHCpack_Operations.Store_Start_System(lq.all);
    PHCpack_Operations.Store_Start_Solutions(sols);
    new_line;
    put("Read ");
    put(DoblDobl_Complex_Solutions.Length_Of(sols),1);
    put_line(" start solutions.");
    new_line;
    PHCpack_Operations.Create_DoblDobl_Homotopy;
    Crude_Path_Trackers.DoblDobl_Track_Paths(true);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a target and a start system
  --   with start solution and then launches the tracking,
  --   in quad double precision.

    lp,lq : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Reading a target system ...");
    get(lp);
    PHCpack_Operations.Store_Target_System(lp.all);
    QuadDobl_System_and_Solutions_io.get(lq,sols);
    PHCpack_Operations.Store_Start_System(lq.all);
    PHCpack_Operations.Store_Start_Solutions(sols);
    new_line;
    put("Read ");
    put(QuadDobl_Complex_Solutions.Length_Of(sols),1);
    put_line(" start solutions.");
    new_line;
    PHCpack_Operations.Create_QuadDobl_Homotopy;
    Crude_Path_Trackers.QuadDobl_Track_Paths(true);
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision
  --   and then launches the proper path tracker.

    ans : character;

  begin
    new_line;
    put_line("Calling the crude path trackers ...");
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_runtrack;
