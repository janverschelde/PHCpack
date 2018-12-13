with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Homotopy_Continuation_Parameters;
with Homotopy_Continuation_Parameters_io;
with Standard_SeriesPade_Tracker;
with DoblDobl_SeriesPade_Tracker;
with QuadDobl_SeriesPade_Tracker;

procedure ts_nxtpadsol is

-- DESCRIPTION :
--   Interactive test on the get_next() method to track paths
--   with a series Pade predictor.

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Allows the user to tune the homotopy continuation parameters,
  --   prompts for a target, start system, and start solutions,
  --   in standard double precision.

    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    target,start : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    new_line;
    put_line("Tuning the homotopy continuation parameters ...");
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    Standard_SeriesPade_Tracker.Init(pars);
    new_line;
    put_line("The stored values of the parameters : ");
    new_line;
    Homotopy_Continuation_Parameters_io.put
      (Standard_SeriesPade_Tracker.Get_Parameters.all);
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and its solutions ...");
    Standard_System_and_Solutions_io.get(start,sols);
    Standard_SeriesPade_Tracker.Init(target,start);
    ls := Standard_Complex_Solutions.Head_Of(sols);
    Standard_SeriesPade_Tracker.Init(ls);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Allows the user to tune the homotopy continuation parameters,
  --   prompts for a target, start system, and start solutions,
  --   in double double precision.

    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    target,start : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    new_line;
    put_line("Tuning the homotopy continuation parameters ...");
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    DoblDobl_SeriesPade_Tracker.Init(pars);
    new_line;
    put_line("The stored values of the parameters : ");
    new_line;
    Homotopy_Continuation_Parameters_io.put
      (DoblDobl_SeriesPade_Tracker.Get_Parameters.all);
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and its solutions ...");
    DoblDobl_System_and_Solutions_io.get(start,sols);
    DoblDobl_SeriesPade_Tracker.Init(target,start);
    ls := DoblDobl_Complex_Solutions.Head_Of(sols);
    DoblDobl_SeriesPade_Tracker.Init(ls);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Allows the user to tune the homotopy continuation parameters,
  --   prompts for a target, start system, and start solutions,
  --   in quad double precision.

    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    target,start : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    new_line;
    put_line("Tuning the homotopy continuation parameters ...");
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    QuadDobl_SeriesPade_Tracker.Init(pars);
    new_line;
    put_line("The stored values of the parameters : ");
    new_line;
    Homotopy_Continuation_Parameters_io.put
      (QuadDobl_SeriesPade_Tracker.Get_Parameters.all);
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system and its solutions ...");
    QuadDobl_System_and_Solutions_io.get(start,sols);
    QuadDobl_SeriesPade_Tracker.Init(target,start);
    ls := QuadDobl_Complex_Solutions.Head_Of(sols);
    QuadDobl_SeriesPade_Tracker.Init(ls);
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision and then launches
  --   the corresponding test.

    ans : character;

  begin
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
end ts_nxtpadsol;
