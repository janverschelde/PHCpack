with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;   use Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;   use DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;   use QuadDobl_System_and_Solutions_io;
with Standard_Path_Tracker;
with DoblDobl_Path_Tracker;
with QuadDobl_Path_Tracker;

procedure ts_nxtsol is

-- DESCRIPTION :
--   Allows tracking of one path via a generator "next"
--   which gives the next solution point on the path.

  procedure Standard_Initialize_Path_Tracker is

    use Standard_Complex_Poly_Systems,Standard_Complex_Solutions;

    tgt_sys,sta_sys : Link_to_Poly_Sys;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading the target system ...");
    get(tgt_sys);
    new_line;
    put_line("Reading the start system and its solutions...");
    get(sta_sys,sols);
    Standard_Path_Tracker.Init(tgt_sys,sta_sys,Head_Of(sols));
  end Standard_Initialize_Path_Tracker;

  procedure DoblDobl_Initialize_Path_Tracker is

    use DoblDobl_Complex_Poly_Systems,DoblDobl_Complex_Solutions;

    tgt_sys,sta_sys : Link_to_Poly_Sys;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading the target system ...");
    get(tgt_sys);
    new_line;
    put_line("Reading the start system and its solutions...");
    get(sta_sys,sols);
    DoblDobl_Path_Tracker.Init(tgt_sys,sta_sys,Head_Of(sols));
  end DoblDobl_Initialize_Path_Tracker;

  procedure QuadDobl_Initialize_Path_Tracker is

    use QuadDobl_Complex_Poly_Systems,QuadDobl_Complex_Solutions;

    tgt_sys,sta_sys : Link_to_Poly_Sys;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading the target system ...");
    get(tgt_sys);
    new_line;
    put_line("Reading the start system and its solutions...");
    get(sta_sys,sols);
    QuadDobl_Path_Tracker.Init(tgt_sys,sta_sys,Head_Of(sols));
  end QuadDobl_Initialize_Path_Tracker;

  procedure Standard_Run_Path_Tracker is
    
    s : Standard_Complex_Solutions.Link_to_Solution
      := Standard_Path_Tracker.get_current;
    ans : character;

  begin
    new_line;
    put_line("The current solution : ");
    Standard_Complex_Solutions_io.put(s.all); new_line;
    loop
      put("Do predictor-corrector step ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      s := Standard_Path_Tracker.get_next;
      put_line("The current solution : ");
      Standard_Complex_Solutions_io.put(s.all); new_line;
    end loop;
  end Standard_Run_Path_Tracker;

  procedure DoblDobl_Run_Path_Tracker is
    
    s : DoblDobl_Complex_Solutions.Link_to_Solution
      := DoblDobl_Path_Tracker.get_current;
    ans : character;

  begin
    new_line;
    put_line("The current solution : ");
    DoblDobl_Complex_Solutions_io.put(s.all); new_line;
    loop
      put("Do predictor-corrector step ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      s := DoblDobl_Path_Tracker.get_next;
      put_line("The current solution : ");
      DoblDobl_Complex_Solutions_io.put(s.all); new_line;
    end loop;
  end DoblDobl_Run_Path_Tracker;

  procedure QuadDobl_Run_Path_Tracker is
    
    s : QuadDobl_Complex_Solutions.Link_to_Solution
      := QuadDobl_Path_Tracker.get_current;
    ans : character;

  begin
    new_line;
    put_line("The current solution : ");
    QuadDobl_Complex_Solutions_io.put(s.all); new_line;
    loop
      put("Do predictor-corrector step ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      s := QuadDobl_Path_Tracker.get_next;
      put_line("The current solution : ");
      QuadDobl_Complex_Solutions_io.put(s.all); new_line;
    end loop;
  end QuadDobl_Run_Path_Tracker;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Tracking a path with a generator ...");
    put_line("  1. in complex standard floating-point arithmetic;");
    put_line("  2. in complex double double arithmetic;");
    put_line("  3. in complex quad double arithmetic;");
    put("Type 1, 2, or 3 to select precision : ");
    Ask_Alternative(ans,"123");
    case ans is
      when '1' =>
        Standard_Initialize_Path_Tracker;
        Standard_Run_Path_Tracker;
        Standard_Path_Tracker.Clear;
      when '2' =>
        DoblDobl_Initialize_Path_Tracker;
        DoblDobl_Run_Path_Tracker;
        DoblDobl_Path_Tracker.Clear;
      when '3' =>
        QuadDobl_Initialize_Path_Tracker;
        QuadDobl_Run_Path_Tracker;
        QuadDobl_Path_Tracker.Clear;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_nxtsol;
