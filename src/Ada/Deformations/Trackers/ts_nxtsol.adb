with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;   use Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;   use DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;   use QuadDobl_System_and_Solutions_io;
with Multprec_Floating_Numbers;
with Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;
with Multprec_System_and_Solutions_io;   use Multprec_System_and_Solutions_io;
with Standard_to_Multprec_Convertors;
with Standard_Path_Tracker;
with DoblDobl_Path_Tracker;
with QuadDobl_Path_Tracker;
with Multprec_Path_Tracker;

procedure ts_nxtsol is

-- DESCRIPTION :
--   Allows tracking of one path via a generator "next"
--   which gives the next solution point on the path.

  procedure Standard_Initialize_Path_Tracker is

  -- DESCRIPTION :
  --   Prompts the user for a target system, a start system,
  --   and start solutions and initializes the standard path tracker.

    use Standard_Complex_Poly_Systems,Standard_Complex_Solutions;

    tgt_sys,sta_sys : Link_to_Poly_Sys;
    sols : Solution_List;
    ans : character;

  begin
    new_line;
    put_line("Reading the target system ...");
    get(tgt_sys);
    new_line;
    put_line("Reading the start system and its solutions...");
    get(sta_sys,sols);
    new_line;
    put("Fixed gamma constant ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' 
     then Standard_Path_Tracker.Init(tgt_sys,sta_sys,true,Head_Of(sols));
     else Standard_Path_Tracker.Init(tgt_sys,sta_sys,false,Head_Of(sols));
    end if;
  end Standard_Initialize_Path_Tracker;

  procedure DoblDobl_Initialize_Path_Tracker is

  -- DESCRIPTION :
  --   Prompts the user for a target system, a start system,
  --   and start solutions and initializes the double double path tracker.

    use DoblDobl_Complex_Poly_Systems,DoblDobl_Complex_Solutions;

    tgt_sys,sta_sys : Link_to_Poly_Sys;
    sols : Solution_List;
    ans : character;

  begin
    new_line;
    put_line("Reading the target system ...");
    get(tgt_sys);
    new_line;
    put_line("Reading the start system and its solutions...");
    get(sta_sys,sols);
    new_line;
    put("Fixed gamma constant ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then DoblDobl_Path_Tracker.Init(tgt_sys,sta_sys,true,Head_Of(sols));
     else DoblDobl_Path_Tracker.Init(tgt_sys,sta_sys,false,Head_Of(sols));
    end if;
  end DoblDobl_Initialize_Path_Tracker;

  procedure QuadDobl_Initialize_Path_Tracker is

  -- DESCRIPTION :
  --   Prompts the user for a target system, a start system,
  --   and start solutions and initializes the quad double path tracker.

    use QuadDobl_Complex_Poly_Systems,QuadDobl_Complex_Solutions;

    tgt_sys,sta_sys : Link_to_Poly_Sys;
    sols : Solution_List;
    ans : character;

  begin
    new_line;
    put_line("Reading the target system ...");
    get(tgt_sys);
    new_line;
    put_line("Reading the start system and its solutions...");
    get(sta_sys,sols);
    new_line;
    put("Fixed gamma constant ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then QuadDobl_Path_Tracker.Init(tgt_sys,sta_sys,true,Head_Of(sols));
     else QuadDobl_Path_Tracker.Init(tgt_sys,sta_sys,false,Head_Of(sols));
    end if;
  end QuadDobl_Initialize_Path_Tracker;

  procedure Multprec_Initialize_Path_Tracker ( deci : in natural32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a target system, a start system,
  --   and start solutions and initializes the multiprecision path tracker.

    use Multprec_Complex_Poly_Systems,Multprec_Complex_Solutions;

    tgt_sys,sta_sys : Link_to_Poly_Sys;
    sols : Solution_List;
    ans : character;
    size : constant natural32
         := Multprec_Floating_Numbers.Decimal_to_Size(deci);

  begin
    new_line;
    put_line("Reading the target system ...");
    get(tgt_sys);
    Standard_to_Multprec_Convertors.Set_Size(tgt_sys.all,size);
    new_line;
    put_line("Reading the start system and its solutions...");
    get(sta_sys,sols);
    Standard_to_Multprec_Convertors.Set_Size(sta_sys.all,size);
    Multprec_Complex_Solutions.Set_Size(sols,size);
    new_line;
    put("Fixed gamma constant ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' 
     then Multprec_Path_Tracker.Init(tgt_sys,sta_sys,true,Head_Of(sols),deci);
     else Multprec_Path_Tracker.Init(tgt_sys,sta_sys,false,Head_Of(sols),deci);
    end if;
  end Multprec_Initialize_Path_Tracker;

  procedure Standard_Run_Path_Tracker is

  -- DESCRIPTION :
  --   Shows the current solution and prompts the user each time to make
  --   one extra predictor-corrector step to compute the next solution
  --   on the path with standard double precision arithmetic.

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

  -- DESCRIPTION :
  --   Shows the current solution and prompts the user each time to make
  --   one extra predictor-corrector step to compute the next solution
  --   on the path with double double precision arithmetic.
    
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

  -- DESCRIPTION :
  --   Shows the current solution and prompts the user each time to make
  --   one extra predictor-corrector step to compute the next solution
  --   on the path with quad double precision arithmetic.
    
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

  procedure Multprec_Run_Path_Tracker is

  -- DESCRIPTION :
  --   Shows the current solution and prompts the user each time to make
  --   one extra predictor-corrector step to compute the next solution
  --   on the path with multiprecision arithmetic.
    
    s : Multprec_Complex_Solutions.Link_to_Solution
      := Multprec_Path_Tracker.get_current;
    ans : character;

  begin
    new_line;
    put_line("The current solution : ");
    Multprec_Complex_Solutions_io.put(s.all); new_line;
    loop
      put("Do predictor-corrector step ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      s := Multprec_Path_Tracker.get_next;
      put_line("The current solution : ");
      Multprec_Complex_Solutions_io.put(s.all); new_line;
    end loop;
  end Multprec_Run_Path_Tracker;

  procedure Main is

    ans : character;
    deci : natural32 := 0;

  begin
    new_line;
    put_line("Tracking a path with a generator ...");
    put_line("  1. in complex standard floating-point arithmetic;");
    put_line("  2. in complex double double arithmetic;");
    put_line("  3. in complex quad double arithmetic;");
    put_line("  4. in complex multiprecision arithmetic.");
    put("Type 1, 2, 3, or 4 to select precision : ");
    Ask_Alternative(ans,"1234");
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
      when '4' =>
        new_line;
        put("Give the number of decimal places in the working precision : ");
        get(deci);
        Multprec_Initialize_Path_Tracker(deci);
        Multprec_Run_Path_Tracker;
        Multprec_Path_Tracker.Clear;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_nxtsol;
