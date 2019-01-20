with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Parameter_Systems;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;
with Solution_Drops;
with Standard_System_and_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Homotopy_Continuation_Parameters;
with Homotopy_Continuation_Parameters_io;
with Series_Path_Trackers;
with Standard_SeriesPade_Tracker;
with DoblDobl_SeriesPade_Tracker;
with QuadDobl_SeriesPade_Tracker;

procedure ts_nxtpadsol is

-- DESCRIPTION :
--   Interactive test on the get_next() method to track paths
--   with a series Pade predictor.

  procedure Standard_Loop
              ( sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   With the homotopy initialized properly,
  --   runs the step-by-step tracker, starting at the solutions in sols,
  --   in standard double precision.

    solsptr : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    ans : character;
    fail : boolean;

  begin
    loop
      ls := Standard_Complex_Solutions.Head_Of(solsptr);
      Standard_SeriesPade_Tracker.Init(ls);
      put_line("Checking the start solution ...");
      Standard_SeriesPade_Tracker.Correct(fail,true);
      if fail then
        put_line("The start solution is NOT okay!?");
      else
        put_line("The start solution is okay.");
        loop
          Standard_SeriesPade_Tracker.Predict_and_Correct(fail,true);
          if fail
           then put_line("Failed to meet the accuracy requirements.  Abort.");
          end if;
          exit when fail;
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end loop;
        ls := Standard_SeriesPade_Tracker.Get_Current_Solution;
        put_line("The solution : ");
        Standard_Complex_Solutions_io.put(ls.all); new_line;
      end if;
      solsptr := Standard_Complex_Solutions.Tail_Of(solsptr);
      exit when Standard_Complex_Solutions.Is_Null(solsptr);
      put("Continue to the next solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Loop;

  procedure DoblDobl_Loop
              ( sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   With the homotopy initialized properly,
  --   runs the step-by-step tracker, starting at the solutions in sols,
  --   in double double precision.

    solsptr : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    ans : character;
    fail : boolean;

  begin
    loop
      ls := DoblDobl_Complex_Solutions.Head_Of(solsptr);
      DoblDobl_SeriesPade_Tracker.Init(ls);
      put_line("Checking the start solution ...");
      DoblDobl_SeriesPade_Tracker.Correct(fail,true);
      if fail then
        put_line("The start solution is NOT okay!?");
      else
        put_line("The start solution is okay.");
        loop
          DoblDobl_SeriesPade_Tracker.Predict_and_Correct(fail,true);
          if fail
           then put_line("Failed to meet the accuracy requirements.  Abort.");
          end if;
          exit when fail;
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end loop;
        ls := DoblDobl_SeriesPade_Tracker.Get_Current_Solution;
        put_line("The solution : ");
        DoblDobl_Complex_Solutions_io.put(ls.all); new_line;
      end if;
      solsptr := DoblDobl_Complex_Solutions.Tail_Of(solsptr);
      exit when DoblDobl_Complex_Solutions.Is_Null(solsptr);
      put("Continue to the next solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end DoblDobl_Loop;

  procedure QuadDobl_Loop
              ( sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   With the homotopy initialized properly,
  --   runs the step-by-step tracker, starting at the solutions in sols,
  --   in quad double precision.

    solsptr : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    ans : character;
    fail : boolean;

  begin
    loop
      ls := QuadDobl_Complex_Solutions.Head_Of(sols);
      QuadDobl_SeriesPade_Tracker.Init(ls);
      put_line("Checking the start solution ...");
      QuadDobl_SeriesPade_Tracker.Correct(fail,true);
      if fail then
        put_line("The start solution is NOT okay!?");
      else
        put_line("The start solution is okay.");
        loop
          QuadDobl_SeriesPade_Tracker.Predict_and_Correct(fail,true);
          if fail
           then put_line("Failed to meet the accuracy requirements.  Abort.");
          end if;
          exit when fail;
          put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end loop;
        ls := QuadDobl_SeriesPade_Tracker.Get_Current_Solution;
        put_line("The solution : ");
        QuadDobl_Complex_Solutions_io.put(ls.all); new_line;
      end if;
      solsptr := QuadDobl_Complex_Solutions.Tail_Of(solsptr);
      exit when QuadDobl_Complex_Solutions.Is_Null(solsptr);
      put("Continue to the next solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end QuadDobl_Loop;
 
  procedure Standard_Main is

  -- DESCRIPTION :
  --   Allows the user to tune the homotopy continuation parameters,
  --   prompts for a target, start system, and start solutions,
  --   in standard double precision.

    nbq,nvr,idx,nbpar : integer32;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    target,start : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols,dropsols : Standard_Complex_Solutions.Solution_List;
    arthom : boolean;
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);

    use Standard_Parameter_Systems;

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
    arthom := Series_Path_Trackers.Prompt_for_Artificial;
    if arthom then
      new_line;
      put_line("Reading the target system ..."); get(target);
      new_line;
      put_line("Reading the start system and its solutions ...");
      Standard_System_and_Solutions_io.get(start,sols);
      Standard_SeriesPade_Tracker.Init(target,start);
      Standard_Loop(sols);
    else
      Read_Parameter_Homotopy(target,sols,nbq,nvr,nbpar);
      declare
        par : Standard_Integer_Vectors.Vector(1..nbpar);
      begin
        par := Define_Parameters(nbq,nvr,nbpar);
        idx := par(1);
      end;
      put("number of equations : "); put(nbq,1); new_line;
      put("number of variables : "); put(nvr,1); new_line;
      put("index of the continuation parameter : "); put(idx,1); new_line;
      Standard_SeriesPade_Tracker.Init(target,idx);
      dropsols := Solution_Drops.Drop(sols,natural32(idx));
      Standard_Complex_Solutions.Set_Continuation_Parameter(dropsols,zero);
      Standard_Loop(dropsols);
    end if;
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
    arthom : boolean;

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
    arthom := Series_Path_Trackers.Prompt_for_Artificial;
    if arthom then
      new_line;
      put_line("Reading the target system ..."); get(target);
      new_line;
      put_line("Reading the start system and its solutions ...");
      DoblDobl_System_and_Solutions_io.get(start,sols);
      DoblDobl_SeriesPade_Tracker.Init(target,start);
      DoblDobl_Loop(sols);
    end if;
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
    arthom : boolean;

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
    arthom := Series_Path_Trackers.Prompt_for_Artificial;
    if arthom then
      new_line;
      put_line("Reading the target system ..."); get(target);
      new_line;
      put_line("Reading the start system and its solutions ...");
      QuadDobl_System_and_Solutions_io.get(start,sols);
      QuadDobl_SeriesPade_Tracker.Init(target,start);
      QuadDobl_Loop(sols);
    end if;
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
