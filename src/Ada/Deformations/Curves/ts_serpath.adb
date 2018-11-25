with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Characters_and_Numbers;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_cv;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;
with Complex_Series_and_Polynomials_io;  use Complex_Series_and_Polynomials_io;
with Series_and_Homotopies;
with Series_and_Trackers;
with Homotopy_Series_Readers;
with Homotopy_Continuation_Parameters;
with Homotopy_Continuation_Parameters_io;
with Drivers_to_Series_Trackers;         use Drivers_to_Series_Trackers;

procedure ts_serpath is

-- DESCRIPTION :
--   Developing path trackers with power series.

  procedure Standard_Test
              ( nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   With a homotopy defined in Standard_Homotopy of nq equations,
  --   and start solutions in sols, test the path tracking,
  --   in standard double precision.

    use Standard_Complex_Solutions;

    h : Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);
    p : Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    ans : character;
    verbose,tofile : boolean;
    file : file_type;
    timer : Timing_Widget;
    prevgamma : Standard_Complex_Numbers.Complex_Number;

  begin
   -- put_line("The homotopy system :"); put_line(h);
   -- put_line("The series system :"); put(s,1); new_line;
    p.gamma := Standard_Homotopy.Accessibility_Constant;
    prevgamma := p.gamma;
    Homotopy_Continuation_Parameters_io.Tune(p);
    if not Standard_Complex_Numbers.Equal(p.gamma,prevgamma)
     then Standard_Reset_Gamma(p.gamma);
    end if;
    Set_Output(file,verbose,tofile);
    if tofile
     then Homotopy_Continuation_Parameters_io.put(file,p); flush(file);
    end if;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      put("Tracking path "); put(i,1); put_line(" ...");
      if tofile then
        Series_and_Trackers.Track_One_Path(file,s,ls.all,p,verbose);
        put(file,"Solution "); put(file,i,1); put_line(file," :");
        put(file,ls.all); new_line(file);
      else
        Series_and_Trackers.Track_One_Path(standard_output,s,ls.all,p,verbose);
        put("Solution "); put(i,1); put_line(" :"); put(ls.all);
        put("Continue to the next path ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end if;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    if tofile then
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(file,p.numdeg,p.dendeg,0,timer);
      Refine_Roots(file,nq,sols);
    else
      put_line("THE SOLUTIONS :");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(standard_output,p.numdeg,p.dendeg,0,timer);
      Refine_Roots(standard_output,nq,sols);
    end if;
    Standard_Complex_Poly_Systems.Clear(h);
    Standard_CSeries_Poly_Systems.Clear(s);
  end Standard_Test;

  procedure DoblDobl_Test
              ( nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   With a homotopy defined in Standard_Homotopy of nq equations,
  --   and start solutions in sols, test the path tracking,
  --   in double double precision.

    use DoblDobl_Complex_Solutions;

    h : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Homotopy_System;
    s : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);
    p : Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    ans : character;
    verbose,tofile : boolean;
    file : file_type;
    timer : Timing_Widget;
    ddgamma : constant DoblDobl_Complex_Numbers.Complex_Number
            := DoblDobl_Homotopy.Accessibility_Constant;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := DoblDobl_Complex_Numbers_cv.DoblDobl_Complex_to_Standard(ddgamma);
    prevgamma : Standard_Complex_Numbers.Complex_Number;

  begin
   -- put_line("The homotopy system :"); put_line(h);
   -- put_line("The series system :"); put(s,1); new_line;
    p.gamma := gamma;
    prevgamma := p.gamma;
    Homotopy_Continuation_Parameters_io.Tune(p);
    if not Standard_Complex_Numbers.Equal(p.gamma,prevgamma)
     then DoblDobl_Reset_Gamma(p.gamma);
    end if;
    Set_Output(file,verbose,tofile);
    if tofile
     then Homotopy_Continuation_Parameters_io.put(file,p); flush(file);
    end if;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      put("Tracking path "); put(i,1); put_line(" ...");
      if tofile then
        Series_and_Trackers.Track_One_Path(file,s,ls.all,p,verbose);
        put(file,"Solution "); put(file,i,1); put_line(file," :");
        put(file,ls.all); new_line(file);
      else
        Series_and_Trackers.Track_One_Path(standard_output,s,ls.all,p,verbose);
        put("Solution "); put(i,1); put_line(" :"); put(ls.all);
        put("Continue to the next path ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end if;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    if tofile then
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(file,p.numdeg,p.dendeg,1,timer);
      Refine_Roots(file,nq,sols);
    else
      put_line("THE SOLUTIONS :");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(standard_output,p.numdeg,p.dendeg,1,timer);
      Refine_Roots(standard_output,nq,sols);
    end if;
    DoblDobl_Complex_Poly_Systems.Clear(h);
    DoblDobl_CSeries_Poly_Systems.Clear(s);
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   With a homotopy defined in Standard_Homotopy of nq equations,
  --   and start solutions in sols, test the path tracking,
  --   in quad double precision.

    use QuadDobl_Complex_Solutions;

    h : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Homotopy_System;
    s : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1,false);
    p : Homotopy_Continuation_Parameters.Parameters
      := Homotopy_Continuation_Parameters.Default_Values;
    tmp : Solution_List := sols;
    len : constant integer32 := integer32(Length_Of(sols));
    ls : Link_to_Solution;
    ans : character;
    verbose,tofile : boolean;
    file : file_type;
    timer : Timing_Widget;
    qdgamma : constant QuadDobl_Complex_Numbers.Complex_Number
            := QuadDobl_Homotopy.Accessibility_Constant;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := QuadDobl_Complex_Numbers_cv.QuadDobl_Complex_to_Standard(qdgamma);
    prevgamma : Standard_Complex_Numbers.Complex_Number;

  begin
   -- put_line("The homotopy system :"); put_line(h);
   -- put_line("The series system :"); put(s,1); new_line;
    p.gamma := gamma;
    prevgamma := p.gamma;
    Homotopy_Continuation_Parameters_io.Tune(p);
    if not Standard_Complex_Numbers.Equal(p.gamma,prevgamma)
     then QuadDobl_Reset_Gamma(p.gamma);
    end if;
    Set_Output(file,verbose,tofile);
    if tofile
     then Homotopy_Continuation_Parameters_io.put(file,p); flush(file);
    end if;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      put("Tracking path "); put(i,1); put_line(" ...");
      if tofile then
        Series_and_Trackers.Track_One_Path(file,s,ls.all,p,verbose);
        put(file,"Solution "); put(file,i,1); put_line(file," :");
        put(file,ls.all); new_line(file);
      else
        Series_and_Trackers.Track_One_Path(standard_output,s,ls.all,p,verbose);
        put("Solution "); put(i,1); put_line(" :"); put(ls.all);
        put("Continue to the next path ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end if;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    if tofile then
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(file,p.numdeg,p.dendeg,2,timer);
      Refine_Roots(file,nq,sols);
    else
      put_line("THE SOLUTIONS :");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(standard_output,p.numdeg,p.dendeg,2,timer);
      Refine_Roots(standard_output,nq,sols);
    end if;
    QuadDobl_Complex_Poly_Systems.Clear(h);
    QuadDobl_CSeries_Poly_Systems.Clear(s);
  end QuadDobl_Test;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in standard double precision.

    nbeq : integer32;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Homotopy_Series_Readers.Standard_Reader(nbeq,sols);
    new_line;
    Standard_Test(nbeq,sols);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in double double precision.

    nbeq : integer32;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Homotopy_Series_Readers.DoblDobl_Reader(nbeq,sols);
    new_line;
    DoblDobl_Test(nbeq,sols);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in quad double precision.

    nbeq : integer32;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Homotopy_Series_Readers.QuadDobl_Reader(nbeq,sols);
    new_line;
    QuadDobl_Test(nbeq,sols);
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision
  --   and then launches the test.

    ans : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
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
end ts_serpath;
