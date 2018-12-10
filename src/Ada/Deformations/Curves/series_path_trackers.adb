with Ada.Calendar;
with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Characters_and_Numbers;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_cv;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
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

package body Series_Path_Trackers is

  procedure Standard_Run
              ( nq : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

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
    monitor,verbose,tofile : boolean;
    file : file_type;
    timer : Timing_Widget;
    prevgamma : Standard_Complex_Numbers.Complex_Number;
    nbrsteps,minnbrsteps,maxnbrsteps : natural32;
    nbrcorrs,minnbrcorrs,maxnbrcorrs : natural32;
    minsize,maxsize,smallest,largest : double_float;
    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;

  begin
   -- put_line("The homotopy system :"); put_line(h);
   -- put_line("The series system :"); put(s,1); new_line;
    p.gamma := Standard_Homotopy.Accessibility_Constant;
    prevgamma := p.gamma;
    Homotopy_Continuation_Parameters_io.Tune(p);
    if not Standard_Complex_Numbers.Equal(p.gamma,prevgamma)
     then Standard_Reset_Gamma(p.gamma);
    end if;
    Set_Output(file,monitor,verbose,tofile);
    if tofile
     then Homotopy_Continuation_Parameters_io.put(file,p); flush(file);
    end if;
    minnbrsteps := p.maxsteps+1; maxnbrsteps := 0;
    minnbrcorrs := (p.maxsteps+1)*p.corsteps+1; maxnbrcorrs := 0;
    smallest := p.maxsize; largest := 0.0;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      if monitor
       then put("Tracking path "); put(i,1); put_line(" ...");
      end if;
      if tofile then
        Series_and_Trackers.Track_One_Path
          (file,s,ls.all,p,nbrsteps,nbrcorrs,minsize,maxsize,verbose);
        if verbose then
          Series_and_Trackers.Write_Path_Statistics
            (file,nbrsteps,nbrcorrs,minsize,maxsize);
        end if;
        put(file,"Solution "); put(file,i,1); put_line(file," :");
        put(file,ls.all); new_line(file);
      else
        Series_and_Trackers.Track_One_Path
          (standard_output,s,ls.all,p,
           nbrsteps,nbrcorrs,minsize,maxsize,verbose);
        if verbose then
          Series_and_Trackers.Write_Path_Statistics
            (standard_output,nbrsteps,nbrcorrs,minsize,maxsize);
        end if;
        put("Solution "); put(i,1); put_line(" :"); put(ls.all);
        put("Continue to the next path ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end if;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      Series_and_Trackers.Update_Counters(minnbrsteps,maxnbrsteps,nbrsteps);
      Series_and_Trackers.Update_Counters(minnbrcorrs,maxnbrcorrs,nbrcorrs);
      Series_and_Trackers.Update_MinMax(smallest,largest,minsize,maxsize);
    end loop;
    tstop(timer);
    if tofile then
      Series_and_Trackers.Write_Total_Path_Statistics
        (file,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs,smallest,largest);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(file,p.numdeg,p.dendeg,0,timer);
      Refine_Roots(file,nq,sols);
      Write_Conclusion(file,start_moment);
    else
      Series_and_Trackers.Write_Total_Path_Statistics
        (standard_output,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs,
         smallest,largest);
      new_line;
      put_line("THE SOLUTIONS :");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(standard_output,p.numdeg,p.dendeg,0,timer);
      Refine_Roots(standard_output,nq,sols);
      Write_Conclusion(standard_output,start_moment);
    end if;
    Standard_Complex_Poly_Systems.Clear(h);
    Standard_CSeries_Poly_Systems.Clear(s);
  end Standard_Run;

  procedure DoblDobl_Run
              ( nq : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

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
    monitor,verbose,tofile : boolean;
    file : file_type;
    timer : Timing_Widget;
    ddgamma : constant DoblDobl_Complex_Numbers.Complex_Number
            := DoblDobl_Homotopy.Accessibility_Constant;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := DoblDobl_Complex_Numbers_cv.DoblDobl_Complex_to_Standard(ddgamma);
    prevgamma : Standard_Complex_Numbers.Complex_Number;
    nbrsteps,minnbrsteps,maxnbrsteps : natural32;
    nbrcorrs,minnbrcorrs,maxnbrcorrs : natural32;
    minsize,maxsize,smallest,largest : double_float;
    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;

  begin
   -- put_line("The homotopy system :"); put_line(h);
   -- put_line("The series system :"); put(s,1); new_line;
    p.gamma := gamma;
    prevgamma := p.gamma;
    Homotopy_Continuation_Parameters_io.Tune(p);
    if not Standard_Complex_Numbers.Equal(p.gamma,prevgamma)
     then DoblDobl_Reset_Gamma(p.gamma);
    end if;
    Set_Output(file,monitor,verbose,tofile);
    if tofile
     then Homotopy_Continuation_Parameters_io.put(file,p); flush(file);
    end if;
    minnbrsteps := p.maxsteps+1; maxnbrsteps := 0;
    minnbrcorrs := (p.maxsteps+1)*p.corsteps+1; maxnbrcorrs := 0;
    smallest := p.maxsize; largest := 0.0;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      if monitor
       then put("Tracking path "); put(i,1); put_line(" ...");
      end if;
      if tofile then
        Series_and_Trackers.Track_One_Path
          (file,s,ls.all,p,nbrsteps,nbrcorrs,minsize,maxsize,verbose);
        if verbose then
          Series_and_Trackers.Write_Path_Statistics
            (file,nbrsteps,nbrcorrs,minsize,maxsize);
        end if;
        put(file,"Solution "); put(file,i,1); put_line(file," :");
        put(file,ls.all); new_line(file);
      else
        Series_and_Trackers.Track_One_Path
          (standard_output,s,ls.all,p,
           nbrsteps,nbrcorrs,minsize,maxsize,verbose);
        if verbose then
          Series_and_Trackers.Write_Path_Statistics
            (standard_output,nbrsteps,nbrcorrs,minsize,maxsize);
        end if;
        put("Solution "); put(i,1); put_line(" :"); put(ls.all);
        put("Continue to the next path ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end if;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      Series_and_Trackers.Update_Counters(minnbrsteps,maxnbrsteps,nbrsteps);
      Series_and_Trackers.Update_Counters(minnbrcorrs,maxnbrcorrs,nbrcorrs);
      Series_and_Trackers.Update_MinMax(smallest,largest,minsize,maxsize);
    end loop;
    tstop(timer);
    if tofile then
      Series_and_Trackers.Write_Total_Path_Statistics
        (file,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs,smallest,largest);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(file,p.numdeg,p.dendeg,1,timer);
      Refine_Roots(file,nq,sols);
      Write_Conclusion(file,start_moment);
    else
      Series_and_Trackers.Write_Total_Path_Statistics
        (standard_output,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs,
         smallest,largest);
      new_line;
      put_line("THE SOLUTIONS :");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(standard_output,p.numdeg,p.dendeg,1,timer);
      Refine_Roots(standard_output,nq,sols);
      Write_Conclusion(standard_output,start_moment);
    end if;
    DoblDobl_Complex_Poly_Systems.Clear(h);
    DoblDobl_CSeries_Poly_Systems.Clear(s);
  end DoblDobl_Run;

  procedure QuadDobl_Run
              ( nq : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

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
    monitor,verbose,tofile : boolean;
    file : file_type;
    timer : Timing_Widget;
    qdgamma : constant QuadDobl_Complex_Numbers.Complex_Number
            := QuadDobl_Homotopy.Accessibility_Constant;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := QuadDobl_Complex_Numbers_cv.QuadDobl_Complex_to_Standard(qdgamma);
    prevgamma : Standard_Complex_Numbers.Complex_Number;
    nbrsteps,minnbrsteps,maxnbrsteps : natural32;
    nbrcorrs,minnbrcorrs,maxnbrcorrs : natural32;
    minsize,maxsize,smallest,largest : double_float;
    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;

  begin
   -- put_line("The homotopy system :"); put_line(h);
   -- put_line("The series system :"); put(s,1); new_line;
    p.gamma := gamma;
    prevgamma := p.gamma;
    Homotopy_Continuation_Parameters_io.Tune(p);
    if not Standard_Complex_Numbers.Equal(p.gamma,prevgamma)
     then QuadDobl_Reset_Gamma(p.gamma);
    end if;
    Set_Output(file,monitor,verbose,tofile);
    if tofile
     then Homotopy_Continuation_Parameters_io.put(file,p); flush(file);
    end if;
    minnbrsteps := p.maxsteps+1; maxnbrsteps := 0;
    minnbrcorrs := (p.maxsteps+1)*p.corsteps+1; maxnbrcorrs := 0;
    smallest := p.maxsize; largest := 0.0;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      if monitor
       then put("Tracking path "); put(i,1); put_line(" ...");
      end if;
      if tofile then
        Series_and_Trackers.Track_One_Path
          (file,s,ls.all,p,nbrsteps,nbrcorrs,minsize,maxsize,verbose);
        if verbose then
          Series_and_Trackers.Write_Path_Statistics
            (file,nbrsteps,nbrcorrs,minsize,maxsize);
        end if;
        put(file,"Solution "); put(file,i,1); put_line(file," :");
        put(file,ls.all); new_line(file);
      else
        Series_and_Trackers.Track_One_Path
          (standard_output,s,ls.all,p,
           nbrsteps,nbrcorrs,minsize,maxsize,verbose);
        if verbose then
          Series_and_Trackers.Write_Path_Statistics
            (standard_output,nbrsteps,nbrcorrs,minsize,maxsize);
        end if;
        put("Solution "); put(i,1); put_line(" :"); put(ls.all);
        put("Continue to the next path ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end if;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
      Series_and_Trackers.Update_Counters(minnbrsteps,maxnbrsteps,nbrsteps);
      Series_and_Trackers.Update_Counters(minnbrcorrs,maxnbrcorrs,nbrcorrs);
      Series_and_Trackers.Update_MinMax(smallest,largest,minsize,maxsize);
    end loop;
    tstop(timer);
    if tofile then
      Series_and_Trackers.Write_Total_Path_Statistics
        (file,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs,smallest,largest);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(file,p.numdeg,p.dendeg,2,timer);
      Refine_Roots(file,nq,sols);
      Write_Conclusion(file,start_moment);
    else
      Series_and_Trackers.Write_Total_Path_Statistics
        (standard_output,minnbrsteps,maxnbrsteps,minnbrcorrs,maxnbrcorrs,
         smallest,largest);
      new_line;
      put_line("THE SOLUTIONS :");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Write_Timer(standard_output,p.numdeg,p.dendeg,2,timer);
      Refine_Roots(standard_output,nq,sols);
      Write_Conclusion(standard_output,start_moment);
    end if;
    QuadDobl_Complex_Poly_Systems.Clear(h);
    QuadDobl_CSeries_Poly_Systems.Clear(s);
  end QuadDobl_Run;

  procedure Standard_Main is

    nbeq : integer32;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Homotopy_Series_Readers.Standard_Reader(nbeq,sols);
    new_line;
    Standard_Run(nbeq,sols);
  end Standard_Main;

  procedure DoblDobl_Main is

    nbeq : integer32;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Homotopy_Series_Readers.DoblDobl_Reader(nbeq,sols);
    new_line;
    DoblDobl_Run(nbeq,sols);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

    nbeq : integer32;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Homotopy_Series_Readers.QuadDobl_Reader(nbeq,sols);
    new_line;
    QuadDobl_Run(nbeq,sols);
  end QuadDobl_Main;

end Series_Path_Trackers;
