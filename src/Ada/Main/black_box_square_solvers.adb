with text_io;                            use text_io;
with Timing_Package,Time_Stamps;         use Timing_Package,Time_Stamps;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;
with Write_Seed_Number;
with Write_Number_of_Tasks;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Scaling;
with DoblDobl_Scaling;
with QuadDobl_Scaling;
with Black_Box_Binomial_Solvers;         use Black_Box_Binomial_Solvers;
with Black_Box_Simplex_Solvers;          use Black_Box_Simplex_Solvers;
with Standard_Monomial_Maps;
with Standard_Monomial_Maps_io;          use Standard_Monomial_Maps_io;
with DoblDobl_Monomial_Maps;
with QuadDobl_Monomial_Maps;
with Black_Box_Root_Counters;            use Black_Box_Root_Counters;
with Standard_BlackBox_Continuations;    use Standard_BlackBox_Continuations;
with DoblDobl_Blackbox_Continuations;    use DoblDobl_Blackbox_Continuations;
with QuadDobl_Blackbox_Continuations;    use QuadDobl_Blackbox_Continuations;
with Greeting_Banners;
with Black_Box_Helpers;

package body Black_Box_Square_Solvers is

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                    append_sols : in boolean; verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
    use Standard_Monomial_Maps;

    timer : Timing_Widget;
    ended_moment : Ada.Calendar.Time;
    outfile : file_type;
    output_to_file : boolean;
    fail : boolean;
    sols : Solution_List;
    maps : Link_to_Array_of_Monomial_Map_Lists;
    rc : natural32;
    roco,hoco,poco,total : duration := 0.0;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solver_cases.Square_Main for Laurent systems,");
      put_line("in double precision ...");
    end if;
    Black_Box_Helpers.Ask_Output_File(outfile,outfilename,output_to_file);
    if output_to_file
     then put(outfile,natural32(p'last),p.all);
    end if;
    tstart(timer);
    if output_to_file
     then Black_Box_Simplex_Solver(outfile,p.all,sols,fail,verbose-1);
     else Black_Box_Simplex_Solver(p.all,sols,fail,verbose-1);
    end if;
    if fail or (Length_Of(sols) = 0) then
      if output_to_file
       then Black_Box_Binomial_Solver(outfile,p.all,maps,fail,verbose-1);
       else Black_Box_Binomial_Solver(p.all,false,maps,fail,verbose-1);
      end if;
    end if;
    if fail or (Length_Of(sols) = 0) then
      declare
        pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
      begin
        Standard_Complex_Laur_Systems.Copy(p.all,pp);
        if output_to_file then
          Black_Box_Root_Counting
            (outfile,integer32(nt),pp,rc,q,sols,roco,hoco,verbose-1);
          if rc /= 0 then
            if nt = 0 then
              Black_Box_Polynomial_Continuation
                (outfile,pp,q,sols,poco,verbose-1);
            else
              Black_Box_Polynomial_Continuation
                (outfile,integer32(nt),pp,q,sols,poco,verbose-1);
            end if;
          end if;
        else
          Black_Box_Root_Counting
            (integer32(nt),false,pp,rc,q,sols,roco,hoco,verbose-1);
          if rc /= 0 then
            if nt = 0 then
              Black_Box_Polynomial_Continuation(pp,q,sols,poco,verbose-1);
            else
              Black_Box_Polynomial_Continuation
                (integer32(nt),pp,q,sols,poco,verbose-1);
            end if;
          end if;
        end if;
      end;
    end if;
    tstop(timer);
    total := Elapsed_User_Time(timer);
    ended_moment := Ada.Calendar.Clock;
    if output_to_file then
      new_line(outfile);
      print_times(outfile,timer,"Solving the polynomial system");
      new_line(outfile);
      Black_Box_Helpers.Timing_Summary(outfile,roco,hoco,poco,total);
      put(outfile,"PHC ran from "); Write_Time_Stamp(outfile,start_moment);
      put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
      put_line(outfile,".");
      Write_Elapsed_Time(outfile,start_moment,ended_moment);
      Write_Number_of_Tasks(outfile,nt);
      Write_Seed_Number(outfile);
      put_line(outfile,Greeting_Banners.Version);
      Close(outfile);
      if maps /= null then
        Append(infilename,maps.all);
      else
        Black_Box_Helpers.Append_Solutions_to_Input_File
          (infilename,sols,append_sols);
      end if;
    else
      new_line;
      put_line("THE SOLUTIONS :");
      if maps /= null
       then put(maps.all);
       else put(Length_Of(sols),natural32(p'last),sols);
      end if;
      Black_Box_Helpers.Timing_Summary(standard_output,roco,hoco,poco,total);
      put("PHC ran from "); Write_Time_Stamp(standard_output,start_moment);
      put(" till "); Write_Time_Stamp(standard_output,ended_moment);
      put_line(".");
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
      Write_Number_of_Tasks(standard_output,nt);
      Write_Seed_Number(standard_output);
      put_line(Greeting_Banners.Version);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                    append_sols : in boolean; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Monomial_Maps;

    timer : Timing_Widget;
    ended_moment : Ada.Calendar.Time;
    outfile : file_type;
    output_to_file : boolean;
    fail : boolean;
    sols : Solution_List;
    maps : Link_to_Array_of_Monomial_Map_Lists;
    rc : natural32;
    roco,hoco,poco,total : duration := 0.0;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solver_cases.Square_Main for Laurent systems,");
      put_line("in double double precision ...");
    end if;
    Black_Box_Helpers.Ask_Output_File(outfile,outfilename,output_to_file);
    if output_to_file
     then put(outfile,p.all);
    end if;
    tstart(timer);
    if output_to_file
     then Black_Box_Simplex_Solver(outfile,p.all,sols,fail,verbose-1);
     else Black_Box_Simplex_Solver(p.all,sols,fail,verbose-1);
    end if;
    if fail or (Length_Of(sols) = 0) then
      if output_to_file
       then Black_Box_Binomial_Solver(outfile,p.all,maps,fail,verbose-1);
       else Black_Box_Binomial_Solver(p.all,false,maps,fail,verbose-1);
      end if;
    end if;
    if fail or (Length_Of(sols) = 0) then
      declare
        pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
      begin
        DoblDobl_Complex_Laur_Systems.Copy(p.all,pp);
        if output_to_file then
          Black_Box_Root_Counting
            (outfile,integer32(nt),pp,rc,q,sols,roco,hoco,verbose-1);
          if rc /= 0 then
            if nt = 0 then
              Black_Box_Polynomial_Continuation
                (outfile,pp,q,sols,poco,verbose-1);
            else
              Black_Box_Polynomial_Continuation
                (outfile,integer32(nt),pp,q,sols,poco,verbose-1);
            end if;
          end if;
        else
          Black_Box_Root_Counting
            (integer32(nt),false,pp,rc,q,sols,roco,hoco,verbose-1);
          if rc /= 0 then
            if nt = 0 then
              Black_Box_Polynomial_Continuation(pp,q,sols,poco,verbose-1);
            else
              Black_Box_Polynomial_Continuation
                (integer32(nt),pp,q,sols,poco,verbose-1);
            end if;
          end if;
        end if;
      end;
    end if;
    tstop(timer);
    total := Elapsed_User_Time(timer);
    ended_moment := Ada.Calendar.Clock;
    if output_to_file then
      new_line(outfile);
      print_times(outfile,timer,"Solving the polynomial system");
      new_line(outfile);
      Black_Box_Helpers.Timing_Summary(outfile,roco,hoco,poco,total);
      put(outfile,"PHC ran from "); Write_Time_Stamp(outfile,start_moment);
      put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
      put_line(outfile,".");
      Write_Elapsed_Time(outfile,start_moment,ended_moment);
      Write_Number_of_Tasks(outfile,nt);
      Write_Seed_Number(outfile);
      put_line(outfile,Greeting_Banners.Version);
      Close(outfile);
     -- currently maps are not computed in double double precision
     -- if maps /= null
     --  then Append(infilename,maps.all);
     --  else
      Black_Box_Helpers.Append_Solutions_to_Input_File
        (infilename,sols,append_sols);
     -- end if;
    else
      new_line;
      put_line("THE SOLUTIONS :");
     -- if maps /= null
     --  then put(maps.all);
     --  else 
      put(Length_Of(sols),natural32(p'last),sols);
     -- end if;
      Black_Box_Helpers.Timing_Summary(standard_output,roco,hoco,poco,total);
      put("PHC ran from "); Write_Time_Stamp(standard_output,start_moment);
      put(" till "); Write_Time_Stamp(standard_output,ended_moment);
      put_line(".");
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
      Write_Number_of_Tasks(standard_output,nt);
      Write_Seed_Number(standard_output);
      put_line(Greeting_Banners.Version);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                    append_sols : in boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Monomial_Maps;

    timer : Timing_Widget;
    ended_moment : Ada.Calendar.Time;
    outfile : file_type;
    output_to_file : boolean;
    fail : boolean;
    sols : Solution_List;
    maps : Link_to_Array_of_Monomial_Map_Lists;
    rc : natural32;
    roco,hoco,poco,total : duration := 0.0;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solver_cases.Square_Main for Laurent systems,");
      put_line("in quad double precision ...");
    end if;
    Black_Box_Helpers.Ask_Output_File(outfile,outfilename,output_to_file);
    if output_to_file
     then put(outfile,p.all);
    end if;
    tstart(timer);
    if output_to_file
     then Black_Box_Simplex_Solver(outfile,p.all,sols,fail,verbose-1);
     else Black_Box_Simplex_Solver(p.all,sols,fail,verbose-1);
    end if;
    if fail or (Length_Of(sols) = 0) then
      if output_to_file
       then Black_Box_Binomial_Solver(outfile,p.all,maps,fail,verbose-1);
       else Black_Box_Binomial_Solver(p.all,false,maps,fail,verbose-1);
      end if;
    end if;
    if fail or (Length_Of(sols) = 0) then
      declare
        pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
      begin
        QuadDobl_Complex_Laur_Systems.Copy(p.all,pp);
        if output_to_file then
          Black_Box_Root_Counting
            (outfile,integer32(nt),pp,rc,q,sols,roco,hoco,verbose-1);
          if rc /= 0 then
            if nt = 0 then
              Black_Box_Polynomial_Continuation
                (outfile,pp,q,sols,poco,verbose-1);
            else
              Black_Box_Polynomial_Continuation
                (outfile,integer32(nt),pp,q,sols,poco,verbose-1);
            end if;
          end if;
        else
          Black_Box_Root_Counting
            (integer32(nt),false,pp,rc,q,sols,roco,hoco,verbose-1);
          if rc /= 0 then
            if nt = 0 then
              Black_Box_Polynomial_Continuation(pp,q,sols,poco,verbose-1);
            else
              Black_Box_Polynomial_Continuation
                (integer32(nt),pp,q,sols,poco,verbose-1);
            end if;
          end if;
        end if;
      end;
    end if;
    tstop(timer);
    total := Elapsed_User_Time(timer);
    ended_moment := Ada.Calendar.Clock;
    if output_to_file then
      new_line(outfile);
      print_times(outfile,timer,"Solving the polynomial system");
      new_line(outfile);
      Black_Box_Helpers.Timing_Summary(outfile,roco,hoco,poco,total);
      put(outfile,"PHC ran from "); Write_Time_Stamp(outfile,start_moment);
      put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
      put_line(outfile,".");
      Write_Elapsed_Time(outfile,start_moment,ended_moment);
      Write_Number_of_Tasks(outfile,nt);
      Write_Seed_Number(outfile);
      put_line(outfile,Greeting_Banners.Version);
      Close(outfile);
     -- currently maps are not computed in quad double precision
     -- if maps /= null
     --  then Append(infilename,maps.all);
     --  else
      Black_Box_Helpers.Append_Solutions_to_Input_File
        (infilename,sols,append_sols);
     -- end if;
    else
      new_line;
      put_line("THE SOLUTIONS :");
     -- if maps /= null
     --  then put(maps.all);
     --  else 
      put(Length_Of(sols),natural32(p'last),sols);
     -- end if;
      Black_Box_Helpers.Timing_Summary(standard_output,roco,hoco,poco,total);
      put("PHC ran from "); Write_Time_Stamp(standard_output,start_moment);
      put(" till "); Write_Time_Stamp(standard_output,ended_moment);
      put_line(".");
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
      Write_Number_of_Tasks(standard_output,nt);
      Write_Seed_Number(standard_output);
      put_line(Greeting_Banners.Version);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                    deflate,append_sols : in boolean;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    timer : timing_widget;
    outfile : file_type;
    rc : natural32;
    sols,sols0 : Solution_List;
    roco,hoco,poco,total : duration := 0.0;
    fail : boolean;
    ended_moment : Ada.Calendar.Time;
    output_to_file : boolean;

  begin
    if verbose > 0 then
      put_line
        ("-> in black_box_solver_cases.Square_Main for polynomial systems,");
      put_line("in double precision ...");
    end if;
    Black_Box_Helpers.Ask_Output_File(outfile,outfilename,output_to_file);
    if output_to_file
     then put(outfile,natural32(p'last),p.all);
    end if;
    tstart(timer);
    if Black_Box_Helpers.Are_Constants_In(p.all) then
      if output_to_file
       then Black_Box_Simplex_Solver(outfile,p.all,sols,fail,verbose-1);
       else Black_Box_Simplex_Solver(p.all,sols,fail,verbose-1);
      end if;
      fail := (fail or (Length_Of(sols) = 0));
    else
      fail := true;
    end if;
    if fail then
      declare
        pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
      begin
        Standard_Complex_Poly_Systems.Copy(p.all,pp);
        if output_to_file then
          if nt >= 2 then
            Pipelined_Root_Counting
              (outfile,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco);
          else
            Black_Box_Root_Counting
              (outfile,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco,
               verbose-1);
          end if;
        else 
          if nt >= 2 then
            Pipelined_Root_Counting
              (integer32(nt),false,pp,false,rc,q,sols,sols0,roco,hoco);
          else
            Black_Box_Root_Counting
              (integer32(nt),false,pp,false,rc,q,sols,sols0,roco,hoco,
               verbose-1);
          end if;
        end if;
       -- put_line("The system after root counting : "); put(pp);
       -- put("rc = "); put(rc,1); new_line;
       -- put("Length_Of(sols) = "); put(Length_Of(sols),1); new_line;
       -- put("Length_Of(sols0) = "); put(Length_Of(sols0),1); new_line;
        if rc /= 0 then
          Standard_Scaling.Scale(pp);
          if output_to_file then
            if nt = 0 then
              Black_Box_Polynomial_Continuation
                (outfile,deflate,pp,q,sols,sols0,poco,verbose-1);
            else
              Black_Box_Polynomial_Continuation
                (outfile,deflate,integer32(nt),pp,q,sols,sols0,poco,verbose-1);
            end if;
          else
           -- put_line("in Black_Box_Solvers.Square_Main ...");
           -- put_line("Calling Black_Box_Polynomial_Continuation ...");
           -- put_line("the start system : "); put(q);
           -- put_line("the target system : "); put(pp);
            if nt = 0 then
              Black_Box_Polynomial_Continuation
                (deflate,pp,q,sols,sols0,poco,verbose-1);
            else
              Black_Box_Polynomial_Continuation
                (deflate,integer32(nt),pp,q,sols,sols0,poco,verbose-1);
            end if;
          end if;
         -- put("Length_Of(sols) = "); put(Length_Of(sols),1); new_line;
         -- put("Length_Of(sols0) = "); put(Length_Of(sols0),1); new_line;
          Push(sols0,sols);
        end if;
      end;
    end if;
    tstop(timer);
    total := Elapsed_User_Time(timer);
    ended_moment := Ada.Calendar.Clock;
    if output_to_file then
      declare
      begin
        new_line(outfile);
        print_times(outfile,timer,"Solving the polynomial system");
        new_line(outfile);
        Black_Box_Helpers.Timing_Summary(outfile,roco,hoco,poco,total);
        put(outfile,"PHC ran from "); Write_Time_Stamp(outfile,start_moment);
        put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
        put_line(outfile,".");
        Write_Elapsed_Time(outfile,start_moment,ended_moment);
        Write_Number_of_Tasks(outfile,nt);
        Write_Seed_Number(outfile);
        put_line(outfile,Greeting_Banners.Version);
        Close(outfile);
        Black_Box_Helpers.Append_Solutions_to_Input_File
          (infilename,sols,append_sols);
      exception
        when others => put_line("exception when output to file !?");
                       raise;
      end;
    else
      new_line;
      put_line("THE SOLUTIONS :");
      put(Length_Of(sols),natural32(p'last),sols);
      Black_Box_Helpers.Timing_Summary(standard_output,roco,hoco,poco,total);
      put("PHC ran from "); Write_Time_Stamp(standard_output,start_moment);
      put(" till "); Write_Time_Stamp(standard_output,ended_moment);
      put_line(".");
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
      Write_Number_of_Tasks(standard_output,nt);
      Write_Seed_Number(standard_output);
      put_line(Greeting_Banners.Version);
    end if;
  exception
    when others =>
      put_line("Exception raised in Black_Box_Solver_Cases.Square_Main");
      put("seed = "); put(Standard_Random_Numbers.Get_Seed,1);
      put("  rc = "); put(rc,1); new_line; raise;
  end Solve;

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                    append_sols : in boolean; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    timer : timing_widget;
    outfile : file_type;
    rc : natural32;
    sols,sols0 : Solution_List;
    roco,hoco,poco,total : duration := 0.0;
    fail : boolean;
    ended_moment : Ada.Calendar.Time;
    output_to_file : boolean;

  begin
    if verbose > 0 then
      put_line
        ("-> in black_box_solver_cases.Square_Main for polynomial systems,");
      put_line("in double double precision ...");
    end if;
    Black_Box_Helpers.Ask_Output_File(outfile,outfilename,output_to_file);
    if output_to_file then
      put(outfile,natural32(p'last),1); new_line(outfile);
      put(outfile,p.all);
    end if;
    tstart(timer);
    if Black_Box_Helpers.Are_Constants_In(p.all) then
      if output_to_file
       then Black_Box_Simplex_Solver(outfile,p.all,sols,fail,verbose-1);
       else Black_Box_Simplex_Solver(p.all,sols,fail,verbose-1);
      end if;
      fail := (fail or (Length_Of(sols) = 0));
    else
      fail := true;
    end if;
    if fail then
      declare
        pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
      begin
        DoblDobl_Complex_Poly_Systems.Copy(p.all,pp);
       -- put_line("The system before root counting : "); put(pp);
        if output_to_file then
          if nt >= 2 then
            Black_Box_Root_Counting
              (outfile,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco,
               verbose-1);
          else
            Black_Box_Root_Counting
              (outfile,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco,
               verbose-1);
          end if;
        else 
          if nt >= 2 then
            Pipelined_Root_Counting
              (integer32(nt),false,pp,false,rc,q,sols,sols0,roco,hoco);
          else
            Black_Box_Root_Counting
              (integer32(nt),false,pp,false,rc,q,sols,sols0,roco,hoco,
               verbose-1);
          end if;
        end if;
       -- put_line("The system after root counting : "); put(pp);
       -- put("rc = "); put(rc,1); new_line;
       -- put("Length_Of(sols) = "); put(Length_Of(sols),1); new_line;
       -- put("Length_Of(sols0) = "); put(Length_Of(sols0),1); new_line;
        if rc /= 0 then
          DoblDobl_Scaling.Scale(pp);
          if output_to_file then
            if nt = 0 then
              Black_Box_Polynomial_Continuation
                (outfile,pp,q,sols,sols0,poco,verbose-1);
            else
              Black_Box_Polynomial_Continuation
                (outfile,integer32(nt),pp,q,sols,sols0,poco,verbose-1);
            end if;
          else
           -- put_line("in Black_Box_Solvers.Square_Main ...");
           -- put_line("Calling Black_Box_Polynomial_Continuation ...");
           -- put_line("the start system : "); put(q);
           -- put_line("the target system : "); put(pp);
            if nt = 0 then
              Black_Box_Polynomial_Continuation
                (pp,q,sols,sols0,poco,verbose-1);
            else
              Black_Box_Polynomial_Continuation
                (integer32(nt),pp,q,sols,sols0,poco,verbose-1);
            end if;
          end if;
         -- put("Length_Of(sols) = "); put(Length_Of(sols),1); new_line;
         -- put("Length_Of(sols0) = "); put(Length_Of(sols0),1); new_line;
          Push(sols0,sols);
        end if;
      end;
    end if;
    tstop(timer);
    total := Elapsed_User_Time(timer);
    ended_moment := Ada.Calendar.Clock;
    if output_to_file then
      declare
      begin
        new_line(outfile);
        print_times(outfile,timer,"Solving the polynomial system");
        new_line(outfile);
        Black_Box_Helpers.Timing_Summary(outfile,roco,hoco,poco,total);
        put(outfile,"PHC ran from "); Write_Time_Stamp(outfile,start_moment);
        put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
        put_line(outfile,".");
        Write_Elapsed_Time(outfile,start_moment,ended_moment);
        Write_Number_of_Tasks(outfile,nt);
        Write_Seed_Number(outfile);
        put_line(outfile,Greeting_Banners.Version);
        Close(outfile);
        Black_Box_Helpers.Append_Solutions_to_Input_File
          (infilename,sols,append_sols);
      exception
        when others => put_line("exception when output to file !?");
                       raise;
      end;
    else
      new_line;
      put_line("THE SOLUTIONS :");
      put(Length_Of(sols),natural32(p'last),sols);
      Black_Box_Helpers.Timing_Summary(standard_output,roco,hoco,poco,total);
      put("PHC ran from "); Write_Time_Stamp(standard_output,start_moment);
      put(" till "); Write_Time_Stamp(standard_output,ended_moment);
      put_line(".");
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
      Write_Number_of_Tasks(standard_output,nt);
      Write_Seed_Number(standard_output);
      put_line(Greeting_Banners.Version);
    end if;
  exception
    when others =>
      put_line("Exception raised in Black_Box_Solvers.Square_Main");
      put("seed = "); put(Standard_Random_Numbers.Get_Seed,1);
      put("  rc = "); put(rc,1); new_line; raise;
  end Solve;

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                    append_sols : in boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    timer : timing_widget;
    outfile : file_type;
    rc : natural32;
    sols,sols0 : Solution_List;
    roco,hoco,poco,total : duration := 0.0;
    fail : boolean;
    ended_moment : Ada.Calendar.Time;
    output_to_file : boolean;

  begin
    if verbose > 0 then
      put_line
        ("-> in black_box_solver_cases.Square_Main for polynomial systems,");
      put_line("in quad double precision ...");
    end if;
    Black_Box_Helpers.Ask_Output_File(outfile,outfilename,output_to_file);
    if output_to_file then
      put(outfile,natural32(p'last),1); new_line(outfile);
      put(outfile,p.all);
    end if;
    tstart(timer);
    if Black_Box_Helpers.Are_Constants_In(p.all) then
      if output_to_file
       then Black_Box_Simplex_Solver(outfile,p.all,sols,fail,verbose-1);
       else Black_Box_Simplex_Solver(p.all,sols,fail,verbose-1);
      end if;
      fail := (fail or (Length_Of(sols) = 0));
    else
      fail := true;
    end if;
    if fail then
      declare
        pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
      begin
        QuadDobl_Complex_Poly_Systems.Copy(p.all,pp);
       -- put_line("The system before root counting : "); put(pp);
        if output_to_file then
          if nt >= 2 then
            Pipelined_Root_Counting
              (outfile,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco);
          else
            Black_Box_Root_Counting
              (outfile,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco,
               verbose-1);
          end if;
        else 
          if nt >= 2 then
            Pipelined_Root_Counting
              (integer32(nt),false,pp,false,rc,q,sols,sols0,roco,hoco);
          else
            Black_Box_Root_Counting
              (integer32(nt),false,pp,false,rc,q,sols,sols0,roco,hoco,
               verbose-1);
          end if;
        end if;
       -- put_line("The system after root counting : "); put(pp);
       -- put("rc = "); put(rc,1); new_line;
       -- put("Length_Of(sols) = "); put(Length_Of(sols),1); new_line;
       -- put("Length_Of(sols0) = "); put(Length_Of(sols0),1); new_line;
        if rc /= 0 then
          QuadDobl_Scaling.Scale(pp);
          if output_to_file then
            if nt = 0 then
              Black_Box_Polynomial_Continuation
                (outfile,pp,q,sols,sols0,poco,verbose-1);
            else
              Black_Box_Polynomial_Continuation
                (outfile,integer32(nt),pp,q,sols,sols0,poco,verbose-1);
            end if;
          else
           -- put_line("in Black_Box_Solvers.Square_Main ...");
           -- put_line("Calling Black_Box_Polynomial_Continuation ...");
           -- put_line("the start system : "); put(q);
           -- put_line("the target system : "); put(pp);
            if nt = 0 then
              Black_Box_Polynomial_Continuation
                (pp,q,sols,sols0,poco,verbose-1);
            else
              Black_Box_Polynomial_Continuation
                (integer32(nt),pp,q,sols,sols0,poco,verbose-1);
            end if;
          end if;
         -- put("Length_Of(sols) = "); put(Length_Of(sols),1); new_line;
         -- put("Length_Of(sols0) = "); put(Length_Of(sols0),1); new_line;
          Push(sols0,sols);
        end if;
      end;
    end if;
    tstop(timer);
    total := Elapsed_User_Time(timer);
    ended_moment := Ada.Calendar.Clock;
    if output_to_file then
      declare
      begin
        new_line(outfile);
        print_times(outfile,timer,"Solving the polynomial system");
        new_line(outfile);
        Black_Box_Helpers.Timing_Summary(outfile,roco,hoco,poco,total);
        put(outfile,"PHC ran from "); Write_Time_Stamp(outfile,start_moment);
        put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
        put_line(outfile,".");
        Write_Elapsed_Time(outfile,start_moment,ended_moment);
        Write_Number_of_Tasks(outfile,nt);
        Write_Seed_Number(outfile);
        put_line(outfile,Greeting_Banners.Version);
        Close(outfile);
        Black_Box_Helpers.Append_Solutions_to_Input_File
          (infilename,sols,append_sols);
      exception
        when others => put_line("exception when output to file !?");
                       raise;
      end;
    else
      new_line;
      put_line("THE SOLUTIONS :");
      put(Length_Of(sols),natural32(p'last),sols);
      Black_Box_Helpers.Timing_Summary(standard_output,roco,hoco,poco,total);
      put("PHC ran from "); Write_Time_Stamp(standard_output,start_moment);
      put(" till "); Write_Time_Stamp(standard_output,ended_moment);
      put_line(".");
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
      Write_Number_of_Tasks(standard_output,nt);
      Write_Seed_Number(standard_output);
      put_line(Greeting_Banners.Version);
    end if;
  exception
    when others =>
      put_line("Exception raised in Black_Box_Solvers.Square_Main");
      put("seed = "); put(Standard_Random_Numbers.Get_Seed,1);
      put("  rc = "); put(rc,1); new_line; raise;
  end Solve;

end Black_Box_Square_Solvers;
