with Ada.Calendar;
with String_Splitters;                   use String_Splitters;
with Timing_Package,Time_Stamps;         use Timing_Package,Time_Stamps;
with Numbers_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Solution_Manipulators;
with DoblDobl_Solution_Manipulators;
with QuadDobl_Solution_Manipulators;
with Black_Box_Solvers;
with Square_and_Embed_Systems;           use Square_and_Embed_Systems;
with Running_Cascades;                   use Running_Cascades;

package body Embeddings_and_Cascades is

  function Lower_Dimension ( nq,nv : in natural32 ) return natural32 is

    lowdim : natural32;

  begin
    if nq >= nv
     then lowdim := 0;  -- allow no negative values for lower bound
     else lowdim := nv-nq;
    end if;
    return lowdim;
  end Lower_Dimension;

  procedure Prompt_for_Top_Dimension
              ( nq,nv : in natural32; topdim,lowdim : out natural32 ) is
  begin
    lowdim := Lower_Dimension(nq,nv);
    loop
      put("The number of equations : "); put(nq,1); new_line;
      put("The number of variables : "); put(nv,1); new_line;
      put("-> the default, largest top dimension is "); put(nv-1,1);
      put_line(" ...");
      put("Give the expected top dimension : ");
      Numbers_io.Read_Natural(topdim);
      exit when (topdim < nv) and (topdim >= lowdim); -- check bounds
      if topdim >= nv then
        put("Error: The top dimension cannot be larger than ");
        put(nv-1,1); put_line(".");
      end if;
      if topdim < lowdim then
        put("Error: The top dimension should be at least "); 
        put(lowdim,1); put_line(".");
      end if;
      put("Please enter a number between "); put(lowdim,1);
      put(" and "); put(nv-1,1); put_line(".");
      put("The default, largest top dimension is ");
      put(nv-1,1); put_line(".");
    end loop;
  end Prompt_for_Top_Dimension;

  procedure Standard_Embed_and_Cascade
              ( nt : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    timer : Timing_Widget;
    rc : natural32;
    rocos : Link_to_String;
    sols : Solution_List;
    topsoltime : duration;
    deflate : boolean;

  begin
    if verbose > 0 then
      put("-> in embeddings_and_cascades.");
      put_line("Standard_Embed_and_Cascade 1 ...");
    end if;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    deflate := (topdim = 0);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(embsys.all,deflate,rc,rocos,sols,verbose-1);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(nt,embsys.all,deflate,rc,rocos,sols,verbose-1);
      tstop(timer);
      Standard_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if rocos /= null then -- the system could be linear ...
      new_line;
      put_line("THE ROOT COUNTS :");
      new_line; put_line(rocos.all);
    end if;
    topsoltime := Elapsed_User_Time(timer);
    new_line;
    put("The CPU time for the solver : ");
    print_hms(standard_output,topsoltime); new_line;
    ended_moment := Ada.Calendar.Clock;
    new_line;
    Write_Elapsed_Time(standard_output,start_moment,ended_moment);
    if not Is_Null(sols) then
      put("Computed "); put(Length_Of(sols),1);
      put_line(" solutions at the top dimension.");
      if topdim > 0 then
        Standard_Run_Cascade
          (nt,topdim,lowdim,embsys.all,sols,filter,factor,verbose-1);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end Standard_Embed_and_Cascade;

  procedure Standard_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    timer : Timing_Widget;
    rc : natural32;
    sols : Solution_List;
    deflate : boolean;

  begin
    if verbose > 0 then
      put("-> in embeddings_and_cascades.");
      put_line("Standard_Embed_and_Cascade 2 ...");
    end if;
    new_line;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    deflate := (topdim = 0);
    put_line(file,embsys.all);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,embsys.all,deflate,rc,sols,verbose-1);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,embsys.all,deflate,rc,sols,verbose-1);
      tstop(timer);
      Standard_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    new_line(file);
    print_times(file,timer,"solving the top dimensional system");
    ended_moment := Ada.Calendar.Clock;
    new_line(file);
    Write_Elapsed_Time(file,start_moment,ended_moment);
    flush(file);
    if not Is_Null(sols) then
      if topdim > 0 then
        Standard_Run_Cascade
          (file,name,nt,topdim,lowdim,embsys.all,sols,filter,factor,verbose-1);
      else
        new_line(file);
        put_line(file,"THE SOLUTIONS :");
        put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end Standard_Embed_and_Cascade;

  procedure Standard_Embed_and_Cascade
              ( nt : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use Standard_Complex_Laurentials;
    use Standard_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    timer : Timing_Widget;
    rc : natural32;
    rocos : Link_to_String;
    sols : Solution_List;
    topsoltime : duration;

  begin
    if verbose > 0 then
      put("-> in embeddings_and_cascades.");
      put_line("Standard_Embed_and_Cascade 3 ...");
    end if;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(embsys.all,rc,rocos,sols,verbose-1);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(nt,embsys.all,rc,rocos,sols,verbose-1);
      tstop(timer);
      Standard_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if rocos /= null then -- the system could be linear ...
      new_line;
      put_line("THE ROOT COUNTS :");
      new_line; put_line(rocos.all);
    end if;
    topsoltime := Elapsed_User_Time(timer);
    new_line;
    put("The CPU time for the solver : ");
    print_hms(standard_output,topsoltime); new_line;
    ended_moment := Ada.Calendar.Clock;
    new_line;
    Write_Elapsed_Time(standard_output,start_moment,ended_moment);
    if not Is_Null(sols) then
      put("Computed "); put(Length_Of(sols),1);
      put_line(" solutions at the top dimension.");
      if topdim > 0 then
        Standard_Run_Cascade
          (nt,topdim,lowdim,embsys.all,sols,filter,factor,verbose-1);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end Standard_Embed_and_Cascade;

  procedure Standard_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use Standard_Complex_Laurentials;
    use Standard_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    timer : Timing_Widget;
    rc : natural32;
    sols : Solution_List;

  begin
    if verbose > 0 then
      put("-> in embeddings_and_cascades.");
      put_line("Standard_Embed_and_Cascade 4 ...");
    end if;
    new_line;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    put_line(file,embsys.all);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,embsys.all,rc,sols,verbose-1);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,embsys.all,rc,sols,verbose-1);
      tstop(timer);
      Standard_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    new_line(file);
    print_times(file,timer,"solving the top dimensional system");
    ended_moment := Ada.Calendar.Clock;
    new_line(file);
    Write_Elapsed_Time(file,start_moment,ended_moment);
    flush(file);
    if not Is_Null(sols) then
      if topdim > 0 then
        Standard_Run_Cascade
          (file,name,nt,topdim,lowdim,embsys.all,sols,filter,factor,verbose-1);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end Standard_Embed_and_Cascade;

  procedure DoblDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    timer : Timing_Widget;
    rc : natural32;
    rocos : Link_to_String;
    sols : Solution_List;
    topsoltime : duration;

  begin
    if verbose > 0 then
      put("-> in embeddings_and_cascades.");
      put_line("DoblDobl_Embed_and_Cascade 1 ...");
    end if;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(embsys.all,rc,rocos,sols,verbose-1);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(nt,embsys.all,rc,rocos,sols,verbose-1);
      tstop(timer);
      DoblDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if rocos /= null then -- the system could be linear
      new_line;
      put_line("THE ROOT COUNTS :");
      new_line; put_line(rocos.all);
    end if;
    topsoltime := Elapsed_User_Time(timer);
    new_line;
    put("The CPU time for the solver : ");
    print_hms(standard_output,topsoltime); new_line;
    ended_moment := Ada.Calendar.Clock;
    new_line;
    Write_Elapsed_Time(standard_output,start_moment,ended_moment);
    if not Is_Null(sols) then
      put("Computed "); put(Length_Of(sols),1);
      put_line(" solutions at the top dimension.");
      if topdim > 0 then
        DoblDobl_Run_Cascade
          (nt,topdim,lowdim,embsys.all,sols,filter,factor,verbose-1);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end DoblDobl_Embed_and_Cascade;

  procedure DoblDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    timer : Timing_Widget;
    rc : natural32;
    sols : Solution_List;

  begin
    if verbose > 0 then
      put("-> in embeddings_and_cascades.");
      put_line("DoblDobl_Embed_and_Cascade 2 ...");
    end if;
    new_line;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    put_line(file,embsys.all);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,embsys.all,rc,sols,verbose-1);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,embsys.all,rc,sols,verbose-1);
      tstop(timer);
      DoblDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    new_line(file);
    print_times(file,timer,"solving the top dimensional system");
    ended_moment := Ada.Calendar.Clock;
    new_line(file);
    Write_Elapsed_Time(file,start_moment,ended_moment);
    flush(file);
    if not Is_Null(sols) then
      if topdim > 0 then
        DoblDobl_Run_Cascade
          (file,name,nt,topdim,lowdim,embsys.all,sols,filter,factor,verbose-1);
      else
        new_line(file);
        put_line(file,"THE SOLUTIONS :");
        put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end DoblDobl_Embed_and_Cascade;

  procedure DoblDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Laurentials;
    use DoblDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    timer : Timing_Widget;
    rc : natural32;
    rocos : Link_to_String;
    sols : Solution_List;
    topsoltime : duration;

  begin
    if verbose > 0 then
      put("-> in embeddings_and_cascades.");
      put_line("DoblDobl_Embed_and_Cascade 3 ...");
    end if;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(embsys.all,rc,rocos,sols,verbose-1);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(nt,embsys.all,rc,rocos,sols,verbose-1);
      tstop(timer);
      DoblDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if rocos /= null then -- the system could be linear ...
      new_line;
      put_line("THE ROOT COUNTS :");
      new_line; put_line(rocos.all);
    end if;
    topsoltime := Elapsed_User_Time(timer);
    new_line;
    put("The CPU time for the solver : ");
    print_hms(standard_output,topsoltime); new_line;
    ended_moment := Ada.Calendar.Clock;
    new_line;
    Write_Elapsed_Time(standard_output,start_moment,ended_moment);
    if not Is_Null(sols) then
      put("Computed "); put(Length_Of(sols),1);
      put_line(" solutions at the top dimension.");
      if topdim > 0 then
        DoblDobl_Run_Cascade
          (nt,topdim,lowdim,embsys.all,sols,filter,factor,verbose-1);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end DoblDobl_Embed_and_Cascade;

  procedure DoblDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Laurentials;
    use DoblDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    timer : Timing_Widget;
    rc : natural32;
    sols : Solution_List;

  begin
    if verbose > 0 then
      put("-> in embeddings_and_cascades.");
      put_line("DoblDobl_Embed_and_Cascade 4 ...");
    end if;
    new_line;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    put_line(file,embsys.all);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,embsys.all,rc,sols,verbose-1);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,embsys.all,rc,sols,verbose-1);
      tstop(timer);
      DoblDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    new_line(file);
    print_times(file,timer,"solving the top dimensional system");
    ended_moment := Ada.Calendar.Clock;
    new_line(file);
    Write_Elapsed_Time(file,start_moment,ended_moment);
    flush(file);
    if not Is_Null(sols) then
      if topdim > 0 then
        DoblDobl_Run_Cascade
          (file,name,nt,topdim,lowdim,embsys.all,sols,filter,factor,verbose-1);
      else
        new_line(file);
        put_line(file,"THE SOLUTIONS :");
        put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end DoblDobl_Embed_and_Cascade;

  procedure QuadDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    timer : Timing_Widget;
    rc : natural32;
    rocos : Link_to_String;
    sols : Solution_List;
    topsoltime : duration;

  begin
    if verbose > 0 then
      put("-> in embeddings_and_cascades.");
      put_line("QuadDobl_Embed_and_Cascade 1 ...");
    end if;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(embsys.all,rc,rocos,sols,verbose-1);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(nt,embsys.all,rc,rocos,sols,verbose-1);
      tstop(timer);
      QuadDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if rocos /= null then -- the system could be linear ...
      new_line;
      put_line("THE ROOT COUNTS :");
      new_line; put_line(rocos.all);
    end if;
    topsoltime := Elapsed_User_Time(timer);
    new_line;
    put("The CPU time for the solver : ");
    print_hms(standard_output,topsoltime); new_line;
    ended_moment := Ada.Calendar.Clock;
    new_line;
    Write_Elapsed_Time(standard_output,start_moment,ended_moment);
    if not Is_Null(sols) then
      put("Computed "); put(Length_Of(sols),1);
      put_line(" solutions at the top dimension.");
      if topdim > 0 then
        QuadDobl_Run_Cascade
          (nt,topdim,lowdim,embsys.all,sols,filter,factor,verbose-1);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end QuadDobl_Embed_and_Cascade;

  procedure QuadDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    timer : Timing_Widget;
    rc : natural32;
    sols : Solution_List;

  begin
    if verbose > 0 then
      put("-> in embeddings_and_cascades.");
      put_line("QuadDobl_Embed_and_Cascade 2 ...");
    end if;
    new_line;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    put_line(file,embsys.all);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,embsys.all,rc,sols,verbose-1);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,embsys.all,rc,sols,verbose-1);
      tstop(timer);
      QuadDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    new_line(file);
    print_times(file,timer,"solving the top dimensional system");
    ended_moment := Ada.Calendar.Clock;
    new_line(file);
    Write_Elapsed_Time(file,start_moment,ended_moment);
    flush(file);
    if not Is_Null(sols) then
      if topdim > 0 then
        QuadDobl_Run_Cascade
          (file,name,nt,topdim,lowdim,embsys.all,sols,filter,factor,verbose-1);
      else
        new_line(file);
        put_line(file,"THE SOLUTIONS :");
        put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end QuadDobl_Embed_and_Cascade;

  procedure QuadDobl_Embed_and_Cascade
              ( nt : in natural32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    timer : Timing_Widget;
    rc : natural32;
    rocos : Link_to_String;
    sols : Solution_List;
    topsoltime : duration;

  begin
    if verbose > 0 then
      put("-> in embeddings_and_cascades.");
      put_line("QuadDobl_Embed_and_Cascade 3 ...");
    end if;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(embsys.all,rc,rocos,sols,verbose-1);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(nt,embsys.all,rc,rocos,sols,verbose-1);
      tstop(timer);
      QuadDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if rocos /= null then
      new_line;
      put_line("THE ROOT COUNTS :");
      new_line; put_line(rocos.all);
    end if;
    topsoltime := Elapsed_User_Time(timer);
    new_line;
    put("The CPU time for the solver : ");
    print_hms(standard_output,topsoltime); new_line;
    ended_moment := Ada.Calendar.Clock;
    new_line;
    Write_Elapsed_Time(standard_output,start_moment,ended_moment);
    if not Is_Null(sols) then
      put("Computed "); put(Length_Of(sols),1);
      put_line(" solutions at the top dimension.");
      if topdim > 0 then
        QuadDobl_Run_Cascade
          (nt,topdim,lowdim,embsys.all,sols,filter,factor,verbose-1);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end QuadDobl_Embed_and_Cascade;

  procedure QuadDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));
    topdim,lowdim : natural32 := 0;
    embsys : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    timer : Timing_Widget;
    rc : natural32;
    sols : Solution_List;

  begin
    if verbose > 0 then
      put("-> in embeddings_and_cascades.");
      put_line("QuadDobl_Embed_and_Cascade 4 ...");
    end if;
    new_line;
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Square_and_Embed(p,topdim,embsys);
    put_line(file,embsys.all);
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,embsys.all,rc,sols,verbose-1);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,embsys.all,rc,sols,verbose-1);
      tstop(timer);
      QuadDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    new_line(file);
    print_times(file,timer,"solving the top dimensional system");
    ended_moment := Ada.Calendar.Clock;
    new_line(file);
    Write_Elapsed_Time(file,start_moment,ended_moment);
    flush(file);
    if not Is_Null(sols) then
      if topdim > 0 then
        QuadDobl_Run_Cascade
          (file,name,nt,topdim,lowdim,embsys.all,sols,filter,factor,verbose-1);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end QuadDobl_Embed_and_Cascade;

  procedure Standard_Solve_with_Callback
              ( nt,topdim,lowdim : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    embsys : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    rc : natural32;
    sols : Standard_Complex_Solutions.Solution_List;
    deflate : boolean;

  begin
    Square_and_Embed(p,topdim,embsys);
    deflate := (topdim = 0);
    if nt = 0 then
      Black_Box_Solvers.Solve(embsys.all,true,deflate,rc,sols);
    else
      Black_Box_Solvers.Solve(nt,embsys.all,true,deflate,rc,sols);
      Standard_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if not Standard_Complex_Solutions.Is_Null(sols) then
      if topdim = 0 then
        Report_Witness_Set(embsys.all,sols,0);
      else
        Standard_Cascade_Callback
          (nt,topdim,lowdim,embsys.all,sols,filter,factor,
           pathcnt,filtcnt,idxfac,Report_Witness_Set);
      end if;
    end if;
  end Standard_Solve_with_Callback;

  procedure Standard_Solve_with_Callback
              ( nt,topdim,lowdim : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    embsys : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    rc : natural32;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      Black_Box_Solvers.Solve(embsys.all,true,rc,sols);
    else
      Black_Box_Solvers.Solve(nt,embsys.all,true,rc,sols);
      Standard_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if not Standard_Complex_Solutions.Is_Null(sols) then
      if topdim = 0 then
        Report_Witness_Set(embsys.all,sols,0);
      else
        Standard_Cascade_Callback
          (nt,topdim,lowdim,embsys.all,sols,filter,factor,
           pathcnt,filtcnt,idxfac,Report_Witness_Set);
      end if;
    end if;
  end Standard_Solve_with_Callback;

  procedure DoblDobl_Solve_with_Callback
              ( nt,topdim,lowdim : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    embsys : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    rc : natural32;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      Black_Box_Solvers.Solve(embsys.all,true,rc,sols);
    else
      Black_Box_Solvers.Solve(nt,embsys.all,true,rc,sols);
      DoblDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if not DoblDobl_Complex_Solutions.Is_Null(sols) then
      if topdim = 0 then
        Report_Witness_Set(embsys.all,sols,0);
      else
        DoblDobl_Cascade_Callback
          (nt,topdim,lowdim,embsys.all,sols,filter,factor,
           pathcnt,filtcnt,idxfac,Report_Witness_Set);
      end if;
    end if;
  end DoblDobl_Solve_with_Callback;

  procedure DoblDobl_Solve_with_Callback
              ( nt,topdim,lowdim : in natural32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    embsys : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    rc : natural32;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      Black_Box_Solvers.Solve(embsys.all,true,rc,sols);
    else
      Black_Box_Solvers.Solve(nt,embsys.all,true,rc,sols);
      DoblDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if not DoblDobl_Complex_Solutions.Is_Null(sols) then
      if topdim = 0 then
        Report_Witness_Set(embsys.all,sols,0);
      else
        DoblDobl_Cascade_Callback
          (nt,topdim,lowdim,embsys.all,sols,filter,factor,
           pathcnt,filtcnt,idxfac,Report_Witness_Set);
      end if;
    end if;
  end DoblDobl_Solve_with_Callback;

  procedure QuadDobl_Solve_with_Callback
              ( nt,topdim,lowdim : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    embsys : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    rc : natural32;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      Black_Box_Solvers.Solve(embsys.all,true,rc,sols);
    else
      Black_Box_Solvers.Solve(nt,embsys.all,true,rc,sols);
      QuadDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if not QuadDobl_Complex_Solutions.Is_Null(sols) then
      if topdim = 0 then
        Report_Witness_Set(embsys.all,sols,0);
      else
        QuadDobl_Cascade_Callback
          (nt,topdim,lowdim,embsys.all,sols,filter,factor,
           pathcnt,filtcnt,idxfac,Report_Witness_Set);
      end if;
    end if;
  end QuadDobl_Solve_with_Callback;

  procedure QuadDobl_Solve_with_Callback
              ( nt,topdim,lowdim : in natural32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    embsys : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    rc : natural32;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Square_and_Embed(p,topdim,embsys);
    if nt = 0 then
      Black_Box_Solvers.Solve(embsys.all,true,rc,sols);
    else
      Black_Box_Solvers.Solve(nt,embsys.all,true,rc,sols);
      QuadDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
    end if;
    if not QuadDobl_Complex_Solutions.Is_Null(sols) then
      if topdim = 0 then
        Report_Witness_Set(embsys.all,sols,0);
      else
        QuadDobl_Cascade_Callback
          (nt,topdim,lowdim,embsys.all,sols,filter,factor,
           pathcnt,filtcnt,idxfac,Report_Witness_Set);
      end if;
    end if;
  end QuadDobl_Solve_with_Callback;

end Embeddings_and_Cascades;
