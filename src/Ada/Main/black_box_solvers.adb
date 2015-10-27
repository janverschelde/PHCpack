with Timing_Package,Time_Stamps;         use Timing_Package,Time_Stamps;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Random_Numbers;
with Write_Seed_Number;
with Standard_Natural_Vectors;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Linear_Poly_Solvers;
with DoblDobl_Linear_Poly_Solvers;
with QuadDobl_Linear_Poly_Solvers;
with Standard_Scaling;
with DoblDobl_Scaling;
with QuadDobl_Scaling;
with Black_Box_Univariate_Solvers;       use Black_Box_Univariate_Solvers;
with Black_Box_Simplex_Solvers;          use Black_Box_Simplex_Solvers;
with Standard_Monomial_Maps;
with Standard_Monomial_Maps_io;          use Standard_Monomial_Maps_io;
with DoblDobl_Monomial_Maps;
with QuadDobl_Monomial_Maps;
with Black_Box_Binomial_Solvers;         use Black_Box_Binomial_Solvers;
with Black_Box_Factorization;            use Black_Box_Factorization;
with Black_Box_Root_Counters;            use Black_Box_Root_Counters;
with Standard_BlackBox_Continuations;    use Standard_BlackBox_Continuations;
with DoblDobl_Blackbox_Continuations;    use DoblDobl_Blackbox_Continuations;
with QuadDobl_Blackbox_Continuations;    use QuadDobl_Blackbox_Continuations;
with Greeting_Banners;

package body Black_Box_Solvers is

  function Is_Constant_In
              ( p : Standard_Complex_Polynomials.Poly ) return boolean is

    use Standard_Complex_Numbers;

    n : constant natural32
      := Standard_Complex_Polynomials.Number_of_Unknowns(p);
    z : constant Standard_Complex_Polynomials.Degrees
      := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    c : constant Complex_Number := Standard_Complex_Polynomials.Coeff(p,z);

  begin
    if REAL_PART(c) = 0.0 and IMAG_PART(c) = 0.0
     then return false;
     else return true;
    end if;
  end Is_Constant_In;

  function Is_Constant_In
              ( p : DoblDobl_Complex_Polynomials.Poly ) return boolean is

    use DoblDobl_Complex_Numbers;

    n : constant natural32
      := DoblDobl_Complex_Polynomials.Number_of_Unknowns(p);
    z : constant DoblDobl_Complex_Polynomials.Degrees
      := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    c : constant Complex_Number := DoblDobl_Complex_Polynomials.Coeff(p,z);
    zero : constant double_double := create(0.0);

  begin
    if REAL_PART(c) = zero and IMAG_PART(c) = zero
     then return false;
     else return true;
    end if;
  end Is_Constant_In;

  function Is_Constant_In
              ( p : QuadDobl_Complex_Polynomials.Poly ) return boolean is

    use QuadDobl_Complex_Numbers;

    n : constant natural32
      := QuadDobl_Complex_Polynomials.Number_of_Unknowns(p);
    z : constant QuadDobl_Complex_Polynomials.Degrees
      := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    c : constant Complex_Number := QuadDobl_Complex_Polynomials.Coeff(p,z);
    zero : constant quad_double := create(0.0);

  begin
    if REAL_PART(c) = zero and IMAG_PART(c) = zero
     then return false;
     else return true;
    end if;
  end Is_Constant_In;

  function Are_Constants_In
              ( p : Standard_Complex_Poly_Systems.Poly_Sys ) return boolean is
  begin
    for i in p'range loop
      if not Is_Constant_In(p(i))
       then return false;
      end if;
    end loop;
    return true;
  end Are_Constants_In;

  function Are_Constants_In
              ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys ) return boolean is
  begin
    for i in p'range loop
      if not Is_Constant_In(p(i))
       then return false;
      end if;
    end loop;
    return true;
  end Are_Constants_In;

  function Are_Constants_In
              ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys ) return boolean is
  begin
    for i in p'range loop
      if not Is_Constant_In(p(i))
       then return false;
      end if;
    end loop;
    return true;
  end Are_Constants_In;

  procedure Timing_Summary
              ( file : in file_type;
                roco,hoco,poco,total : in duration ) is

    b0 : constant string :=
     "  ---------------------------------------------------------------------";
    b1 : constant string :=
     "  |                    TIMING INFORMATION SUMMARY                     |";
    b2 : constant string :=
     "  |   root counts  |  start system  |  continuation  |   total time   |";

  begin
    put_line(file,b0);
    put_line(file,b1);
    put_line(file,b0);
    put_line(file,b2);
    put_line(file,b0);
    put(file,"  | ");
    print_hms(file,roco); put(file," | ");
    print_hms(file,hoco); put(file," | ");
    print_hms(file,poco); put(file," | ");
    print_hms(file,total); put_line(file," |");
    put_line(file,b0);
  end Timing_Summary;

  procedure Append_Solutions_to_Input_File
              ( infilename : in string;
                sols : in Standard_Complex_Solutions.Solution_list;
                append_sols : in boolean ) is

    use Standard_Complex_Solutions;

    infile : file_type;

  begin
    if not Is_Null(sols) and append_sols then
      Open_Append_File(infile,infilename);
      new_line(infile);
      put_line(infile,"THE SOLUTIONS :");
      put(infile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Close(infile);
    end if;
  end Append_Solutions_to_Input_File;

  procedure Append_Solutions_to_Input_File
              ( infilename : in string;
                sols : in DoblDobl_Complex_Solutions.Solution_list;
                append_sols : in boolean ) is

    use DoblDobl_Complex_Solutions;

    infile : file_type;

  begin
    if not Is_Null(sols) and append_sols then
      Open_Append_File(infile,infilename);
      new_line(infile);
      put_line(infile,"THE SOLUTIONS :");
      put(infile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Close(infile);
    end if;
  end Append_Solutions_to_Input_File;

  procedure Append_Solutions_to_Input_File
              ( infilename : in string;
                sols : in QuadDobl_Complex_Solutions.Solution_list;
                append_sols : in boolean ) is

    use QuadDobl_Complex_Solutions;

    infile : file_type;

  begin
    if not Is_Null(sols) and append_sols then
      Open_Append_File(infile,infilename);
      new_line(infile);
      put_line(infile,"THE SOLUTIONS :");
      put(infile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Close(infile);
    end if;
  end Append_Solutions_to_Input_File;

  procedure Ask_Output_File
              ( outfile : out file_type; outfilename : in string;
                output_to_file : out boolean ) is

    ans : character := 'y';

  begin
    if outfilename = "" then
      new_line;
      put("Do you want the output to file ? (y/n) ");
      Ask_Yes_or_No(ans);
    end if;
    if ans = 'y'
     then Create_Output_File(outfile,outfilename);
    end if;
    output_to_file := (ans = 'y');
  end Ask_Output_file;

  procedure Single_Main
              ( infilename,outfilename : in string;
                p : in Standard_Complex_Polynomials.Poly;
                append_sols : in boolean ) is

    use Standard_Complex_Solutions;

    n : constant natural32
      := Standard_Complex_Polynomials.Number_of_Unknowns(p);
    sols : Solution_List;
    outfile : file_type;
    output_to_file : boolean;

  begin
   -- put("we have one single polynomial in ");
   -- put(n,1); put_line(" variable(s)...");
    Ask_Output_File(outfile,outfilename,output_to_file);
    if n = 1 then
      if output_to_file then
        Black_Box_Durand_Kerner(outfile,p,sols);
        Append_Solutions_to_Input_file(infilename,sols,append_sols);
      else
        Black_Box_Durand_Kerner(p,sols);
        new_line;
        put_line("THE SOLUTIONS :");
        put(Length_Of(sols),1,sols);
      end if;
    elsif n > 1 then
      Standard_Black_Box_Factorization(infilename,outfile,p);
    else
      if output_to_file then
        put(outfile,"Number of unknowns = "); put(outfile,n,1);
        put_line(outfile,"...");
      else 
        put("Number of unknowns = "); put(n,1); put_line("...");
      end if;
    end if;
  end Single_Main;

  procedure Single_Main
              ( infilename,outfilename : in string;
                p : in DoblDobl_Complex_Polynomials.Poly;
                append_sols : in boolean ) is

    use DoblDobl_Complex_Solutions;

    n : constant natural32
      := DoblDobl_Complex_Polynomials.Number_of_Unknowns(p);
    sols : Solution_List;
    outfile : file_type;
    output_to_file : boolean;

  begin
   -- put("we have one single polynomial in ");
   -- put(n,1); put_line(" variable(s)...");
    Ask_Output_File(outfile,outfilename,output_to_file);
    if n = 1 then
      if output_to_file then
        Black_Box_Durand_Kerner(outfile,p,sols);
        Append_Solutions_to_Input_file(infilename,sols,append_sols);
      else
        Black_Box_Durand_Kerner(p,sols);
        new_line;
        put_line("THE SOLUTIONS :");
        put(Length_Of(sols),1,sols);
      end if;
    elsif n > 1 then
      DoblDobl_Black_Box_Factorization(infilename,outfile,p);
    else
      if output_to_file then
        put(outfile,"Number of unknowns = "); put(outfile,n,1);
        put_line(outfile,"...");
      else 
        put("Number of unknowns = "); put(n,1); put_line("...");
      end if;
    end if;
  end Single_Main;

  procedure Single_Main
              ( infilename,outfilename : in string;
                p : in QuadDobl_Complex_Polynomials.Poly;
                append_sols : in boolean ) is

    use QuadDobl_Complex_Solutions;

    n : constant natural32
      := QuadDobl_Complex_Polynomials.Number_of_Unknowns(p);
    sols : Solution_List;
    outfile : file_type;
    output_to_file : boolean;

  begin
   -- put("we have one single polynomial in ");
   -- put(n,1); put_line(" variable(s)...");
    Ask_Output_File(outfile,outfilename,output_to_file);
    if n = 1 then
      if output_to_file then
        Black_Box_Durand_Kerner(outfile,p,sols);
        Append_Solutions_to_Input_file(infilename,sols,append_sols);
      else
        Black_Box_Durand_Kerner(p,sols);
        new_line;
        put_line("THE SOLUTIONS :");
        put(Length_Of(sols),1,sols);
      end if;
    elsif n > 1 then
      QuadDobl_Black_Box_Factorization(infilename,outfile,p);
    else
      if output_to_file then
        put(outfile,"Number of unknowns = "); put(outfile,n,1);
        put_line(outfile,"...");
      else 
        put("Number of unknowns = "); put(n,1); put_line("...");
      end if;
    end if;
  end Single_Main;

  procedure Linear_Main
               ( infilename,outfilename : in string;
                 p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                 n : in natural32; append_sols : in boolean;
                 fail : out boolean ) is

    use Standard_Complex_Solutions;

    outfile : file_type;
    s : Solution(integer32(n));
    ls : Link_to_Solution;
    sols : Solution_List;
    timer : timing_widget;
    output_to_file : boolean;

  begin
    tstart(timer);
    Standard_Linear_Poly_Solvers.Solve(p.all,s,fail);
    tstop(timer);
    if not fail then
      ls := new Solution'(s);
      Construct(ls,sols);
      Ask_Output_File(outfile,outfilename,output_to_file);
      if output_to_file then
        put(outfile,natural32(p'last),p.all);
        Append_Solutions_to_Input_File(infilename,sols,append_sols);
        new_line(outfile);
        put_line(outfile,"THE SOLUTIONS :");
        put(outfile,1,n,sols);
        new_line(outfile);
        print_times(outfile,timer,"Solving the polynomial system");
        Close(outfile);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(1,natural32(s.n),sols);
      end if;
    end if;
  end Linear_Main;

  procedure Linear_Main
               ( infilename,outfilename : in string;
                 p : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 n : in natural32; append_sols : in boolean;
                 fail : out boolean ) is

    use DoblDobl_Complex_Solutions;

    outfile : file_type;
    s : Solution(integer32(n));
    ls : Link_to_Solution;
    sols : Solution_List;
    timer : timing_widget;
    output_to_file : boolean;

  begin
    tstart(timer);
    DoblDobl_Linear_Poly_Solvers.Solve(p.all,s,fail);
    tstop(timer);
    if not fail then
      ls := new Solution'(s);
      Construct(ls,sols);
      Ask_Output_File(outfile,outfilename,output_to_file);
      if output_to_file then
        put(outfile,p.all);
        Append_Solutions_to_Input_File(infilename,sols,append_sols);
        new_line(outfile);
        put_line(outfile,"THE SOLUTIONS :");
        put(outfile,1,n,sols);
        new_line(outfile);
        print_times(outfile,timer,"Solving the polynomial system");
        Close(outfile);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(1,natural32(s.n),sols);
      end if;
    end if;
  end Linear_Main;

  procedure Linear_Main
               ( infilename,outfilename : in string;
                 p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                 n : in natural32; append_sols : in boolean;
                 fail : out boolean ) is

    use QuadDobl_Complex_Solutions;

    outfile : file_type;
    s : Solution(integer32(n));
    ls : Link_to_Solution;
    sols : Solution_List;
    timer : timing_widget;
    output_to_file : boolean;

  begin
    tstart(timer);
    QuadDobl_Linear_Poly_Solvers.Solve(p.all,s,fail);
    tstop(timer);
    if not fail then
      ls := new Solution'(s);
      Construct(ls,sols);
      Ask_Output_File(outfile,outfilename,output_to_file);
      if output_to_file then
        put(outfile,p.all);
        Append_Solutions_to_Input_File(infilename,sols,append_sols);
        new_line(outfile);
        put_line(outfile,"THE SOLUTIONS :");
        put(outfile,1,n,sols);
        new_line(outfile);
        print_times(outfile,timer,"Solving the polynomial system");
        Close(outfile);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(1,natural32(s.n),sols);
      end if;
    end if;
  end Linear_Main;

  procedure Square_Main
              ( nt : in natural32; infilename,outfilename : in string;
                start_moment : in Ada.Calendar.Time;
                p : in Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                append_sols : in boolean ) is

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
    Ask_Output_File(outfile,outfilename,output_to_file);
    if output_to_file
     then put(outfile,natural32(p'last),p.all);
    end if;
    tstart(timer);
    if output_to_file
     then Black_Box_Simplex_Solver(outfile,p.all,sols,fail);
     else Black_Box_Simplex_Solver(p.all,sols,fail);
    end if;
    if fail or (Length_Of(sols) = 0) then
      if output_to_file
       then Black_Box_Binomial_Solver(outfile,p.all,maps,fail);
       else Black_Box_Binomial_Solver(p.all,false,maps,fail);
      end if;
    end if;
    if fail or (Length_Of(sols) = 0) then
      declare
        pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
      begin
        Standard_Complex_Laur_Systems.Copy(p.all,pp);
        if output_to_file then
          Black_Box_Root_Counting
            (outfile,integer32(nt),pp,rc,q,sols,roco,hoco);
          if rc /= 0 then
            if nt = 0
             then Black_Box_Polynomial_Continuation(outfile,pp,q,sols,poco);
             else Black_Box_Polynomial_Continuation
                    (outfile,integer32(nt),pp,q,sols,poco);
            end if;
          end if;
        else
          Black_Box_Root_Counting(integer32(nt),false,pp,rc,q,sols,roco,hoco);
          if rc /= 0 then
            if nt = 0
             then Black_Box_Polynomial_Continuation(pp,q,sols,poco);
             else Black_Box_Polynomial_Continuation
                    (integer32(nt),pp,q,sols,poco);
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
      Timing_Summary(outfile,roco,hoco,poco,total);
      put(outfile,"PHC ran from "); Write_Time_Stamp(outfile,start_moment);
      put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
      put_line(outfile,".");
      Write_Elapsed_Time(outfile,start_moment,ended_moment);
      Write_Seed_Number(outfile);
      put_line(outfile,Greeting_Banners.Version);
      Close(outfile);
      if maps /= null
       then Append(infilename,maps.all);
       else Append_Solutions_to_Input_File(infilename,sols,append_sols);
      end if;
    else
      new_line;
      put_line("THE SOLUTIONS :");
      if maps /= null
       then put(maps.all);
       else put(Length_Of(sols),natural32(p'last),sols);
      end if;
      Timing_Summary(standard_output,roco,hoco,poco,total);
      put("PHC ran from "); Write_Time_Stamp(standard_output,start_moment);
      put(" till "); Write_Time_Stamp(standard_output,ended_moment);
      put_line(".");
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
      Write_Seed_Number(standard_output);
      put_line(Greeting_Banners.Version);
    end if;
  end Square_Main;

  procedure Square_Main
              ( nt : in natural32; infilename,outfilename : in string;
                start_moment : in Ada.Calendar.Time;
                p : in DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                append_sols : in boolean ) is

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
    Ask_Output_File(outfile,outfilename,output_to_file);
    if output_to_file
     then put(outfile,p.all);
    end if;
    tstart(timer);
    if output_to_file
     then Black_Box_Simplex_Solver(outfile,p.all,sols,fail);
     else Black_Box_Simplex_Solver(p.all,sols,fail);
    end if;
    if fail or (Length_Of(sols) = 0) then
      if output_to_file
       then Black_Box_Binomial_Solver(outfile,p.all,maps,fail);
       else Black_Box_Binomial_Solver(p.all,false,maps,fail);
      end if;
    end if;
    if fail or (Length_Of(sols) = 0) then
      declare
        pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
      begin
        DoblDobl_Complex_Laur_Systems.Copy(p.all,pp);
        if output_to_file then
          Black_Box_Root_Counting
            (outfile,integer32(nt),pp,rc,q,sols,roco,hoco);
          if rc /= 0 then
            if nt = 0
             then Black_Box_Polynomial_Continuation(outfile,pp,q,sols,poco);
             else Black_Box_Polynomial_Continuation
                    (outfile,integer32(nt),pp,q,sols,poco);
            end if;
          end if;
        else
          Black_Box_Root_Counting(integer32(nt),false,pp,rc,q,sols,roco,hoco);
          if rc /= 0 then
            if nt = 0
             then Black_Box_Polynomial_Continuation(pp,q,sols,poco);
             else Black_Box_Polynomial_Continuation
                    (integer32(nt),pp,q,sols,poco);
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
      Timing_Summary(outfile,roco,hoco,poco,total);
      put(outfile,"PHC ran from "); Write_Time_Stamp(outfile,start_moment);
      put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
      put_line(outfile,".");
      Write_Elapsed_Time(outfile,start_moment,ended_moment);
      Write_Seed_Number(outfile);
      put_line(outfile,Greeting_Banners.Version);
      Close(outfile);
     -- currently maps are not computed in double double precision
     -- if maps /= null
     --  then Append(infilename,maps.all);
     --  else
      Append_Solutions_to_Input_File(infilename,sols,append_sols);
     -- end if;
    else
      new_line;
      put_line("THE SOLUTIONS :");
     -- if maps /= null
     --  then put(maps.all);
     --  else 
      put(Length_Of(sols),natural32(p'last),sols);
     -- end if;
      Timing_Summary(standard_output,roco,hoco,poco,total);
      put("PHC ran from "); Write_Time_Stamp(standard_output,start_moment);
      put(" till "); Write_Time_Stamp(standard_output,ended_moment);
      put_line(".");
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
      Write_Seed_Number(standard_output);
      put_line(Greeting_Banners.Version);
    end if;
  end Square_Main;

  procedure Square_Main
              ( nt : in natural32; infilename,outfilename : in string;
                start_moment : in Ada.Calendar.Time;
                p : in QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                append_sols : in boolean ) is

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
    Ask_Output_File(outfile,outfilename,output_to_file);
    if output_to_file
     then put(outfile,p.all);
    end if;
    tstart(timer);
    if output_to_file
     then Black_Box_Simplex_Solver(outfile,p.all,sols,fail);
     else Black_Box_Simplex_Solver(p.all,sols,fail);
    end if;
    if fail or (Length_Of(sols) = 0) then
      if output_to_file
       then Black_Box_Binomial_Solver(outfile,p.all,maps,fail);
       else Black_Box_Binomial_Solver(p.all,false,maps,fail);
      end if;
    end if;
    if fail or (Length_Of(sols) = 0) then
      declare
        pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
      begin
        QuadDobl_Complex_Laur_Systems.Copy(p.all,pp);
        if output_to_file then
          Black_Box_Root_Counting
            (outfile,integer32(nt),pp,rc,q,sols,roco,hoco);
          if rc /= 0 then
            if nt = 0
             then Black_Box_Polynomial_Continuation(outfile,pp,q,sols,poco);
             else Black_Box_Polynomial_Continuation
                    (outfile,integer32(nt),pp,q,sols,poco);
            end if;
          end if;
        else
          Black_Box_Root_Counting(integer32(nt),false,pp,rc,q,sols,roco,hoco);
          if rc /= 0 then
            if nt = 0
             then Black_Box_Polynomial_Continuation(pp,q,sols,poco);
             else Black_Box_Polynomial_Continuation
                    (integer32(nt),pp,q,sols,poco);
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
      Timing_Summary(outfile,roco,hoco,poco,total);
      put(outfile,"PHC ran from "); Write_Time_Stamp(outfile,start_moment);
      put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
      put_line(outfile,".");
      Write_Elapsed_Time(outfile,start_moment,ended_moment);
      Write_Seed_Number(outfile);
      put_line(outfile,Greeting_Banners.Version);
      Close(outfile);
     -- currently maps are not computed in quad double precision
     -- if maps /= null
     --  then Append(infilename,maps.all);
     --  else
      Append_Solutions_to_Input_File(infilename,sols,append_sols);
     -- end if;
    else
      new_line;
      put_line("THE SOLUTIONS :");
     -- if maps /= null
     --  then put(maps.all);
     --  else 
      put(Length_Of(sols),natural32(p'last),sols);
     -- end if;
      Timing_Summary(standard_output,roco,hoco,poco,total);
      put("PHC ran from "); Write_Time_Stamp(standard_output,start_moment);
      put(" till "); Write_Time_Stamp(standard_output,ended_moment);
      put_line(".");
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
      Write_Seed_Number(standard_output);
      put_line(Greeting_Banners.Version);
    end if;
  end Square_Main;

  procedure Square_Main
              ( nt : in natural32; infilename,outfilename : in string;
                start_moment : in Ada.Calendar.Time;
                p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                append_sols : in boolean ) is

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
   -- put_line("in Black_Box_Solvers.Square_Main for regular poly system");
    Ask_Output_File(outfile,outfilename,output_to_file);
    if output_to_file
     then put(outfile,natural32(p'last),p.all);
    end if;
    tstart(timer);
    if Are_Constants_In(p.all) then
      if output_to_file
       then Black_Box_Simplex_Solver(outfile,p.all,sols,fail);
       else Black_Box_Simplex_Solver(p.all,sols,fail);
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
       -- put_line("The system before root counting : "); put(pp);
        if output_to_file then
          Black_Box_Root_Counting
            (outfile,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco);
        else 
          Black_Box_Root_Counting
            (integer32(nt),false,pp,false,rc,q,sols,sols0,roco,hoco);
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
                (outfile,pp,q,sols,sols0,poco);
            else
              Black_Box_Polynomial_Continuation
                (outfile,integer32(nt),pp,q,sols,sols0,poco);
            end if;
          else
           -- put_line("in Black_Box_Solvers.Square_Main ...");
           -- put_line("Calling Black_Box_Polynomial_Continuation ...");
           -- put_line("the start system : "); put(q);
           -- put_line("the target system : "); put(pp);
            if nt = 0
             then Black_Box_Polynomial_Continuation(pp,q,sols,sols0,poco);
             else Black_Box_Polynomial_Continuation
                    (integer32(nt),pp,q,sols,sols0,poco);
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
        Timing_Summary(outfile,roco,hoco,poco,total);
        put(outfile,"PHC ran from "); Write_Time_Stamp(outfile,start_moment);
        put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
        put_line(outfile,".");
        Write_Elapsed_Time(outfile,start_moment,ended_moment);
        Write_Seed_Number(outfile);
        put_line(outfile,Greeting_Banners.Version);
        Close(outfile);
        Append_Solutions_to_Input_File(infilename,sols,append_sols);
      exception
        when others => put_line("exception when output to file !?");
                       raise;
      end;
    else
      new_line;
      put_line("THE SOLUTIONS :");
      put(Length_Of(sols),natural32(p'last),sols);
      Timing_Summary(standard_output,roco,hoco,poco,total);
      put("PHC ran from "); Write_Time_Stamp(standard_output,start_moment);
      put(" till "); Write_Time_Stamp(standard_output,ended_moment);
      put_line(".");
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
      Write_Seed_Number(standard_output);
      put_line(Greeting_Banners.Version);
    end if;
  exception
    when others =>
      put_line("Exception raised in Black_Box_Solvers.Square_Main");
      put("seed = "); put(Standard_Random_Numbers.Get_Seed,1);
      put("  rc = "); put(rc,1); new_line; raise;
  end Square_Main;

  procedure Square_Main
              ( nt : in natural32; infilename,outfilename : in string;
                start_moment : in Ada.Calendar.Time;
                p : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                append_sols : in boolean ) is

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
   -- put_line("in Black_Box_Solvers.Square_Main for dobldobl poly system");
    Ask_Output_File(outfile,outfilename,output_to_file);
    if output_to_file then
      put(outfile,natural32(p'last),1); new_line(outfile);
      put(outfile,p.all);
    end if;
    tstart(timer);
    if Are_Constants_In(p.all) then
      if output_to_file
       then Black_Box_Simplex_Solver(outfile,p.all,sols,fail);
       else Black_Box_Simplex_Solver(p.all,sols,fail);
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
          Black_Box_Root_Counting
            (outfile,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco);
        else 
          Black_Box_Root_Counting
            (integer32(nt),false,pp,false,rc,q,sols,sols0,roco,hoco);
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
                (outfile,pp,q,sols,sols0,poco);
            else
              Black_Box_Polynomial_Continuation
                (outfile,integer32(nt),pp,q,sols,sols0,poco);
            end if;
          else
           -- put_line("in Black_Box_Solvers.Square_Main ...");
           -- put_line("Calling Black_Box_Polynomial_Continuation ...");
           -- put_line("the start system : "); put(q);
           -- put_line("the target system : "); put(pp);
            if nt = 0
             then Black_Box_Polynomial_Continuation(pp,q,sols,sols0,poco);
             else Black_Box_Polynomial_Continuation
                    (integer32(nt),pp,q,sols,sols0,poco);
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
        Timing_Summary(outfile,roco,hoco,poco,total);
        put(outfile,"PHC ran from "); Write_Time_Stamp(outfile,start_moment);
        put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
        put_line(outfile,".");
        Write_Elapsed_Time(outfile,start_moment,ended_moment);
        Write_Seed_Number(outfile);
        put_line(outfile,Greeting_Banners.Version);
        Close(outfile);
        Append_Solutions_to_Input_File(infilename,sols,append_sols);
      exception
        when others => put_line("exception when output to file !?");
                       raise;
      end;
    else
      new_line;
      put_line("THE SOLUTIONS :");
      put(Length_Of(sols),natural32(p'last),sols);
      Timing_Summary(standard_output,roco,hoco,poco,total);
      put("PHC ran from "); Write_Time_Stamp(standard_output,start_moment);
      put(" till "); Write_Time_Stamp(standard_output,ended_moment);
      put_line(".");
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
      Write_Seed_Number(standard_output);
      put_line(Greeting_Banners.Version);
    end if;
  exception
    when others =>
      put_line("Exception raised in Black_Box_Solvers.Square_Main");
      put("seed = "); put(Standard_Random_Numbers.Get_Seed,1);
      put("  rc = "); put(rc,1); new_line; raise;
  end Square_Main;

  procedure Square_Main
              ( nt : in natural32; infilename,outfilename : in string;
                start_moment : in Ada.Calendar.Time;
                p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                append_sols : in boolean ) is

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
   -- put_line("in Black_Box_Solvers.Square_Main for regular poly system");
    Ask_Output_File(outfile,outfilename,output_to_file);
    if output_to_file then
      put(outfile,natural32(p'last),1); new_line(outfile);
      put(outfile,p.all);
    end if;
    tstart(timer);
    if Are_Constants_In(p.all) then
      if output_to_file
       then Black_Box_Simplex_Solver(outfile,p.all,sols,fail);
       else Black_Box_Simplex_Solver(p.all,sols,fail);
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
          Black_Box_Root_Counting
            (outfile,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco);
        else 
          Black_Box_Root_Counting
            (integer32(nt),false,pp,false,rc,q,sols,sols0,roco,hoco);
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
                (outfile,pp,q,sols,sols0,poco);
            else
              Black_Box_Polynomial_Continuation
                (outfile,integer32(nt),pp,q,sols,sols0,poco);
            end if;
          else
           -- put_line("in Black_Box_Solvers.Square_Main ...");
           -- put_line("Calling Black_Box_Polynomial_Continuation ...");
           -- put_line("the start system : "); put(q);
           -- put_line("the target system : "); put(pp);
            if nt = 0
             then Black_Box_Polynomial_Continuation(pp,q,sols,sols0,poco);
             else Black_Box_Polynomial_Continuation
                    (integer32(nt),pp,q,sols,sols0,poco);
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
        Timing_Summary(outfile,roco,hoco,poco,total);
        put(outfile,"PHC ran from "); Write_Time_Stamp(outfile,start_moment);
        put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
        put_line(outfile,".");
        Write_Elapsed_Time(outfile,start_moment,ended_moment);
        Write_Seed_Number(outfile);
        put_line(outfile,Greeting_Banners.Version);
        Close(outfile);
        Append_Solutions_to_Input_File(infilename,sols,append_sols);
      exception
        when others => put_line("exception when output to file !?");
                       raise;
      end;
    else
      new_line;
      put_line("THE SOLUTIONS :");
      put(Length_Of(sols),natural32(p'last),sols);
      Timing_Summary(standard_output,roco,hoco,poco,total);
      put("PHC ran from "); Write_Time_Stamp(standard_output,start_moment);
      put(" till "); Write_Time_Stamp(standard_output,ended_moment);
      put_line(".");
      Write_Elapsed_Time(standard_output,start_moment,ended_moment);
      Write_Seed_Number(standard_output);
      put_line(Greeting_Banners.Version);
    end if;
  exception
    when others =>
      put_line("Exception raised in Black_Box_Solvers.Square_Main");
      put("seed = "); put(Standard_Random_Numbers.Get_Seed,1);
      put("  rc = "); put(rc,1); new_line; raise;
  end Square_Main;

  procedure Solve ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if p'first = p'last then
      n := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(Standard_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      Standard_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          Standard_Complex_Poly_Systems.Copy(p,pp);
          Black_Box_Root_Counting(0,silent,pp,false,rc,q,sols,sols0,roco,hoco);
          if rc /= 0 then
            Standard_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation(pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if p'first = p'last then
      n := DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(DoblDobl_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      DoblDobl_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          DoblDobl_Complex_Poly_Systems.Copy(p,pp);
          Black_Box_Root_Counting(0,silent,pp,false,rc,q,sols,sols0,roco,hoco);
          if rc /= 0 then
            DoblDobl_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation(pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if p'first = p'last then
      n := QuadDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(QuadDobl_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      QuadDobl_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          QuadDobl_Complex_Poly_Systems.Copy(p,pp);
          Black_Box_Root_Counting(0,silent,pp,false,rc,q,sols,sols0,roco,hoco);
          if rc /= 0 then
            QuadDobl_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation(pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if p'first = p'last then
      n := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(Standard_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      Standard_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          Standard_Complex_Poly_Systems.Copy(p,pp);
          Black_Box_Root_Counting(file,0,pp,false,rc,q,sols,sols0,roco,hoco);
          if rc /= 0 then
            Standard_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation(file,pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if p'first = p'last then
      n := DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(DoblDobl_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      DoblDobl_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          DoblDobl_Complex_Poly_Systems.Copy(p,pp);
          Black_Box_Root_Counting(file,0,pp,false,rc,q,sols,sols0,roco,hoco);
          if rc /= 0 then
            DoblDobl_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation(file,pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if p'first = p'last then
      n := QuadDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(QuadDobl_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      QuadDobl_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          QuadDobl_Complex_Poly_Systems.Copy(p,pp);
          Black_Box_Root_Counting(file,0,pp,false,rc,q,sols,sols0,roco,hoco);
          if rc /= 0 then
            QuadDobl_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation(file,pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,silent,pp,rc,q,sols,roco,hoco);
      if rc /= 0
       then Black_Box_Polynomial_Continuation(pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,silent,pp,rc,q,sols,roco,hoco);
      if rc /= 0
       then Black_Box_Polynomial_Continuation(pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,silent,pp,rc,q,sols,roco,hoco);
      if rc /= 0
       then Black_Box_Polynomial_Continuation(pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(file,0,pp,rc,q,sols,roco,hoco);
      if rc /= 0
       then Black_Box_Polynomial_Continuation(file,pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(file,0,pp,rc,q,sols,roco,hoco);
      if rc /= 0
       then Black_Box_Polynomial_Continuation(file,pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(file,0,pp,rc,q,sols,roco,hoco);
      if rc /= 0
       then Black_Box_Polynomial_Continuation(file,pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if p'first = p'last then
      n := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(Standard_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      Standard_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          Standard_Complex_Poly_Systems.Copy(p,pp);
          Black_Box_Root_Counting
            (integer32(nt),silent,pp,false,rc,q,sols,sols0,roco,hoco);
          if rc /= 0 then
            Standard_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation
              (integer32(nt),pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if p'first = p'last then
      n := DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(DoblDobl_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      DoblDobl_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          DoblDobl_Complex_Poly_Systems.Copy(p,pp);
          Black_Box_Root_Counting
            (integer32(nt),silent,pp,false,rc,q,sols,sols0,roco,hoco);
          if rc /= 0 then
            DoblDobl_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation
              (integer32(nt),pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if p'first = p'last then
      n := QuadDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(QuadDobl_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      QuadDobl_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          QuadDobl_Complex_Poly_Systems.Copy(p,pp);
          Black_Box_Root_Counting
            (integer32(nt),silent,pp,false,rc,q,sols,sols0,roco,hoco);
          if rc /= 0 then
            QuadDobl_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation
              (integer32(nt),pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List ) is
 
    use Standard_Complex_Solutions;

    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if p'first = p'last then
      n := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(Standard_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      Standard_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          Standard_Complex_Poly_Systems.Copy(p,pp);
          Black_Box_Root_Counting
            (file,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco);
          if rc /= 0 then
            Standard_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation
              (file,integer32(nt),pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is
 
    use DoblDobl_Complex_Solutions;

    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if p'first = p'last then
      n := DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(DoblDobl_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      DoblDobl_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          DoblDobl_Complex_Poly_Systems.Copy(p,pp);
          Black_Box_Root_Counting
            (file,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco);
          if rc /= 0 then
            DoblDobl_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation
              (file,integer32(nt),pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is
 
    use QuadDobl_Complex_Solutions;

    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if p'first = p'last then
      n := QuadDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(QuadDobl_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      QuadDobl_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          QuadDobl_Complex_Poly_Systems.Copy(p,pp);
          Black_Box_Root_Counting
            (file,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco);
          if rc /= 0 then
            QuadDobl_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation
              (file,integer32(nt),pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(integer32(nt),silent,pp,rc,q,sols,roco,hoco);
      if rc /= 0 then
        Black_Box_Polynomial_Continuation(integer32(nt),pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(integer32(nt),silent,pp,rc,q,sols,roco,hoco);
      if rc /= 0 then
        Black_Box_Polynomial_Continuation(integer32(nt),pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(integer32(nt),silent,pp,rc,q,sols,roco,hoco);
      if rc /= 0 then
        Black_Box_Polynomial_Continuation(integer32(nt),pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(file,integer32(nt),pp,rc,q,sols,roco,hoco);
      if rc /= 0 then
        Black_Box_Polynomial_Continuation(file,integer32(nt),pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(file,integer32(nt),pp,rc,q,sols,roco,hoco);
      if rc /= 0 then
        Black_Box_Polynomial_Continuation(file,integer32(nt),pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(file,integer32(nt),pp,rc,q,sols,roco,hoco);
      if rc /= 0 then
        Black_Box_Polynomial_Continuation(file,integer32(nt),pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

end Black_Box_Solvers;
