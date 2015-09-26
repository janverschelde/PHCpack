with Ada.Calendar;
with text_io;                            use text_io;
with Timing_Package,Time_Stamps;         use Timing_Package,Time_Stamps;
with Communications_with_User;
with File_Scanning;                      use File_Scanning;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Write_Seed_Number;
with Standard_Integer_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials; 
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with Standard_Poly_Laur_Convertors;
with Standard_System_Readers;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Binomial_Varieties;
with Standard_Binomial_Varieties_io;
with Standard_Monomial_Maps;             use Standard_Monomial_Maps;
with Standard_Monomial_Maps_io;          use Standard_Monomial_Maps_io;
with Black_Box_Binomial_Solvers;         use Black_Box_Binomial_Solvers;
with Greeting_Banners;
with Black_Box_Solvers;                  use Black_Box_Solvers;
with bablsolve;

procedure bablphc ( nt : in natural32; infilename,outfilename : in string ) is

-- NOTE about the difference between Laurent and ordinary polynomials :
--   For Laurent binomial systems (the genuine ones with negative powers),
--   a stable mixed volume or an affine solution set does not make sense.
--   For ordinary binomial systems of positive dimension, we compute all
--   affine maps, or a stable mixed volume is computed.

  start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;

  procedure Write_Toric_Binomial_Solutions
               ( file : in file_type; d : in natural32;
                 M : in Standard_Integer_Matrices.Matrix;
                 c : in Solution_List ) is

  -- DESCRIPTION :
  --   Writes the solutions of the binomial system to file.

    tmp : Solution_List := c;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Standard_Binomial_Varieties_io.Write_Header(file,natural32(M'last(1)),d);
      Standard_Binomial_Varieties_io.Write_Solution(file,d,M,ls.v);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Toric_Binomial_Solutions;

  procedure Append_Toric_Binomial_Solutions_to_Input_File
              ( name : in string; d : in natural32;
                M : in Standard_Integer_Matrices.Matrix;
                c : in Solution_List ) is

    file : file_type;

  begin
    if not Is_Null(c) then
      Communications_with_User.Open_Append_File(file,name);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      Write_Toric_Binomial_Solutions(file,d,M,c);
      close(file);
    end if;
  end Append_Toric_Binomial_Solutions_to_Input_File;

  procedure Append_Affine_Binomial_Solutions_to_Input_File
              ( name : in string;
                c : in Link_to_Array_of_Monomial_Map_Lists ) is

    file : file_type;

  begin
    if c /= null then
      Communications_with_User.Open_Append_File(file,name);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,c.all);
      close(file);
    end if;
  end Append_Affine_Binomial_Solutions_to_Input_File;

  procedure Toric_Binomial_Solver
              ( p : in Laur_Sys; append_sols : in boolean;
                fail : out boolean ) is

  -- DESCRIPTION :
  --   If p is a binomial system, then it is solved.

    timer : Timing_Widget;
    d : natural32;
    M : Standard_Integer_Matrices.Link_to_Matrix;
    c : Solution_List;
    outfile : file_type;
    to_file : boolean;
    ended_moment : Ada.Calendar.Time;

  begin
    tstart(timer);
    Standard_Binomial_Varieties.Black_Box_Solver(p,fail,integer32(d),M,c);
    tstop(timer);
    fail := (d = 0) or fail;  -- binomial varieties not for isolated case!
    if not fail then
      if append_sols then
        Append_Toric_Binomial_Solutions_to_Input_File(infilename,d,M.all,c);
      end if;
      Black_Box_Solvers.Ask_Output_File(outfile,outfilename,to_file);
      ended_moment := Ada.Calendar.Clock;
      if to_file then
        Write_Toric_Binomial_Solutions(outfile,d,M.all,c);
        new_line(outfile);
        print_times(outfile,timer,"solving the binomial system");
        new_line(outfile);
        put(outfile,"PHC ran from ");
        Write_Time_Stamp(outfile,start_moment);
        put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
        put_line(outfile,".");
        Write_Elapsed_Time(outfile,start_moment,ended_moment);
        Write_Seed_Number(outfile);
        put_line(outfile,Greeting_Banners.Version);
        close(outfile);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        Write_Toric_Binomial_Solutions(standard_output,d,M.all,c);
        new_line;
        print_times(standard_output,timer,"solving the binomial system");
        new_line;
        put("PHC ran from ");
        Write_Time_Stamp(standard_output,start_moment);
        put(" till "); Write_Time_Stamp(standard_output,ended_moment);
        put_line(".");
        Write_Elapsed_Time(standard_output,start_moment,ended_moment);
        Write_Seed_Number(standard_output);
        put_line(Greeting_Banners.Version);
      end if;
    end if;
  end Toric_Binomial_Solver;

  procedure Affine_Binomial_Solver
              ( p : in Laur_Sys; append_sols : in boolean;
                fail : out boolean ) is

  -- DESCRIPTION :
  --   If p is a binomial system, then it is solved.

    timer : Timing_Widget;
    sols : Link_to_Array_of_Monomial_Map_Lists;
    outfile : file_type;
    to_file : boolean;
    ended_moment : Ada.Calendar.Time;

  begin
    Black_Box_Solvers.Ask_Output_File(outfile,outfilename,to_file);
    if to_file then
      put(outfile,p'last,1); put(outfile," ");
      put(outfile,Number_of_Unknowns(p(p'first)),1); new_line(outfile);
      put(outfile,p);
    end if;
    tstart(timer);
    Black_Box_Binomial_Solver(p,false,sols,fail);
    tstop(timer);
    if not fail and sols /= null then
      if append_sols
       then Append_Affine_Binomial_Solutions_to_Input_File(infilename,sols);
      end if;
      ended_moment := Ada.Calendar.Clock;
      if to_file then
        new_line(outfile);
        put_line(outfile,"THE SOLUTIONS :");
        put(outfile,sols.all);
        new_line(outfile);
        Show_Degrees(outfile,sols.all);
        new_line(outfile);
        print_times(outfile,timer,"solving the binomial system");
        new_line(outfile);
        put(outfile,"PHC ran from ");
        Write_Time_Stamp(outfile,start_moment);
        put(outfile," till "); Write_Time_Stamp(outfile,ended_moment);
        put_line(outfile,".");
        Write_Elapsed_Time(outfile,start_moment,ended_moment);
        Write_Seed_Number(outfile);
        close(outfile);
      else
        new_line;
        put_line("THE SOLUTIONS :");
        put(standard_output,sols.all);
        new_line;
        Show_Degrees(sols.all);
        new_line;
        print_times(standard_output,timer,"solving the binomial system");
        new_line;
        put("PHC ran from ");
        Write_Time_Stamp(standard_output,start_moment);
        put(" till "); Write_Time_Stamp(standard_output,ended_moment);
        put_line(".");
        Write_Elapsed_Time(standard_output,start_moment,ended_moment);
        Write_Seed_Number(outfile);
      end if;
    end if;
  end Affine_Binomial_Solver;

  procedure Toric_Binomial_Solver
              ( p : in Poly_Sys; append_sols : in boolean;
                fail : out boolean ) is

  -- DESCRIPTION :
  --   If p is a binomial system, then it is solved.

    q : constant Laur_Sys(p'range)
      := Standard_Poly_Laur_Convertors.Polynomial_to_Laurent_System(p);

  begin
   -- Toric_Binomial_Solver(q,append_sols,fail);
    Affine_Binomial_Solver(q,append_sols,fail);
  end Toric_Binomial_Solver;

  procedure Solve ( p : in Link_to_Poly_Sys; append_sols : in boolean ) is

  -- DESCRIPTION :
  --   Runs the blackbox solver for a polynomial system.

    n : constant natural32
      := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
    fail : boolean;

  begin
   -- put_line("in bablphc.Solve for regular polynomial system");
    if p'last = p'first then
      Single_Main(infilename,outfilename,p(p'first),append_sols);
    elsif p'last = integer32(n) then
      Linear_Main(infilename,outfilename,p,n,append_sols,fail);
      if fail
       then Square_Main(nt,infilename,outfilename,start_moment,p,append_sols);
      end if;
    else
      Toric_Binomial_Solver(p.all,append_sols,fail);
      if fail
       then bablsolve(p.all);
      end if;
    end if;
  end Solve;

  procedure Solve ( p : in Link_to_Laur_Sys; append_sols : in boolean ) is

    fail : boolean;

  begin
   -- put_line("in bablphc.Solve for Laurent polynomial system");
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(p.all) then
     -- put_line("calling Toric_Binomial_Solver ...");
      Toric_Binomial_Solver(p.all,append_sols,fail);
    else
      Affine_Binomial_Solver(p.all,append_sols,fail);
    end if;
    if fail then
      Square_Main(nt,infilename,outfilename,start_moment,p,append_sols);
    end if;
  end Solve;

  procedure Main is

    infile : file_type;
    append_sols : boolean := false;
    p : Link_to_Poly_Sys;
    q : Link_to_Laur_Sys;

  begin
    Standard_System_Readers.Read_System(infile,infilename,q);
    if q = null then
      put_line(Greeting_Banners.welcome & ".");
      put("Running the blackbox solver");
      if nt = 0
       then put(", no tasking");
       else put(", with "); put(nt,1); put(" tasks");
      end if;
      put_line(", in double precision.");
      new_line; get(q);
    else
      Scan_and_Skip(infile,"SOLUTIONS",append_sols);
      append_sols := not append_sols;
      close(infile);
    end if;
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(q.all) then
      Solve(q,append_sols);
    else
      declare
        use Standard_Laur_Poly_Convertors;
        t : constant Poly_Sys(q'range)
          := Positive_Laurent_Polynomial_System(q.all);
      begin
        p := new Poly_Sys'(t);
        Solve(p,append_sols);
      end;
    end if;
  end Main;

begin
  Main;
end bablphc;
