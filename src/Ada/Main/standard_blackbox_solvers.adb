with Timing_Package,Time_Stamps;         use Timing_Package,Time_Stamps;
with Communications_with_User;
with File_Scanning;                      use File_Scanning;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Write_Seed_Number;
with Write_Number_of_Tasks;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with Standard_Poly_Laur_Convertors;
with Standard_System_Readers;
with Standard_Binomial_Varieties;
with Standard_Binomial_Varieties_io;
with Standard_Monomial_Maps_io;          use Standard_Monomial_Maps_io;
with Black_Box_Binomial_Solvers;         use Black_Box_Binomial_Solvers;
with Greeting_Banners;
with Black_Box_Helpers;
with Black_Box_Linear_Solvers;
with Black_Box_Single_Solvers;
with Black_Box_Square_Solvers;
with bablsolve;

package body Standard_Blackbox_Solvers is

  procedure Write_Toric_Binomial_Solutions
               ( file : in file_type; d : in natural32;
                 M : in Standard_Integer_Matrices.Matrix;
                 c : in Solution_List ) is

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
              ( nt : in natural32; start_moment : in Ada.Calendar.Time;
                p : in Laur_Sys; append_sols : in boolean;
                infilename,outfilename : in string; outfile : out file_type;
                to_file,fail : out boolean; v : in integer32 := 0 ) is

    timer : Timing_Widget;
    d : natural32;
    M : Standard_Integer_Matrices.Link_to_Matrix;
    c : Solution_List;
    ended_moment : Ada.Calendar.Time;

  begin
    if v > 0
     then put_line("-> in bablphc.Toric_Binomial_Solver 1 ...");
    end if;
    tstart(timer);
    Standard_Binomial_Varieties.Black_Box_Solver(p,fail,integer32(d),M,c);
    tstop(timer);
    fail := (d = 0) or fail;  -- binomial varieties not for isolated case!
    if not fail then
      if append_sols then
        Append_Toric_Binomial_Solutions_to_Input_File(infilename,d,M.all,c);
      end if;
      Black_Box_Helpers.Ask_Output_File(outfile,outfilename,to_file);
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
        Write_Number_of_Tasks(outfile,nt);
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
        Write_Number_of_Tasks(standard_output,nt);
        Write_Seed_Number(standard_output);
        put_line(Greeting_Banners.Version);
      end if;
    end if;
  end Toric_Binomial_Solver;

  procedure Affine_Binomial_Solver
              ( nt : in natural32; start_moment : in Ada.Calendar.Time;
                p : in Laur_Sys; append_sols : in boolean;
                infilename,outfilename : in string;
                outfile : out file_type; outnewname : out Link_to_String;
                to_file,fail : out boolean; v : in integer32 := 0 ) is

    timer : Timing_Widget;
    sols : Link_to_Array_of_Monomial_Map_Lists;
    ended_moment : Ada.Calendar.Time;

  begin
    if v > 0
     then put_line("-> in bablphc.Affine_Binomial_Solver ...");
    end if;
    Black_Box_Helpers.Ask_Output_File
      (outfile,outfilename,to_file,outnewname);
    if to_file then
      put(outfile,p'last,1); put(outfile," ");
      put(outfile,Number_of_Unknowns(p(p'first)),1); new_line(outfile);
      put(outfile,p);
    end if;
    tstart(timer);
    Black_Box_Binomial_Solver(p,false,sols,fail,v-1);
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
        Write_Number_of_Tasks(outfile,nt);
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
        Write_Number_of_Tasks(standard_output,nt);
        Write_Seed_Number(standard_output);
      end if;
    end if;
  end Affine_Binomial_Solver;

  procedure Toric_Binomial_Solver
              ( nt : in natural32; start_moment : in Ada.Calendar.Time;
                p : in Poly_Sys; append_sols : in boolean;
                infilename,outfilename : in string; outfile : out file_type;
                outnewname : out Link_to_String;
                to_file,fail : out boolean; v : in integer32 := 0 ) is

    q : constant Laur_Sys(p'range)
      := Standard_Poly_Laur_Convertors.Polynomial_to_Laurent_System(p);

  begin
    if v > 0
     then put_line("-> in bablphc.Toric_Binomial_Solver 2 ...");
    end if;
   -- Toric_Binomial_Solver(q,append_sols,fail);
    Affine_Binomial_Solver
      (nt,start_moment,q,append_sols,infilename,outfilename,outfile,
       outnewname,to_file,fail,v-1);
  end Toric_Binomial_Solver;

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in Link_to_Poly_Sys; append_sols : in boolean;
                    v : in integer32 := 0 ) is

    n : constant natural32
      := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
    to_file,fail : boolean;
    outfile : file_type;
    outnewname : Link_to_String;

  begin
    if v > 0
     then put_line("-> in standard_blackbox_solvers.Solve 1 ...");
    end if;
    if p'last = p'first then
      Black_Box_Single_Solvers.Solve
        (infilename,outfilename,p(p'first),append_sols,v-1);
    elsif p'last = integer32(n) then
      Black_Box_Linear_Solvers.Solve
        (infilename,outfilename,p,n,append_sols,fail,v-1);
      if fail then
        Black_Box_Square_Solvers.Solve
          (nt,infilename,outfilename,start_moment,p,true,append_sols,v-1);
      end if;
    else
      Toric_Binomial_Solver
        (nt,start_moment,p.all,append_sols,infilename,outfilename,
         outfile,outnewname,to_file,fail,v-1);
      if fail then
        if outnewname = null
         then bablsolve(p.all,outfilename,outfile,to_file,v-1);
         else bablsolve(p.all,outnewname.all,outfile,to_file,v-1);
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in Link_to_Laur_Sys; append_sols : in boolean;
                    v : in integer32 := 0 ) is

    to_file,fail : boolean;
    outfile : file_type;
    outnewname : Link_to_String;

  begin
    if v > 0
     then put_line("-> in standard_blackbox_solvers.Solve 2 ...");
    end if;
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(p.all) then
     -- put_line("calling Toric_Binomial_Solver ...");
      Toric_Binomial_Solver
        (nt,start_moment,p.all,append_sols,infilename,outfilename,
         outfile,to_file,fail,v-1);
    else
      Affine_Binomial_Solver
        (nt,start_moment,p.all,append_sols,infilename,outfilename,
         outfile,outnewname,to_file,fail,v-1);
    end if;
    if fail then
      Black_Box_Square_Solvers.Solve
        (nt,infilename,outfilename,start_moment,p,append_sols,v-1);
    end if;
  end Solve;

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   verbose : in integer32 := 0 ) is

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    infile : file_type;
    append_sols : boolean := false;
    p : Link_to_Poly_Sys;
    q : Link_to_Laur_Sys;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in standard_blackbox_solvers.Main ...");
    end if;
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
      Solve(nt,infilename,outfilename,start_moment,q,append_sols,verbose-1);
    else
      declare
        use Standard_Laur_Poly_Convertors;
        t : constant Poly_Sys(q'range)
          := Positive_Laurent_Polynomial_System(q.all);
      begin
        p := new Poly_Sys'(t);
        Solve(nt,infilename,outfilename,start_moment,p,append_sols,verbose-1);
      end;
    end if;
  end Main;

end Standard_Blackbox_Solvers;
