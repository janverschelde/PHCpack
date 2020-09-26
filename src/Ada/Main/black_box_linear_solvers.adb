with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Linear_Poly_Solvers;
with DoblDobl_Linear_Poly_Solvers;
with QuadDobl_Linear_Poly_Solvers;
with Black_Box_Helpers;

package body Black_Box_Linear_Solvers is

  procedure Solve ( infilename,outfilename : in string;
                    p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                    n : in natural32; append_sols : in boolean;
                    fail : out boolean; verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    outfile : file_type;
    s : Solution(integer32(n));
    ls : Link_to_Solution;
    sols : Solution_List;
    timer : timing_widget;
    output_to_file : boolean;

  begin
    if verbose > 0 then
      put_line("-> in black_box_linear_solvers.Solve 1 ...");
    end if;
    tstart(timer);
    Standard_Linear_Poly_Solvers.Solve(p.all,s,fail);
    tstop(timer);
    if not fail then
      ls := new Solution'(s);
      Construct(ls,sols);
      Black_Box_Helpers.Ask_Output_File(outfile,outfilename,output_to_file);
      if output_to_file then
        put(outfile,natural32(p'last),p.all);
        Black_Box_Helpers.Append_Solutions_to_Input_File
          (infilename,sols,append_sols);
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
  end Solve;

  procedure Solve ( infilename,outfilename : in string;
                    p : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                    n : in natural32; append_sols : in boolean;
                    fail : out boolean; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    outfile : file_type;
    s : Solution(integer32(n));
    ls : Link_to_Solution;
    sols : Solution_List;
    timer : timing_widget;
    output_to_file : boolean;

  begin
    if verbose > 0 then
      put_line("-> in black_box_linear_solvers.Solve 2 ...");
    end if;
    tstart(timer);
    DoblDobl_Linear_Poly_Solvers.Solve(p.all,s,fail);
    tstop(timer);
    if not fail then
      ls := new Solution'(s);
      Construct(ls,sols);
      Black_Box_Helpers.Ask_Output_File(outfile,outfilename,output_to_file);
      if output_to_file then
        put(outfile,p.all);
        Black_Box_Helpers.Append_Solutions_to_Input_File
          (infilename,sols,append_sols);
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
  end Solve;

  procedure Solve ( infilename,outfilename : in string;
                    p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                    n : in natural32; append_sols : in boolean;
                    fail : out boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    outfile : file_type;
    s : Solution(integer32(n));
    ls : Link_to_Solution;
    sols : Solution_List;
    timer : timing_widget;
    output_to_file : boolean;

  begin
    if verbose > 0 then
      put_line("-> in black_box_linear_solvers.Solve 3 ...");
    end if;
    tstart(timer);
    QuadDobl_Linear_Poly_Solvers.Solve(p.all,s,fail);
    tstop(timer);
    if not fail then
      ls := new Solution'(s);
      Construct(ls,sols);
      Black_Box_Helpers.Ask_Output_File(outfile,outfilename,output_to_file);
      if output_to_file then
        put(outfile,p.all);
        Black_Box_Helpers.Append_Solutions_to_Input_File
          (infilename,sols,append_sols);
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
  end Solve;

end Black_Box_Linear_Solvers;
