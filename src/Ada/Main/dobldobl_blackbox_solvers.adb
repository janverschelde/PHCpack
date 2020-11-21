with text_io;                            use text_io;
with File_Scanning;                      use File_Scanning;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with DoblDobl_Laur_Poly_Convertors;
with DoblDobl_Polynomial_Convertors;     use DoblDobl_Polynomial_Convertors;
with DoblDobl_System_Readers;
with Greeting_Banners;
with Black_Box_Linear_Solvers;
with Black_Box_Single_Solvers;
with Black_Box_Square_Solvers;
with bablsolve;

package body DoblDobl_BlackBox_Solvers is

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                    append_sols : in boolean; v : in integer32 := 0 ) is

    n : constant natural32
      := DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
    outfile : file_type;
    fail : boolean;

  begin
    if v > 0
     then put_line("-> in dobldobl_blackbox_solvers.Solve 1 ...");
    end if;
    if p'last = p'first then
      Black_Box_Single_Solvers.Solve
        (infilename,outfilename,p(p'first),append_sols,v-1);
    elsif p'last = integer32(n) then
      Black_Box_Linear_Solvers.Solve
        (infilename,outfilename,p,n,append_sols,fail,v-1);
      if fail then
        Black_Box_Square_Solvers.Solve
          (nt,infilename,outfilename,start_moment,p,append_sols,v-1);
      end if;
    else
      declare
        sp : Standard_Complex_Poly_Systems.Poly_Sys(p'range)
           := DoblDobl_Complex_to_Standard_Poly_Sys(p.all);
      begin
        bablsolve(sp,outfilename,outfile,false,v-1);
        Standard_Complex_Poly_Systems.Clear(sp);
      end;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32; infilename,outfilename : in string;
                    start_moment : in Ada.Calendar.Time;
                    p : in DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                    append_sols : in boolean; v : in integer32 := 0 ) is
  begin
    if v > 0
     then put_line("-> in dobldobl_blackbox_solvers.Solve 2 ...");
    end if;
    Black_Box_Square_Solvers.Solve
      (nt,infilename,outfilename,start_moment,p,append_sols,v-1);
  end Solve;

  procedure Main ( nt : in natural32; infilename,outfilename : in string;
                   verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Laur_Systems;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    infile : file_type;
    append_sols : boolean := false;
    p : Link_to_Poly_Sys;
    q : Link_to_Laur_Sys;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in bablphc2.Main ...");
    end if;
    DoblDobl_System_Readers.Read_System(infile,infilename,q);
    if q = null then
      put_line(Greeting_Banners.welcome & ".");
      put("Running the blackbox solver");
      if nt = 0
       then put(", no tasking");
       else put(", with "); put(nt,1); put(" tasks");
      end if;
      put_line(", in double double precision.");
      new_line; get(q);
    else
      Scan_and_Skip(infile,"SOLUTIONS",append_sols);
      append_sols := not append_sols;
      close(infile);
    end if;
    if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(q.all) then
      Solve(nt,infilename,outfilename,start_moment,q,append_sols,verbose-1);
    else
      declare
        use DoblDobl_Laur_Poly_Convertors;
        t : constant Poly_Sys(q'range)
          := Positive_Laurent_Polynomial_System(q.all);
      begin
        p := new Poly_Sys'(t);
        Solve(nt,infilename,outfilename,start_moment,p,append_sols,verbose-1);
      end;
    end if;
  end Main;

end DoblDobl_BlackBox_Solvers;
