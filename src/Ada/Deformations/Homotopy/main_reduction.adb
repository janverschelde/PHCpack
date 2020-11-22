with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Numbers_io;                         use Numbers_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_System_Readers;
with DoblDobl_System_Readers;
with QuadDobl_System_Readers;
with Total_Degree_Start_Systems;         use Total_Degree_Start_Systems;
with Reduction_of_Polynomial_Systems;    use Reduction_of_Polynomial_Systems;
with Reduction_of_Nonsquare_Systems;     use Reduction_of_Nonsquare_Systems;

package body Main_Reduction is

  procedure Display_Info is

    i : array(1..12)of string(1..65);

  begin
    i( 1):="The goal of reduction is to rewrite the system into an equivalent";
    i( 2):="one  (i.e.:  the  same  finite  solutions) that has a lower total";
    i( 3):="degree, so  that  fewer  solution  paths  need  to  be  followed.";
    i( 4):="Sometimes  reduction  can  already detect whether a system has no";
    i( 5):="solutions or an infinite number of solutions.                    ";
    i( 6):="  We distinguish between linear  and  nonlinear  reduction.   The";
    i( 7):="first  type  performs  row-reduction on the coefficient matrix of";
    i( 8):="the system.  By nonlinear reduction, highest-degree monomials are";
    i( 9):="eliminated   by  replacing  a  polynomial  in  the  system  by  a";
    i(10):="Subtraction-polynomial.  This second type is more  powerful,  but";
    i(11):="also  more  expensive.   Bounds  have  to  be  set  to  limit the";
    i(12):="combinatorial enumeration.                                       ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Display_Info;

  procedure Display_Menu ( exit_opt : in boolean; ans : in out character ) is

  -- DESCRIPTION :
  --   Shows the general menu for reduction methods,
  --   prompting for a selection.

    m : array(0..3) of string(1..65);

  begin
    m(0):="  0 : No Reduction            : leave the menu                   ";
    m(1):="  1 : Linear Reduction        : triangulate coefficient matrix   ";
    m(2):="  2 : Sparse Linear Reduction : diagonalize coefficient matrix   ";
    m(3):="  3 : Nonlinear Reduction     : S-polynomial combinations        ";
    loop
      new_line;
      put_line("MENU for Reducing Polynomial Systems :");
      if exit_opt then
        for i in m'range loop
          put_line(m(i));
        end loop;
        put("Type 0, 1, 2, or 3 to select reduction, or i for info : ");
        Ask_Alternative(ans,"0123i");
      else
        for i in 1..m'last loop
          put_line(m(i));
        end loop;
        put("Type 1 , 2, or 3 to select reduction, or i for info : ");
        Ask_Alternative(ans,"123i");
      end if;
      if ans = 'i'
       then new_line; Display_Info;
      end if;
      exit when ans /= 'i';
    end loop;
  end Display_Menu;

  procedure Linear_Display_Menu
              ( exit_opt : in boolean; ans : in out character ) is

  -- DESCRIPTION :
  --   Shows the restricted menu for linear reduction methods,
  --   prompting for a selection.

    m : array(0..2) of string(1..65);

  begin
    m(0):="  0 : No Reduction            : leave the menu                   ";
    m(1):="  1 : Linear Reduction        : triangulate coefficient matrix   ";
    m(2):="  2 : Sparse Linear Reduction : diagonalize coefficient matrix   ";
    loop
      new_line;
      put_line("MENU for Reducing Polynomial Systems :");
      if exit_opt then
        for i in m'range loop
          put_line(m(i));
        end loop;
        put("Type 0, 1 or 2 to select reduction : ");
        Ask_Alternative(ans,"012");
      else
        for i in 1..m'last loop
          put_line(m(i));
        end loop;
        put("Type 1 or 2 to select reduction : ");
        Ask_Alternative(ans,"12");
      end if;
      if ans = 'i'
       then new_line; Display_Info;
      end if;
      exit when ans /= 'i';
    end loop;
  end Linear_Display_Menu;

  procedure Rationalize
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Shortens the exponent vectors when an unknown disappears.

    n : constant integer32 := p'length;
    to_rationalize : boolean;

    procedure rationalize ( p : in out Poly; k : in integer32 ) is

      procedure rat_term ( t : in out Term; cont : out boolean ) is

        tmp : constant Standard_Complex_Polynomials.Degrees
            := new Standard_Natural_Vectors.Vector'(1..n-1 => 0);

      begin
        for i in t.dg'range loop
          if i < k then
            tmp(i) := t.dg(i);
          elsif i > k then
            tmp(i-1) := t.dg(i);
          end if;
        end loop;
        Standard_Natural_Vectors.Clear
          (Standard_Natural_Vectors.Link_to_Vector(t.dg));
        t.dg := tmp;
        cont := true;
      end rat_term;
      procedure rat_terms is new Changing_Iterator(rat_term);

    begin
      rat_terms(p);
    end rationalize;

  begin
    for i in reverse 1..n loop
      to_rationalize := true;   -- suppose the i-th unknown may disappear
      for j in 1..n loop
        if Degree(p(j),i) > 0
         then to_rationalize := false;
        end if;
        exit when not to_rationalize;
      end loop;
      if to_rationalize then
        for j in 1..n loop
          rationalize(p(j),i);
        end loop;
      end if;
    end loop;
  end Rationalize;

  procedure Write_Diagnostics
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 diagonal,inconsistent,infinite : in boolean;
                 d : out natural32 ) is

  -- DESCRIPTION :
  --   Writes diagnostics after reduction.

    b : natural32;

  begin
    if not diagonal then
      if inconsistent then
        put_line("  Inconsistent system: no solutions.");
        put_line(file,"  Inconsistent system: no solutions.");
      elsif infinite then
        put("  Probably an infinite number");
        put_line(" of solutions.");
        put(file,"  Probably an infinite number");
        put_line(file," of solutions.");
      else
        put("  The total degree is ");
        put(file,"  The total degree is ");
        b := Total_Degree(p);
        put(b,1); put(file,b,1);
        put_line("."); put_line(file,".");
        d := b;
      end if;
    else
      put_line("  No initial terms could be eliminated.");
      put_line(file,"  No initial terms could be eliminated.");
    end if;
  end Write_Diagnostics;

  procedure Write_Diagnostics
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 diagonal,inconsistent,infinite : in boolean;
                 d : out natural32 ) is

  -- DESCRIPTION :
  --   Writes diagnostics after reduction.

    b : natural32;

  begin
    if not diagonal then
      if inconsistent then
        put_line("  Inconsistent system: no solutions.");
        put_line(file,"  Inconsistent system: no solutions.");
      elsif infinite then
        put("  Probably an infinite number");
        put_line(" of solutions.");
        put(file,"  Probably an infinite number");
        put_line(file," of solutions.");
      else
        put("  The total degree is ");
        put(file,"  The total degree is ");
        b := Total_Degree(p);
        put(b,1); put(file,b,1);
        put_line("."); put_line(file,".");
        d := b;
      end if;
    else
      put_line("  No initial terms could be eliminated.");
      put_line(file,"  No initial terms could be eliminated.");
    end if;
  end Write_Diagnostics;

  procedure Write_Diagnostics
               ( file : in file_type;
                 p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 diagonal,inconsistent,infinite : in boolean;
                 d : out natural32 ) is

  -- DESCRIPTION :
  --   Writes diagnostics after reduction.

    b : natural32;

  begin
    if not diagonal then
      if inconsistent then
        put_line("  Inconsistent system: no solutions.");
        put_line(file,"  Inconsistent system: no solutions.");
      elsif infinite then
        put("  Probably an infinite number");
        put_line(" of solutions.");
        put(file,"  Probably an infinite number");
        put_line(file," of solutions.");
      else
        put("  The total degree is ");
        put(file,"  The total degree is ");
        b := Total_Degree(p);
        put(b,1); put(file,b,1);
        put_line("."); put_line(file,".");
        d := b;
      end if;
    else
      put_line("  No initial terms could be eliminated.");
      put_line(file,"  No initial terms could be eliminated.");
    end if;
  end Write_Diagnostics;

  procedure Write_Results
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                timer : in Timing_Widget; banner : in string ) is
  begin
    new_line(file);
    put_line(file,"THE REDUCED SYSTEM :");
    put(file,natural32(p'last),p);
    new_line(file);
    print_times(file,timer,banner);
  end Write_Results;

  procedure Write_Results
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                timer : in Timing_Widget; banner : in string ) is
  begin
    new_line(file);
    put_line(file,"THE REDUCED SYSTEM :");
    put(file,natural32(p'last),p);
    new_line(file);
    print_times(file,timer,banner);
  end Write_Results;

  procedure Write_Results
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                timer : in Timing_Widget; banner : in string ) is
  begin
    new_line(file);
    put_line(file,"THE REDUCED SYSTEM :");
    put(file,natural32(p'last),p);
    new_line(file);
    print_times(file,timer,banner);
  end Write_Results;

  procedure Linear_Reduction
              ( file : in file_type;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                d : out natural32 ) is

    timer : Timing_Widget;
    diagonal,inconsistent,infinite : boolean := false;

  begin
    new_line(file);
    put_line(file,"LINEAR REDUCTION : ");
    tstart(timer);
    Reduce(p,diagonal,inconsistent,infinite);
    tstop(timer);
    Write_Diagnostics(file,p,diagonal,inconsistent,infinite,d);
    Write_Results(file,p,timer,"Linear Reduction");
  end Linear_Reduction;

  procedure Linear_Reduction
              ( file : in file_type;
                p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                d : out natural32 ) is

    timer : Timing_Widget;
    diagonal,inconsistent,infinite : boolean := false;

  begin
    new_line(file);
    put_line(file,"LINEAR REDUCTION : ");
    tstart(timer);
    Reduce(p,diagonal,inconsistent,infinite);
    tstop(timer);
    Write_Diagnostics(file,p,diagonal,inconsistent,infinite,d);
    Write_Results(file,p,timer,"Linear Reduction");
  end Linear_Reduction;

  procedure Linear_Reduction
              ( file : in file_type;
                p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                d : out natural32 ) is

    timer : Timing_Widget;
    diagonal,inconsistent,infinite : boolean := false;

  begin
    new_line(file);
    put_line(file,"LINEAR REDUCTION : ");
    tstart(timer);
    Reduce(p,diagonal,inconsistent,infinite);
    tstop(timer);
    Write_Diagnostics(file,p,diagonal,inconsistent,infinite,d);
    Write_Results(file,p,timer,"Linear Reduction");
  end Linear_Reduction;

  procedure Sparse_Linear_Reduction
              ( file : in file_type;
                p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                d : out natural32 ) is

    timer : Timing_Widget;
    diagonal : constant boolean := false;
    inconsistent,infinite : boolean := false;

  begin
    new_line(file);
    put_line(file,"SPARSE LINEAR REDUCTION : ");
    tstart(timer);
    Sparse_Reduce(p,inconsistent,infinite);
    tstop(timer);
    Write_Diagnostics(file,p,diagonal,inconsistent,infinite,d);
    Write_Results(file,p,timer,"Sparse Reduction");
  end Sparse_Linear_Reduction;

  procedure Sparse_Linear_Reduction
              ( file : in file_type;
                p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                d : out natural32 ) is

    timer : Timing_Widget;
    diagonal : constant boolean := false;
    inconsistent,infinite : boolean := false;

  begin
    new_line(file);
    put_line(file,"SPARSE LINEAR REDUCTION : ");
    tstart(timer);
    Sparse_Reduce(p,inconsistent,infinite);
    tstop(timer);
    Write_Diagnostics(file,p,diagonal,inconsistent,infinite,d);
    Write_Results(file,p,timer,"Sparse Reduction");
  end Sparse_Linear_Reduction;

  procedure Sparse_Linear_Reduction
              ( file : in file_type;
                p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                d : out natural32 ) is

    timer : Timing_Widget;
    diagonal : constant boolean := false;
    inconsistent,infinite : boolean := false;

  begin
    new_line(file);
    put_line(file,"SPARSE LINEAR REDUCTION : ");
    tstart(timer);
    Sparse_Reduce(p,inconsistent,infinite);
    tstop(timer);
    Write_Diagnostics(file,p,diagonal,inconsistent,infinite,d);
    Write_Results(file,p,timer,"Sparse Reduction");
  end Sparse_Linear_Reduction;

  procedure Nonlinear_Reduction
                ( file : in file_type;
                  p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                  d : out natural32 ) is

    use Standard_Complex_Poly_Systems;

    res : Poly_Sys(p'range);
   -- sparse : boolean;
    cnt_eq,cnt_sp,cnt_rp,max_eq,max_sp,max_rp : natural32;

    timer : Timing_Widget;
    b : natural32;

  begin
    Rationalize(p);
    Clear(res); Copy(p,res);
    new_line(file);
    put_line(file,"NONLINEAR REDUCTION :");
    put("  Give the limit on #equal degree replacements : ");
    Read_Natural(max_eq); cnt_eq := 0;
    put("  Give the limit on #computed S-polynomials : ");
    Read_Natural(max_sp); cnt_sp := 0;
    put("  Give the limit on #computed R-polynomials : ");
    Read_Natural(max_rp); cnt_rp := 0;
    put(file,"  The limit on #equal degree replacements : ");
    put(file,max_eq,1); new_line(file);
    put(file,"  The limit on #computed S-polynomials : ");
    put(file,max_sp,1); new_line(file);
    put(file,"  The limit on #computed R-polynomials : ");
    put(file,max_rp,1); new_line(file);
   -- sparse := false;
    tstart(timer);
   -- if sparse
   --  then Sparse_Reduce(p,res,cnt_eq,max_eq);
   --  else
    Reduce(p,res,cnt_eq,max_eq,cnt_sp,max_sp,cnt_rp,max_rp);
   -- end if;
    tstop(timer);
    Clear(p); Copy(res,p); Clear(res);
    b := Total_Degree(p);
    put("The total degree is "); put(b,1);
    put(file,"The total degree is "); put(file,b,1);
    put_line("."); put_line(file,".");
    d := b;
    new_line(file); 
    put_line(file,"Amount of arithmetic work");
    put(file,"   #equal replacements     : ");
    put(file,cnt_eq,4); new_line(file);
    put(file,"   #computed S-polynomials : ");
    put(file,cnt_sp,4); new_line(file);
    put(file,"   #computed R-polynomials : ");
    put(file,cnt_rp,4); new_line(file);
    new_line(file);
    put_line(file,"The reduced system :");
    put(file,p'length,p);
    new_line(file);
    print_times(file,timer,"Nonlinear Reduction");
    new_line(file);
  end Nonlinear_Reduction;

  procedure Overconstrained_Reduction
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

    ans : character;
    n : constant integer32 := p'length;
    m : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    res : Poly_Sys(1..m);

  begin
    new_line;
    put_line("MENU for Overconstrained Reduction :");
    put_line("  0. No reduction : leave the menu.");
    put_line("  1. Add a random combination of the remainder to the first.");
    put_line("  2. Reduce the first equations with the remainder.");
    put("Type 0, 1, or 2 to select reduction : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '1' => res := Random_Square(p);
      when '2' => res := Reduced_Square(p);
      when others => null;
    end case;
    if ans /= '0' then
      for i in 1..m loop
        Copy(res(i),p(i)); Clear(res(i));
      end loop;
      for i in m+1..n loop
        Clear(p(i));
      end loop;
    end if;
  end Overconstrained_Reduction;

  procedure Reduce ( file : in file_type;
                     p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                     d : out natural32; exit_option : in boolean;
                     verbose : in integer32 := 0 ) is

    n : constant natural32 := natural32(p'length);
    ans : character := '0';

  begin
    if verbose > 0
     then put_line("-> in main_reduction.Reduce 1 ...");
    end if;
    Display_Menu(exit_option,ans);
    case ans is
      when '1' => Linear_Reduction(file,p,d);
      when '2' => Sparse_Linear_Reduction(file,p,d);
      when '3' => Sparse_Linear_Reduction(file,p,d);
                  Nonlinear_Reduction(file,p,d);
      when others => null;
    end case;
    if Number_of_Unknowns(p(p'first)) < n
     then Overconstrained_Reduction(p);
    end if;
  end Reduce;

  procedure Reduce ( file : in file_type;
                     p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                     d : out natural32; exit_option : in boolean;
                     verbose : in integer32 := 0 ) is

    ans : character := '0';

  begin
    if verbose > 0
     then put_line("-> in main_reduction.Reduce 2 ...");
    end if;
    Linear_Display_Menu(exit_option,ans);
    case ans is
      when '1' => Linear_Reduction(file,p,d);
      when '2' => Sparse_Linear_Reduction(file,p,d);
      when others => null;
    end case;
  end Reduce;

  procedure Reduce ( file : in file_type;
                     p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                     d : out natural32; exit_option : in boolean;
                     verbose : in integer32 := 0 ) is

    ans : character := '0';

  begin
    if verbose > 0
     then put_line("-> in main_reduction.Reduce 3 ...");
    end if;
    Linear_Display_Menu(exit_option,ans);
    case ans is
      when '1' => Linear_Reduction(file,p,d);
      when '2' => Sparse_Linear_Reduction(file,p,d);
      when others => null;
    end case;
  end Reduce;

  procedure Standard_Main ( infilename,outfilename : in string;
                            verbose : in integer32 := 0 ) is

    use Standard_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    infile,outfile : file_type;
    d : natural32;
    ans : character;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in main_reduction.Standard_Main ...");
    end if;
    Standard_System_Readers.Read_System(infile,infilename,lp);
    if lp = null
     then new_line; get(lp);
    end if;
    Create_Output_File(outfile,outfilename);
    put(outfile,natural32(lp'last),lp.all); new_line(outfile);
    Main_Reduction.Reduce(outfile,lp.all,d,false);
    Close(outfile);
    new_line;
    put("Do you want the reduced system on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      declare
        redfile : file_type;
      begin
        put_line("Reading the name of the output file.");
        Read_Name_and_Create_File(redfile);
        put(redfile,natural32(lp'last),lp.all);
        Close(redfile);
      end;
    end if;
  end Standard_Main;

  procedure DoblDobl_Main ( infilename,outfilename : in string;
                            verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    infile,outfile : file_type;
    d : natural32;
    ans : character;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in main_reduction.DoblDobl_Main ...");
    end if;
    DoblDobl_System_Readers.Read_System(infile,infilename,lp);
    if lp = null
     then new_line; get(lp);
    end if;
    Create_Output_File(outfile,outfilename);
    put(outfile,natural32(lp'last),lp.all); new_line(outfile);
    Main_Reduction.Reduce(outfile,lp.all,d,false);
    Close(outfile);
    new_line;
    put("Do you want the reduced system on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      declare
        redfile : file_type;
      begin
        put_line("Reading the name of the output file.");
        Read_Name_and_Create_File(redfile);
        put(redfile,natural32(lp'last),lp.all);
        Close(redfile);
      end;
    end if;
  end DoblDobl_Main;

  procedure QuadDobl_Main ( infilename,outfilename : in string;
                            verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    infile,outfile : file_type;
    d : natural32;
    ans : character;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in mainred4.Main ...");
    end if;
    QuadDobl_System_Readers.Read_System(infile,infilename,lp);
    if lp = null
     then new_line; get(lp);
    end if;
    Create_Output_File(outfile,outfilename);
    put(outfile,natural32(lp'last),lp.all); new_line(outfile);
    Main_Reduction.Reduce(outfile,lp.all,d,false);
    Close(outfile);
    new_line;
    put("Do you want the reduced system on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      declare
        redfile : file_type;
      begin
        put_line("Reading the name of the output file.");
        Read_Name_and_Create_File(redfile);
        put(redfile,natural32(lp'last),lp.all);
        Close(redfile);
      end;
    end if;
  end QuadDobl_Main;

end Main_Reduction;
