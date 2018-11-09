with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_System_and_Solutions_io;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_VecVecs;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_VecVecs;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_VecVecs;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;
with Complex_Series_and_Polynomials;
with Complex_Series_and_Polynomials_io;
with Series_and_Solutions;
with Power_Series_Methods;              use Power_Series_Methods;
with Regular_Newton_Puiseux;

procedure mainseries ( precision : in character;
                       infilename,outfilename : in string ) is

  procedure Run_Newton
             ( file : in file_type; echelon : in boolean;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Series_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   terms in a power series solution to p.
  --   Newton's method is performed in standard double precision.

  -- ON ENTRY :
  --   file    to write the output to;
  --   echelon indicates whether to use echelon Newton;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       initial terms in series expansion solutions.

    nbrit,maxdeg : integer32 := 0;
    ans : character;
    verbose : boolean;
    timer : Timing_Widget;

  begin
    put("Give the number of steps in Newton's method : "); get(nbrit);
    put("Give the maximal degree of the series : "); get(maxdeg);
    put("Do you want extra diagnostic output in every Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    if echelon then
      new_line;
      put_line("Echelon Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_Echelon_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    else
      new_line;
      put_line("SVD Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_SVD_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"power series Newton in double precision");
  end Run_Newton;

  procedure Run_Newton
             ( file : in file_type; echelon : in boolean;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in DoblDobl_Complex_Series_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   terms in a power series solution to p.
  --   Newton's method is performed in double double precision.

  -- ON ENTRY :
  --   file    to write the output to;
  --   echelon indicates whether to use echelon Newton or not;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       initial terms in series expansion solutions.

    maxdeg,nbrit : integer32 := 0;
    ans : character;
    verbose : boolean;
    timer : Timing_Widget;

  begin
    put("Give the number of steps in Newton's method : "); get(nbrit);
    put("Give the maximal degree of the series : "); get(maxdeg);
    put("Do you want extra diagnostic output in every Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    if echelon then
      new_line;
      put_line("Echelon Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_Echelon_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    else
      new_line;
      put_line("SVD Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_SVD_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"power series Newton in double double precision");
  end Run_Newton;

  procedure Run_Newton
             ( file : in file_type; echelon : in boolean;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in QuadDobl_Complex_Series_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   terms in a power series solution to p.
  --   Newton's method is performed in quad double precision.

  -- ON ENTRY :
  --   file    to write the output to;
  --   echelon indicates whether to use echelon Newton;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       initial terms in series expansion solutions.

    maxdeg,nbrit : integer32 := 0;
    ans : character;
    verbose : boolean;
    timer : Timing_Widget;

  begin
    put("Give the number of steps in Newton's method : "); get(nbrit);
    put("Give the maximal degree of the series : "); get(maxdeg);
    put("Do you want extra diagnostic output in every Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    if echelon then
      new_line;
      put_line("Echelon Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_Echelon_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    else
      new_line;
      put_line("SVD Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_SVD_Newton(file,maxdeg,nbrit,p,s,verbose);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"power series Newton in quad double precision");
  end Run_Newton;

  procedure Run_Newton_at_Constant
             ( file : in file_type; nq,idx : in integer32;
               p : in Standard_Complex_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in standard double precision,
  --   starting at the constant term, at a solution.

  -- ON ENTRY :
  --   file    to write the output to;
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions.

    use Standard_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : constant Standard_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : Standard_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    Run_Newton(file,false,srp,srv);
    Standard_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton_at_Constant;

  procedure Run_Newton_at_Constant
             ( file : in file_type; nq,idx : in integer32;
               p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               s : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in double double precision,
  --   starting at the constant term, at a solution.

  -- ON ENTRY :
  --   file    to write the output to;
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions.

    use DoblDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : constant DoblDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : DoblDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    Run_Newton(file,false,srp,srv);
    DoblDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton_at_Constant;

  procedure Run_Newton_at_Constant
             ( file : in file_type; nq,idx : in integer32;
               p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               s : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in quad double precision,
  --   starting at the constant term, at a solution.

  -- ON ENTRY :
  --   file    to write the output to;
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions.

    use QuadDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : constant QuadDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    Run_Newton(file,false,srp,srv);
    QuadDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton_at_Constant;

  procedure Standard_Main_at_Constant is

  -- DESCRIPTION :
  --   Performs a test in standard double precision,
  --   prompting the user for a system and its solutions.

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    sols : Solution_List;

  begin
    if infilename = "" then
      new_line;
      put_line("Reading the file name for a system and solutions ...");
      Standard_System_and_Solutions_io.get(lp,sols);
    else
      Open_Input_File(infile,infilename);
      Standard_System_and_Solutions_io.get(infile,lp,sols);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    put("Read "); put(integer32(Length_Of(sols)),1); put_line(" solutions.");
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    new_line;
    put("Give the index of the parameter : "); get(idx);
    Run_Newton_at_Constant(outfile,nq,idx,lp.all,sols);
  end Standard_Main_at_Constant;

  procedure DoblDobl_Main_at_Constant is

  -- DESCRIPTION :
  --   Performs a test in double double precision,
  --   prompting the user for a system and its solutions.

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    sols : Solution_List;

  begin
    if infilename = "" then
      new_line;
      put_line("Reading the file name for a system and solutions ...");
      DoblDobl_System_and_Solutions_io.get(lp,sols);
    else
      Open_Input_File(infile,infilename);
      DoblDobl_System_and_Solutions_io.get(infile,lp,sols);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    put("Read "); put(integer32(Length_Of(sols)),1); put_line(" solutions.");
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    new_line;
    put("Give the index of the parameter : "); get(idx);
    Run_Newton_at_Constant(outfile,nq,idx,lp.all,sols);
  end DoblDobl_Main_at_Constant;

  procedure QuadDobl_Main_at_Constant is

  -- DESCRIPTION :
  --   Performs a test in quad double precision,
  --   prompting the user for a system and its solutions.

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    sols : Solution_List;

  begin
    if infilename = "" then
      new_line;
      put_line("Reading the file name for a system and solutions ...");
      QuadDobl_System_and_Solutions_io.get(lp,sols);
    else
      Open_Input_File(infile,infilename);
      QuadDobl_System_and_Solutions_io.get(infile,lp,sols);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    put("Read "); put(integer32(Length_Of(sols)),1); put_line(" solutions.");
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    new_line;
    put("Give the index of the parameter : "); get(idx);
    Run_Newton_at_Constant(outfile,nq,idx,lp.all,sols);
  end QuadDobl_Main_at_Constant;

  procedure Standard_Main_at_Series is

  -- DESCRIPTION :
  --   Performs a test in standard double precision,
  --   prompting the user for a system and the leading terms
  --   of a series as the start for Newton's method.

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    srv : Standard_Complex_Series_Vectors.Link_to_Vector;

  begin
    if infilename = "" then
      new_line;
      put_line("Reading a polynomial system ...");
      get(lp);
    else
      Open_Input_File(infile,infilename);
      get(infile,lp);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    new_line;
    put("Give the index of the parameter : "); get(idx);
    new_line;
    put_line("Reading a series to start Newton's method at ...");
    Complex_Series_and_Polynomials_io.get(srv);
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    declare
      s : Standard_Complex_Series_VecVecs.VecVec(1..1);
      srp : Standard_CSeries_Poly_Systems.Poly_Sys(lp'range)
          := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    begin
      s(1) := srv;
      Run_Newton(outfile,true,srp,s);
    end;
  end Standard_Main_at_Series;

  procedure DoblDobl_Main_at_Series is

  -- DESCRIPTION :
  --   Performs a test in double double precision,
  --   prompting the user for a system and the leading terms
  --   of a series as the start for Newton's method.

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    srv : DoblDobl_Complex_Series_Vectors.Link_to_Vector;

  begin
    if infilename = "" then
      new_line;
      put_line("Reading a polynomial system ...");
      get(lp);
    else
      Open_Input_File(infile,infilename);
      get(infile,lp);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    new_line;
    put("Give the index of the parameter : "); get(idx);
    new_line;
    put_line("Reading a series to start Newton's method at ...");
    Complex_Series_and_Polynomials_io.get(srv);
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    declare
      s : DoblDobl_Complex_Series_VecVecs.VecVec(1..1);
      srp : DoblDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
          := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    begin
      s(1) := srv;
      Run_Newton(outfile,true,srp,s);
    end;
  end DoblDobl_Main_at_Series;

  procedure QuadDobl_Main_at_Series is

  -- DESCRIPTION :
  --   Performs a test in quad double precision,
  --   prompting the user for a system and the leading terms
  --   of a series as the start for Newton's method.

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    srv : QuadDobl_Complex_Series_Vectors.Link_to_Vector;

  begin
    if infilename = "" then
      new_line;
      put_line("Reading a polynomial system ...");
      get(lp);
    else
      Open_Input_File(infile,infilename);
      get(infile,lp);
      close(infile);
    end if;
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    new_line;
    put("Give the index of the parameter : "); get(idx);
    new_line;
    put_line("Reading a series to start Newton's method at ...");
    Complex_Series_and_Polynomials_io.get(srv);
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    declare
      s : QuadDobl_Complex_Series_VecVecs.VecVec(1..1);
      srp : QuadDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
          := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    begin
      s(1) := srv;
      Run_Newton(outfile,true,srp,s);
    end;
  end QuadDobl_Main_at_Series;

  procedure Nonzero_Precision_Main is

  -- DESCRIPTION :
  --   This main procedure must be called when precision /= '0',
  --   when the working precision was set at the command line.

    ans : character;

  begin 
    new_line;
    put("Compute tropisms with polyhedral methods ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Using as lifting the powers of the first variable,");
      put_line("assuming coefficients are sufficiently generic ...");
      case precision is
        when '1' =>
          put_line("The working precision is double precision");
          Regular_Newton_Puiseux.Standard_Main;
        when '2' =>
          put_line("The working precision is double double precision.");
          Regular_Newton_Puiseux.DoblDobl_Main;
        when '4' =>
          put_line("The working precision is quad double precision.");
          Regular_Newton_Puiseux.QuadDobl_Main;
        when others => null;
      end case;
    else
      new_line;
      put("Start Newton's method at a constant term ? (y/n) ");
      Ask_Yes_or_No(ans);
      case precision is
        when '1' =>
          put_line("The working precision is double precision");
          if ans = 'y'
           then Standard_Main_at_Constant;
           else Standard_Main_at_Series;
          end if;
        when '2' =>
          put_line("The working precision is double double precision.");
          if ans = 'y'
           then DoblDobl_Main_at_Constant;
           else DoblDobl_Main_at_Series;
          end if;
        when '4' =>
          put_line("The working precision is quad double precision.");
          if ans = 'y'
           then QuadDobl_Main_at_Constant;
           else QuadDobl_Main_at_Series;
          end if;
        when others => null;
      end case;
    end if;
  end Nonzero_Precision_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user to select the working precision
  --   if precision /= '0' and then launches the corresponding procedure.

    prc,ans : character;

  begin
    if precision /= '0' then
      Nonzero_Precision_Main;
    else
      new_line;
      put_line("MENU to select the working precision :");
      put_line("  0. standard double precision;");
      put_line("  1. double double precision;");
      put_line("  2. quad double precision.");
      put("Type 0, 1, or 2 to select the working precision : ");
      Ask_Alternative(prc,"012");
      new_line;
      put("Compute tropisms with polyhedral methods ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put_line("Using as lifting the powers of the first variable,");
        put_line("assuming coefficients are sufficiently generic ...");
        case prc is
          when '0' => Regular_Newton_Puiseux.Standard_Main;
          when '1' => Regular_Newton_Puiseux.DoblDobl_Main;
          when '2' => Regular_Newton_Puiseux.QuadDobl_Main;
          when others => null;
        end case;
      else
        new_line;
        put("Start Newton's method at constant ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          case prc is
            when '0' => Standard_Main_at_Constant;
            when '1' => DoblDobl_Main_at_Constant;
            when '2' => QuadDobl_Main_at_Constant;
            when others => null;
          end case;
        else
          case prc is
            when '0' => Standard_Main_at_Series;
            when '1' => DoblDobl_Main_at_Series;
            when '2' => QuadDobl_Main_at_Series;
            when others => null;
          end case;
        end if;
      end if;
    end if;
  end Main;

begin
  Main;
end Mainseries;
