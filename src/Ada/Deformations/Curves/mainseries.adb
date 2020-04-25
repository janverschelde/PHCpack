with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
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
with Series_Path_Trackers;
with Interactive_Pade_Trackers;
with Track_Path_Convolutions;

procedure mainseries ( precision : in character;
                       infilename,outfilename : in string;
                       verbose : in integer32 := 0 ) is

  procedure Run_Newton
             ( file : in file_type; echelon : in boolean;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Series_VecVecs.VecVec;
               vrb : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   terms in a power series solution to p.
  --   Newton's method is performed in standard double precision.

  -- ON ENTRY :
  --   file    to write the output to;
  --   echelon indicates whether to use echelon Newton;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       initial terms in series expansion solutions;
  --   vrb     the verbose level.

    nbrit,maxdeg : integer32 := 0;
    ans : character;
    verbose : boolean;
    timer : Timing_Widget;

  begin
    if vrb > 0
     then put_line("-> in mainseries.Run_Newton 1 ...");
    end if;
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
               s : in DoblDobl_Complex_Series_VecVecs.VecVec;
               vrb : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   terms in a power series solution to p.
  --   Newton's method is performed in double double precision.

  -- ON ENTRY :
  --   file    to write the output to;
  --   echelon indicates whether to use echelon Newton or not;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       initial terms in series expansion solutions;
  --   vrb     the verbose level.

    maxdeg,nbrit : integer32 := 0;
    ans : character;
    verbose : boolean;
    timer : Timing_Widget;

  begin
    if vrb > 0
     then put_line("-> in mainseries.Run_Newton 2 ...");
    end if;
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
               s : in QuadDobl_Complex_Series_VecVecs.VecVec;
               vrb : in integer32 := 0 ) is

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
    if vrb > 0
     then put_line("-> in mainseries.Run_Newton 3 ...");
    end if;
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
             ( file : in file_type; idx : in integer32;
               p : in Standard_Complex_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Solutions.Solution_List;
               vrb : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in standard double precision,
  --   starting at the constant term, at a solution.

  -- ON ENTRY :
  --   file    to write the output to;
  --   idx     index to the series parameter;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions;
  --   vrb     the verbose level.

    use Standard_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : constant Standard_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : Standard_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    if vrb > 0
     then put_line("-> in mainseries.Run_Newton_at_Constant 1 ...");
    end if;
    Run_Newton(file,false,srp,srv,vrb-1);
    Standard_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton_at_Constant;

  procedure Run_Newton_at_Constant
             ( file : in file_type; idx : in integer32;
               p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               s : in DoblDobl_Complex_Solutions.Solution_List;
               vrb : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in double double precision,
  --   starting at the constant term, at a solution.

  -- ON ENTRY :
  --   file    to write the output to;
  --   idx     index to the series parameter;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions;
  --   vrb     the verbose level.

    use DoblDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : constant DoblDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : DoblDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    if vrb > 0
     then put_line("-> in mainseries.Run_Newton_at_Constant 2 ...");
    end if;
    Run_Newton(file,false,srp,srv,vrb-1);
    DoblDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton_at_Constant;

  procedure Run_Newton_at_Constant
             ( file : in file_type; idx : in integer32;
               p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               s : in QuadDobl_Complex_Solutions.Solution_List;
               vrb : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in quad double precision,
  --   starting at the constant term, at a solution.

  -- ON ENTRY :
  --   file    to write the output to;
  --   idx     index to the series parameter;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions;
  --   vrb     the verbose level.

    use QuadDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    srv : constant QuadDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    if vrb > 0
     then put_line("-> in mainseries.Run_Newton_at_Constant 3 ...");
    end if;
    Run_Newton(file,false,srp,srv,vrb-1);
    QuadDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton_at_Constant;

  procedure Standard_Main_at_Constant ( vrb : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Performs a test in standard double precision,
  --   prompting the user for a system and its solutions.
  --   The vrb on input is the verbose level.

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    sols : Solution_List;

  begin
    if vrb > 0
     then put_line("-> in mainseries.Standard_Main_at_Constant ...");
    end if;
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
    Run_Newton_at_Constant(outfile,idx,lp.all,sols,vrb-1);
  end Standard_Main_at_Constant;

  procedure DoblDobl_Main_at_Constant ( vrb : in integer32 := 0 ) is

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
    if vrb > 0
     then put_line("-> in mainseries.DoblDobl_Main_at_Constant ...");
    end if;
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
    Run_Newton_at_Constant(outfile,idx,lp.all,sols,vrb-1);
  end DoblDobl_Main_at_Constant;

  procedure QuadDobl_Main_at_Constant ( vrb : in integer32 := 0 ) is

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
    if vrb > 0
     then put_line("-> in mainseries.QuadDobl_Main_at_Constant ...");
    end if;
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
    Run_Newton_at_Constant(outfile,idx,lp.all,sols,vrb-1);
  end QuadDobl_Main_at_Constant;

  procedure Standard_Main_at_Series ( vrb : in integer32 := 0 ) is

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
    if vrb > 0
     then put_line("-> in mainseries.Standard_Main_at_Series ...");
    end if;
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
      Run_Newton(outfile,true,srp,s,vrb-1);
      Standard_CSeries_Poly_Systems.Clear(srp);
    end;
  end Standard_Main_at_Series;

  procedure DoblDobl_Main_at_Series ( vrb : in integer32 := 0 ) is

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
    if vrb > 0
     then put_line("-> in mainseries.DoblDobl_Main_at_Constant ...");
    end if;
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
      Run_Newton(outfile,true,srp,s,vrb-1);
      DoblDobl_CSeries_Poly_Systems.Clear(srp);
    end;
  end DoblDobl_Main_at_Series;

  procedure QuadDobl_Main_at_Series ( vrb : in integer32 := 0 ) is

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
    if vrb > 0
     then put_line("-> in mainseries.QuadDobl_Main_at_Constant ...");
    end if;
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
      Run_Newton(outfile,true,srp,s,vrb-1);
      QuadDobl_CSeries_Poly_Systems.Clear(srp);
    end;
  end QuadDobl_Main_at_Series;

  procedure Nonzero_Precision_Main
              ( valprc : in character; vrb : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   This main procedure is called when the precision is determined,
  --   either interactively or via the value of precision.
  --   The value for the precision in valprc is not '0', it is
  --   '1' for double, '2' for double double, or '4' for quad double.
  --   The procedure can be called immediately if the precision is
  --   set at the command line.

    ans : character;

  begin 
    if vrb > 0
     then put_line("-> in mainseries.Nonzero_Precision_Main ...");
    end if;
    new_line;
    put_line("MENU for power series methods :");
    put_line("  1. apply polyhedral methods for tropisms;");
    put_line("  2. run Newton's method starting at a series or a point;");
    put_line("  3. track paths with Pade approximants as predictor;");
    put_line("  4. run a faster Newton-Fabry-Pade-Hesse path tracker.");
    put("Type 1, 2, 3, or 4 to select the method : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' =>
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
      when '2' =>
        new_line;
        put("Start Newton's method at a constant term ? (y/n) ");
        Ask_Yes_or_No(ans);
        case precision is
          when '1' =>
            put_line("The working precision is double precision");
            if ans = 'y'
             then Standard_Main_at_Constant(vrb-1);
             else Standard_Main_at_Series(vrb-1);
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
             then QuadDobl_Main_at_Constant(vrb-1);
             else QuadDobl_Main_at_Series(vrb-1);
            end if;
          when others => null;
        end case;
      when '3' =>
        new_line;
        put("Step-by-step interactive execution ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          case valprc is
            when '1' => Interactive_Pade_Trackers.Standard_Main(vrb-1);
            when '2' => Interactive_Pade_Trackers.DoblDobl_Main(vrb-1);
            when '4' => Interactive_Pade_Trackers.QuadDobl_Main(vrb-1);
            when others => null;
          end case;
        else
          case valprc is
            when '1' => Series_Path_Trackers.Standard_Main(vrb-1);
            when '2' => Series_Path_Trackers.DoblDobl_Main(vrb-1);
            when '4' => Series_Path_Trackers.QuadDobl_Main(vrb-1);
            when others => null;
          end case;
       end if;
      when '4' =>
        case valprc is
          when '1' => Track_Path_Convolutions.Standard_Main(vrb-1);
          when '2' => Track_Path_Convolutions.DoblDobl_Main(vrb-1);
          when '4' => Track_Path_Convolutions.QuadDobl_Main(vrb-1);
          when others => null;
        end case;
      when others => null;
    end case;
  end Nonzero_Precision_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user to select the working precision
  --   if precision /= '0' and then launches Nonzero_Precision_Main.

    prc : character;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in mainseries.Main ...");
    end if;
    if precision /= '0' then
      Nonzero_Precision_Main(precision,verbose-1);
    else
      new_line;
      put_line("MENU to select the working precision :");
      put_line("  0. standard double precision;");
      put_line("  1. double double precision;");
      put_line("  2. quad double precision.");
      put("Type 0, 1, or 2 to select the working precision : ");
      Ask_Alternative(prc,"012");
      case prc is
        when '0' => Nonzero_Precision_Main('1',verbose-1);
        when '1' => Nonzero_Precision_Main('2',verbose-1);
        when '2' => Nonzero_Precision_Main('4',verbose-1);
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end Mainseries;
