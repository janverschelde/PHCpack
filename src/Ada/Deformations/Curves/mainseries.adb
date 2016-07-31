with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;
with QuadDobl_System_and_Solutions_io;
with Standard_Dense_Series_VecVecs;
with DoblDobl_Dense_Series_VecVecs;
with QuadDobl_Dense_Series_VecVecs;
with Standard_Series_Poly_Systems;
with DoblDobl_Series_Poly_Systems;
with QuadDobl_Series_Poly_Systems;
with Series_and_Polynomials;
with Series_and_Solutions;
with Power_Series_Methods;              use Power_Series_Methods;

procedure mainseries ( precision : in character;
                       infilename,outfilename : in string ) is

  procedure Run_Newton
             ( file : in file_type; nq,idx : in integer32;
               p : in Standard_Complex_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in standard double precision.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions.

    use Standard_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    dim : constant integer32 := Head_Of(s).n - 1;
    srv : constant Standard_Dense_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : Standard_Series_Poly_Systems.Poly_Sys(p'range)
        := Series_and_Polynomials.System_to_Series_System(p,idx);
    nbrit : integer32 := 0;
    ans : character;
    verbose : boolean;
    timer : Timing_Widget;

  begin
    put("Give the number of steps in Newton's method : "); get(nbrit);
    put("Do you want extra diagnostic output in every Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
   -- if nq = dim then
   --   put_line("LU Newton will be applied.  See the output file ...");
   --   new_line;
   --   tstart(timer);
   --   Run_LU_Newton(file,nbrit,srp,srv,verbose);
   --   tstop(timer);
   -- else
      put_line("SVD Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_SVD_Newton(file,nbrit,srp,srv,verbose);
      tstop(timer);
   -- end if;
    Standard_Series_Poly_Systems.Clear(srp);
    new_line(file);
    print_times(file,timer,"power series Newton in double precision");
  end Run_Newton;

  procedure Run_Newton
             ( file : in file_type; nq,idx : in integer32;
               p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               s : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in double double precision.

  -- ON ENTRY :
  --   file    the output file for the results;
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions.

    use DoblDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    dim : constant integer32 := Head_Of(s).n - 1;
    srv : constant DoblDobl_Dense_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : DoblDobl_Series_Poly_Systems.Poly_Sys(p'range)
        := Series_and_Polynomials.System_to_Series_System(p,idx);
    nbrit : integer32 := 0;
    ans : character;
    verbose : boolean;
    timer : Timing_Widget;

  begin
    put("Give the number of steps in Newton's method : "); get(nbrit);
    put("Do you want extra diagnostic output in every Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
   -- if nq = dim then
   --   put_line("LU Newton will be applied.  See the output file ...");
   --   new_line;
   --   tstart(timer);
   --   Run_LU_Newton(file,nbrit,srp,srv,verbose);
   --   tstop(timer);
   -- else
      put_line("SVD Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_SVD_Newton(file,nbrit,srp,srv,verbose);
      tstop(timer);
   -- end if;
    DoblDobl_Series_Poly_Systems.Clear(srp);
    new_line(file);
    print_times(file,timer,"power series Newton in double double precision");
  end Run_Newton;

  procedure Run_Newton
             ( file : in file_type; nq,idx : in integer32;
               p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               s : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in quad double precision.

  -- ON ENTRY :
  --   file    the output file to write the results;
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions.

    use QuadDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    dim : constant integer32 := Head_Of(s).n - 1;
    srv : constant QuadDobl_Dense_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : QuadDobl_Series_Poly_Systems.Poly_Sys(p'range)
        := Series_and_Polynomials.System_to_Series_System(p,idx);
    nbrit : integer32 := 0;
    ans : character;
    verbose : boolean;
    timer : Timing_Widget;

  begin
    put("Give the number of steps in Newton's method : "); get(nbrit);
    put("Do you want extra diagnostic output in every Newton step ? (y/n) ");
    Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
   -- if nq = dim then
   --   put_line("LU Newton will be applied.  See the output file ...");
   --   new_line;
   --   tstart(timer);
   --   Run_LU_Newton(file,nbrit,srp,srv,verbose);
   --   tstop(timer);
   -- else
      put_line("SVD Newton will be applied.  See the output file ...");
      new_line;
      tstart(timer);
      Run_SVD_Newton(file,nbrit,srp,srv,verbose);
      tstop(timer);
   -- end if;
    QuadDobl_Series_Poly_Systems.Clear(srp);
    new_line(file);
    print_times(file,timer,"power series Newton in quad double precision");
  end Run_Newton;

  procedure Standard_Main is

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
    Run_Newton(outfile,nq,idx,lp.all,sols);
  end Standard_Main;

  procedure DoblDobl_Main is

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
    Run_Newton(outfile,nq,idx,lp.all,sols);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

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
    Run_Newton(outfile,nq,idx,lp.all,sols);
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user to select the working precision
  --   if precision /= '0' and then launches the corresponding procedure.

    prc : character;

  begin
    if precision = '0' then
      new_line;
      put_line("MENU to select the working precision :");
      put_line("  0. standard double precision;");
      put_line("  1. double double precision;");
      put_line("  2. quad double precision.");
      put("Type 0, 1, or 2 to select the working precision : ");
      Ask_Alternative(prc,"012");
      case prc is
        when '0' => Standard_Main;
        when '1' => DoblDobl_Main;
        when '2' => QuadDobl_Main;
        when others => null;
      end case;
    else
      case precision is
        when '1' =>
          put_line("The working precision is double precision");
          Standard_Main;
        when '2' =>
          put_line("The working precision is double double precision.");
          DoblDobl_Main;
        when '4' =>
          put_line("The working precision is quad double precision.");
          QuadDobl_Main;
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end Mainseries;
