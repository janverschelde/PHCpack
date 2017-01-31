with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Random_Numbers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with QuadDobl_Random_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Standard_PolySys_Container;
with Standard_Solutions_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_Solutions_Container;
with PHCpack_Operations;

procedure ts_ademanypaths is

-- DESCRIPTION :
--   Procedure to test the development of Newton's method and path trackers
--   using algorithmic differentiation for the evaluation.

  function Prompt_for_Verbose return integer32 is

  -- DESCRIPTION :
  --   Asks the user if extra output is wanted and returns 0 or 1.

    ans : character;

  begin
    new_line;
    put("Extra output before and after the computations ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then return 1;
     else return 0;
    end if;
  end Prompt_for_Verbose;

  procedure Standard_ADE_Newton ( verbose : in integer32 ) is

    return_of_call : integer32;

    function newton ( v : integer32 ) return integer32;
    pragma import(C, newton, "standard_adenewton");

  begin
    return_of_call := newton(verbose);
  end Standard_ADE_Newton;

  procedure DoblDobl_ADE_Newton ( verbose : in integer32 ) is

    return_of_call : integer32;

    function newton ( v : integer32 ) return integer32;
    pragma import(C, newton, "dobldobl_adenewton");

  begin
    return_of_call := newton(verbose);
  end DoblDobl_ADE_Newton;

  procedure QuadDobl_ADE_Newton ( verbose : in integer32 ) is

    return_of_call : integer32;

    function newton ( v : integer32 ) return integer32;
    pragma import(C, newton, "quaddobl_adenewton");

  begin
    return_of_call := newton(verbose);
  end QuadDobl_ADE_Newton;

  procedure Standard_ADE_Track_One
              ( verbose : in integer32; 
                gamma : in Standard_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   With the homotopy as defined by the containers,
  --   tracks one path in standard double precision.

  -- ON ENTRY :
  --   verbose  >0 if intermediate output needs to be written to screen;
  --   gamma    the gamma constant in the homotopy.

    use Standard_Complex_Numbers;

    return_of_call : integer32;
    regam : constant double_float := REAL_PART(gamma);
    imgam : constant double_float := IMAG_PART(gamma);

    function one_track ( v : integer32; r,i : double_float ) return integer32;
    pragma import(C, one_track, "standard_adeonepath");

  begin
    return_of_call := one_track(verbose,regam,imgam);
  end Standard_ADE_Track_One;

  procedure DoblDobl_ADE_Track_One
              ( verbose : in integer32; 
                gamma : in Standard_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   With the homotopy as defined by the containers,
  --   tracks one path in double double precision.

  -- ON ENTRY :
  --   verbose  >0 if intermediate output needs to be written to screen;
  --   gamma    the gamma constant in the homotopy.

    use Standard_Complex_Numbers;

    return_of_call : integer32;
    regam : constant double_float := REAL_PART(gamma);
    imgam : constant double_float := IMAG_PART(gamma);

    function one_track ( v : integer32; r,i : double_float ) return integer32;
    pragma import(C, one_track, "dobldobl_adeonepath");

  begin
    return_of_call := one_track(verbose,regam,imgam);
  end DoblDobl_ADE_Track_One;

  procedure QuadDobl_ADE_Track_One
              ( verbose : in integer32; 
                gamma : in Standard_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   With the homotopy as defined by the containers,
  --   tracks one path in quad double precision.

  -- ON ENTRY :
  --   verbose  >0 if intermediate output needs to be written to screen;
  --   gamma    the gamma constant in the homotopy.

    use Standard_Complex_Numbers;

    return_of_call : integer32;
    regam : constant double_float := REAL_PART(gamma);
    imgam : constant double_float := IMAG_PART(gamma);

    function one_track ( v : integer32; r,i : double_float ) return integer32;
    pragma import(C, one_track, "quaddobl_adeonepath");

  begin
    return_of_call := one_track(verbose,regam,imgam);
  end QuadDobl_ADE_Track_One;

  procedure Standard_ADE_Track_Many
              ( verbose : in integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   With the homotopy as defined by the containers,
  --   calls the path trackers in standard double precision.

  -- ON ENTRY :
  --   verbose  >0 if intermediate output needs to be written to screen;
  --   gamma    the gamma constant in the homotopy.

    use Standard_Complex_Numbers;

    return_of_call : integer32;
    regam : constant double_float := REAL_PART(gamma);
    imgam : constant double_float := IMAG_PART(gamma);

    function many_track ( v : integer32; r,i : double_float ) return integer32;
    pragma import(C, many_track, "standard_ademanypaths");

  begin
    return_of_call := many_track(verbose,regam,imgam);
  end Standard_ADE_Track_Many;

  procedure DoblDobl_ADE_Track_Many
              ( verbose : in integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   With the homotopy as defined by the containers,
  --   calls the path trackers in double double precision.

  -- ON ENTRY :
  --   verbose  >0 if intermediate output needs to be written to screen;
  --   gamma    the gamma constant in the homotopy.

    use Standard_Complex_Numbers;

    return_of_call : integer32;
    regam : constant double_float := REAL_PART(gamma);
    imgam : constant double_float := IMAG_PART(gamma);

    function many_track ( v : integer32; r,i : double_float ) return integer32;
    pragma import(C, many_track, "dobldobl_ademanypaths");

  begin
    return_of_call := many_track(verbose,regam,imgam);
  end DoblDobl_ADE_Track_Many;

  procedure QuadDobl_ADE_Track_Many
              ( verbose : in integer32;
                gamma : in Standard_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   With the homotopy as defined by the containers,
  --   calls the path trackers in quad double precision.

  -- ON ENTRY :
  --   verbose  >0 if intermediate output needs to be written to screen;
  --   gamma    the gamma constant in the homotopy.

    use Standard_Complex_Numbers;

    return_of_call : integer32;
    regam : constant double_float := REAL_PART(gamma);
    imgam : constant double_float := IMAG_PART(gamma);

    function many_track ( v : integer32; r,i : double_float ) return integer32;
    pragma import(C, many_track, "quaddobl_ademanypaths");

  begin
    return_of_call := many_track(verbose,regam,imgam);
  end QuadDobl_ADE_Track_Many;

  procedure Standard_Newton is

  -- DESCRIPTION :
  --   Reads a polynomial system and a corresponding list of solutions
  --   in standard double precision.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;

  begin
    new_line;
    put_line("Reading a system with solutions ...");
    Standard_System_and_Solutions_io.get(p,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    put_line("Initializing the systems container ...");
    Standard_PolySys_Container.Initialize(p.all);
    put_line("Initializing the solutions container ...");
    Standard_Solutions_Container.Initialize(sols);
    verbose := Prompt_for_Verbose;
    Standard_ADE_Newton(verbose);
    newtsols := Standard_Solutions_Container.Retrieve;
    put_line("The solutions after Newton's method :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
  end Standard_Newton;

  procedure DoblDobl_Newton is

  -- DESCRIPTION :
  --   Reads a polynomial system and a corresponding list of solutions
  --   in double double precision.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;

  begin
    new_line;
    put_line("Reading a system with solutions ...");
    DoblDobl_System_and_Solutions_io.get(p,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    put_line("Initializing the systems container ...");
    DoblDobl_PolySys_Container.Initialize(p.all);
    put_line("Initializing the solutions container ...");
    DoblDobl_Solutions_Container.Initialize(sols);
    verbose := Prompt_for_Verbose;
    DoblDobl_ADE_Newton(verbose);
    newtsols := DoblDobl_Solutions_Container.Retrieve;
    put_line("The solutions after Newton's method :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
  end DoblDobl_Newton;

  procedure QuadDobl_Newton is

  -- DESCRIPTION :
  --   Reads a polynomial system and a corresponding list of solutions
  --   in quad double precision.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;

  begin
    new_line;
    put_line("Reading a system with solutions ...");
    QuadDobl_System_and_Solutions_io.get(p,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    put_line("Initializing the systems container ...");
    QuadDobl_PolySys_Container.Initialize(p.all);
    put_line("Initializing the solutions container ...");
    QuadDobl_Solutions_Container.Initialize(sols);
    verbose := Prompt_for_Verbose;
    QuadDobl_ADE_Newton(verbose);
    newtsols := QuadDobl_Solutions_Container.Retrieve;
    put_line("The solutions after Newton's method :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
  end QuadDobl_Newton;

  procedure Standard_Track_one_Path is

  -- DESCRIPTION :
  --   Reads a target and start system, with a solution in standard double
  --   precision and calls the path tracker.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    target,start : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;

  begin
    new_line;
    put_line("Reading a target system ...");
    get(target);
    new_line;
    put_line("Reading a start system with solutions ...");
    Standard_System_and_Solutions_io.get(start,sols);
    new_line;
    put("gamma = "); put(gamma); new_line;
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    put_line("Initializing the data for container copies ...");
    PHCpack_Operations.Store_Target_System(target.all);
    PHCpack_Operations.Store_Start_System(start.all);
    put_line("Initializing the solutions container ...");
    Standard_Solutions_Container.Initialize(sols);
    verbose := Prompt_for_Verbose;
    Standard_ADE_Track_One(verbose,gamma);
    newtsols := Standard_Solutions_Container.Retrieve;
    put_line("The solutions after path tracking :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
  end Standard_Track_one_Path;

  procedure DoblDobl_Track_one_Path is

  -- DESCRIPTION :
  --   Reads a target and start system, with a solution in double double
  --   precision and calls the path tracker.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    target,start : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;

  begin
    new_line;
    put_line("Reading a target system ...");
    get(target);
    new_line;
    put_line("Reading a start system with solutions ...");
    DoblDobl_System_and_Solutions_io.get(start,sols);
    new_line;
    put("gamma = "); put(gamma); new_line;
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    put_line("Initializing the data for container copies ...");
    PHCpack_Operations.Store_Target_System(target.all);
    PHCpack_Operations.Store_Start_System(start.all);
    put_line("Initializing the solutions container ...");
    DoblDobl_Solutions_Container.Initialize(sols);
    verbose := Prompt_for_Verbose;
    DoblDobl_ADE_Track_One(verbose,gamma);
    newtsols := DoblDobl_Solutions_Container.Retrieve;
    put_line("The solutions after path tracking :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
  end DoblDobl_Track_one_Path;

  procedure QuadDobl_Track_one_Path is

  -- DESCRIPTION :
  --   Reads a target and start system, with a solution in quad double
  --   precision and calls the path tracker.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    target,start : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;

  begin
    new_line;
    put_line("Reading a target system ...");
    get(target);
    new_line;
    put_line("Reading a start system with solutions ...");
    QuadDobl_System_and_Solutions_io.get(start,sols);
    new_line;
    put("gamma = "); put(gamma); new_line;
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    put_line("Initializing the data for container copies ...");
    PHCpack_Operations.Store_Target_System(target.all);
    PHCpack_Operations.Store_Start_System(start.all);
    put_line("Initializing the solutions container ...");
    QuadDobl_Solutions_Container.Initialize(sols);
    verbose := Prompt_for_Verbose;
    QuadDobl_ADE_Track_One(verbose,gamma);
    newtsols := QuadDobl_Solutions_Container.Retrieve;
    put_line("The solutions after path tracking :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
  end QuadDobl_Track_one_Path;

  procedure Standard_Track_many_Paths is

  -- DESCRIPTION :
  --   Reads a target and start system, with solutions in standard double
  --   precision, and calls the path tracker.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    target,start : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;

  begin
    new_line;
    put_line("Reading a target system ...");
    get(target);
    new_line;
    put_line("Reading a start system with solutions ...");
    Standard_System_and_Solutions_io.get(start,sols);
    new_line;
    put("gamma = "); put(gamma); new_line;
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    put_line("Initializing the data for container copies ...");
    PHCpack_Operations.Store_Target_System(target.all);
    PHCpack_Operations.Store_Start_System(start.all);
    put_line("Initializing the solutions container ...");
    Standard_Solutions_Container.Initialize(sols);
    verbose := Prompt_for_Verbose;
    tstart(timer);
    Standard_ADE_Track_Many(verbose,gamma);
    tstop(timer);
    newtsols := Standard_Solutions_Container.Retrieve;
    put_line("The solutions after path tracking :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
    new_line;
    print_times
      (standard_output,timer,"tracking many paths in double precision");
  end Standard_Track_many_Paths;

  procedure DoblDobl_Track_many_Paths is

  -- DESCRIPTION :
  --   Reads a target and start system, with solutions in double double
  --   precision, and calls the path tracker.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    target,start : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;

  begin
    new_line;
    put_line("Reading a target system ...");
    get(target);
    new_line;
    put_line("Reading a start system with solutions ...");
    DoblDobl_System_and_Solutions_io.get(start,sols);
    new_line;
    put("gamma = "); put(gamma); new_line;
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    put_line("Initializing the data for container copies ...");
    PHCpack_Operations.Store_Target_System(target.all);
    PHCpack_Operations.Store_Start_System(start.all);
    put_line("Initializing the solutions container ...");
    DoblDobl_Solutions_Container.Initialize(sols);
    verbose := Prompt_for_Verbose;
    tstart(timer);
    DoblDobl_ADE_Track_Many(verbose,gamma);
    tstop(timer);
    newtsols := DoblDobl_Solutions_Container.Retrieve;
    put_line("The solutions after path tracking :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
    new_line;
    print_times
      (standard_output,timer,"tracking many paths in double double precision");
  end DoblDobl_Track_many_Paths;

  procedure QuadDobl_Track_many_Paths is

  -- DESCRIPTION :
  --   Reads a target and start system, with solutions in double double
  --   precision, and calls the path tracker.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    target,start : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    verbose : integer32;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;

  begin
    new_line;
    put_line("Reading a target system ...");
    get(target);
    new_line;
    put_line("Reading a start system with solutions ...");
    QuadDobl_System_and_Solutions_io.get(start,sols);
    new_line;
    put("gamma = "); put(gamma); new_line;
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    put_line("Initializing the data for container copies ...");
    PHCpack_Operations.Store_Target_System(target.all);
    PHCpack_Operations.Store_Start_System(start.all);
    put_line("Initializing the solutions container ...");
    QuadDobl_Solutions_Container.Initialize(sols);
    verbose := Prompt_for_Verbose;
    tstart(timer);
    QuadDobl_ADE_Track_Many(verbose,gamma);
    tstop(timer);
    newtsols := QuadDobl_Solutions_Container.Retrieve;
    put_line("The solutions after path tracking :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
    new_line;
    print_times
      (standard_output,timer,"tracking many paths in quad double precision");
  end QuadDobl_Track_many_Paths;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precsion and then calls the test.

    tst,prc : character;

  begin
    new_line;
    put_line("MENU to test the algorithmic differentiation :");
    put_line("  0. running Newton's method;");
    put_line("  1. tracking one path;");
    put_line("  2. tracking many paths;");
    put("Type 0, 1, or 2 to select a test : ");
    Ask_Alternative(tst,"012");
    new_line;
    put_line("MENU to test the ADE code to track many paths :");
    put_line("  0. in standard double precision;");
    put_line("  1. in double double precision;");
    put_line("  2. in quad double precision.");
    put("Type 0, 1, or 2 to make your choice : ");
    Ask_Alternative(prc,"012");
    case tst is
      when '0' =>
        case prc is
          when '0' => Standard_Newton;
          when '1' => DoblDobl_Newton;
          when '2' => QuadDobl_Newton;
          when others => null;
        end case;
      when '1' =>
        case prc is
          when '0' => Standard_Track_one_Path;
          when '1' => DoblDobl_Track_one_Path;
          when '2' => QuadDobl_Track_one_Path;
          when others => null;
        end case;
      when '2' =>
        case prc is
          when '0' => Standard_Track_many_Paths;
          when '1' => DoblDobl_Track_many_Paths;
          when '2' => QuadDobl_Track_many_Paths;
          when others => null;
        end case;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_ademanypaths;
