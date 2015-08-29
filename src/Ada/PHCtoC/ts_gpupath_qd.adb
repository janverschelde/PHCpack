with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with QuadDobl_PolySys_Container;
with QuadDobl_Solutions_Container;
with PHCpack_Operations;
with QuadDobl_Accelerated_Trackers;      use QuadDobl_Accelerated_Trackers;

procedure ts_gpupath_qd is

-- DESCRIPTION :
--   Procedure to test the development of the acceleration of Newton's method
--   and the acceleration of path trackers in quad double precision.

  function Prompt_for_Mode return integer32 is

  -- DESCRIPTION :
  --   Shows the user the menu and prompts for a mode of execution,
  --   returned as 0, 1, or 2.

    ans : character;

  begin
    new_line;
    put_line("MENU for the mode of execution : ");
    put_line("  0. Both the CPU and GPU will execute.");
    put_line("  1. Only the CPU will execute.");
    put_line("  2. Only the GPU will execute.");
    put("Type 0, 1, or 2 : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => return 0;
      when '1' => return 1;
      when '2' => return 2;
      when others => return -1;
    end case;
  end Prompt_for_Mode;

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

  procedure QuadDobl_Newton is

  -- DESCRIPTION :
  --   Reads a polynomial system and a corresponding list of solutions
  --   in standard double precision.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    execmode,verbose : integer32;

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
    execmode := Prompt_for_Mode;
    verbose := Prompt_for_Verbose;
    QuadDobl_GPU_Newton(execmode,verbose);
    newtsols := QuadDobl_Solutions_Container.Retrieve;
    put_line("The solutions after Newton's method :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
  end QuadDobl_Newton;

  procedure QuadDobl_Track_one_Path is

  -- DESCRIPTION :
  --   Reads a target and start system, with a solution in quad double
  --   precision and calls the accelerated path tracker.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    target,start : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    execmode,verbose : integer32;
    gamma : constant Complex_Number
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
    execmode := Prompt_for_Mode;
    verbose := Prompt_for_Verbose;
    QuadDobl_GPU_Track_One(execmode,verbose,gamma);
    newtsols := QuadDobl_Solutions_Container.Retrieve;
    put_line("The solutions after path tracking :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
  end QuadDobl_Track_one_Path;

  procedure QuadDobl_Track_many_Paths is

  -- DESCRIPTION :
  --   Reads a target and start system, with solutions in quad double
  --   precision, and calls the accelerated path tracker.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    target,start : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    execmode,verbose : integer32;
    gamma : constant Complex_Number
          := Standard_Random_Numbers.Random1;

  begin
    new_line;
    put_line("Reading a target system ...");
    get(target);
    new_line;
    put_line("Reading a start system with solutions ...");
    QuadDobl_System_and_Solutions_io.get(start,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    put_line("Initializing the data for container copies ...");
    PHCpack_Operations.Store_Target_System(target.all);
    PHCpack_Operations.Store_Start_System(start.all);
    put_line("Initializing the solutions container ...");
    QuadDobl_Solutions_Container.Initialize(sols);
    execmode := Prompt_for_Mode;
    verbose := Prompt_for_Verbose;
    QuadDobl_GPU_Track_Many(execmode,verbose,gamma);
    newtsols := QuadDobl_Solutions_Container.Retrieve;
    put_line("The solutions after path tracking :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
  end QuadDobl_Track_many_Paths;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test GPU acceleration in quad double precision :");
    put_line("  0. apply Newton's method on one single solution;");
    put_line("  1. track one single solution path from start to target;");
    put_line("  2. track many solution paths from start to target.");
    put("Type 0, 1, or 2 to make your choice : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => QuadDobl_Newton;
      when '1' => QuadDobl_Track_one_Path;
      when '2' => QuadDobl_Track_Many_Paths;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_gpupath_qd;
