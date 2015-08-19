with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with PHCpack_Operations;
with DoblDobl_Solutions_Container;

procedure ts_gpumanypaths_dd is

-- DESCRIPTION :
--   Procedure to test the development of the acceleration of the tracking
--   of many solution paths, in double double precision.
--   The Ada procedure remains in control.
--   Data is passed to the C++ code via the systems and solutions containers.

  procedure DoblDobl_GPU_Track ( execmode,verbose : integer32 ) is

  -- DESCRIPTION :
  --   Calls the accelerated path tracker double double precision.
  --   The first input parameter is the mode of execution: 0, 1, or 2.
  --   If verbose > 0, then additional output is written to screen.

    return_of_call : integer32;

    function track ( m,v : integer32 ) return integer32;
    pragma import(C, track, "gpumanypaths_dd");

  begin
    return_of_call := track(execmode,verbose);
  end DoblDobl_GPU_Track;

  function Prompt_for_Mode return integer32 is

  -- DESCRIPTION :
  --   Shows the user the menu and prompts for a mode of execution,
  --   returned as 0, 1, or 2.

    ans : character;

  begin
    new_line;
    put_line("MENU for the mode of execution : ");
    put_line("  0. Both the CPU and GPU will run the path tracker.");
    put_line("  1. Only the CPU will execute the path tracker.");
    put_line("  2. Only the GPU will execute the path tracker.");
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

  procedure DoblDobl_Track is

  -- DESCRIPTION :
  --   Reads a target and start system, with solutions in double double
  --   precision.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    target,start : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;
    execmode,verbose : integer32;

  begin
    new_line;
    put_line("Reading a target system ...");
    get(target);
    new_line;
    put_line("Reading a start system with solutions ...");
    DoblDobl_System_and_Solutions_io.get(start,sols);
    new_line;
    put("Read "); put(Length_Of(sols),1); put_line(" solutions.");
    put_line("Initializing the data for container copies ...");
    PHCpack_Operations.Store_Target_System(target.all);
    PHCpack_Operations.Store_Start_System(start.all);
    put_line("Initializing the solutions container ...");
    DoblDobl_Solutions_Container.Initialize(sols);
    execmode := Prompt_for_Mode;
    verbose := Prompt_for_Verbose;
    DoblDobl_GPU_Track(execmode,verbose);
    newtsols := DoblDobl_Solutions_Container.Retrieve;
    put_line("The solutions after path tracking :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
  end DoblDobl_Track;

  procedure Main is
  begin
    DoblDobl_Track;
  end Main;

begin
  Main;
end ts_gpumanypaths_dd;
