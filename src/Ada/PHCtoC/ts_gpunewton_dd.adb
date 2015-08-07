with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with DoblDobl_PolySys_Container;
with DoblDobl_Solutions_Container;

procedure ts_gpunewton_dd is

-- DESCRIPTION :
--   Procedure to test the development of the acceleration of Newton's method
--   as called from an Ada procedure.  The Ada procedure remains in control.
--   Data is passed to the C++ code via the systems and solutions containers.

  procedure DoblDobl_GPU_Newton is

  -- DESCRIPTION :
  --   Calls the accelerated Newton's method in double double precision.

    return_of_call : integer32;

    function newton return integer32;
    pragma import(C, newton, "gpunewton_dd");

  begin
    return_of_call := newton;
  end DoblDobl_GPU_Newton;

  procedure DoblDobl_Newton is

  -- DESCRIPTION :
  --   Reads a polynomial system and a corresponding list of solutions
  --   in double double precision.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;

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
    DoblDobl_GPU_Newton;
    newtsols := DoblDobl_Solutions_Container.retrieve;
    put_line("The updated solutions :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
  end DoblDobl_Newton;

  procedure Main is
  begin
    DoblDobl_Newton;
  end Main;

begin
  Main;
end ts_gpunewton_dd;
