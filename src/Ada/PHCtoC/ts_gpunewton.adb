with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with Standard_PolySys_Container;
with Standard_Solutions_Container;

procedure ts_gpunewton is

-- DESCRIPTION :
--   Procedure to test the development of the acceleration of Newton's method
--   as called from an Ada procedure.  The Ada procedure remains in control.
--   Data is passed to the C++ code via the systems and solutions containers.

  procedure Standard_GPU_Newton is

  -- DESCRIPTION :
  --   Calls the accelerated Newton's method in standard double precision.

    return_of_call : integer32;

    function newton return integer32;
    pragma import(C, newton, "gpunewton_d");

  begin
    return_of_call := newton;
  end Standard_GPU_Newton;

  procedure Standard_Newton is

  -- DESCRIPTION :
  --   Reads a polynomial system and a corresponding list of solutions
  --   in standard double precision.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;

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
    Standard_GPU_Newton;
    newtsols := Standard_Solutions_Container.Retrieve;
    put_line("The solutions after Newton's method :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
  end Standard_Newton;

  procedure Main is
  begin
    Standard_Newton;
  end Main;

begin
  Main;
end ts_gpunewton;
