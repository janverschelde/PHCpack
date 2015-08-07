with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with QuadDobl_PolySys_Container;
with QuadDobl_Solutions_Container;

procedure ts_gpunewton_qd is

-- DESCRIPTION :
--   Procedure to test the development of the acceleration of Newton's method
--   as called from an Ada procedure.  The Ada procedure remains in control.
--   Data is passed to the C++ code via the systems and solutions containers.

  procedure QuadDobl_GPU_Newton is

  -- DESCRIPTION :
  --   Calls the accelerated Newton's method in quad double precision.

    return_of_call : integer32;

    function newton return integer32;
    pragma import(C, newton, "gpunewton_qd");

  begin
    return_of_call := newton;
  end QuadDobl_GPU_Newton;

  procedure QuadDobl_Newton is

  -- DESCRIPTION :
  --   Reads a polynomial system and a corresponding list of solutions
  --   in quad double precision.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols,newtsols : Solution_List;

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
    QuadDobl_GPU_Newton;
    newtsols := QuadDobl_Solutions_Container.Retrieve;
    put_line("The updated solutions :");
    put(standard_output,Length_Of(newtsols),
        natural32(Head_Of(newtsols).n),newtsols);
  end QuadDobl_Newton;

  procedure Main is
  begin
    QuadDobl_Newton;
  end Main;

begin
  Main;
end ts_gpunewton_qd;
