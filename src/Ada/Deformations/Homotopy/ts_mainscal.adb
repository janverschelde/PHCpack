with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Multprec_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Scaling_Methods;

procedure ts_mainscal is

-- DESCRIPTION :
--   Calls the main procedures to scale a polynomial system.

  procedure Standard_Scale is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and calls
  --   the driver to scale in standard double precision.

    use Standard_Complex_Vectors;
    use Standard_Complex_Poly_Systems;

    file : file_type;
    lp : Link_to_Poly_Sys;
    scvc : Link_to_Vector;
    bas : natural32 := 2;

  begin
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Scaling_Methods.Main(file,lp.all,bas,scvc);
  end Standard_Scale;

  procedure DoblDobl_Scale is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and calls
  --   the driver to scale in double double precision.

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Poly_Systems;

    file : file_type;
    lp : Link_to_Poly_Sys;
    scvc : Link_to_Vector;
    bas : natural32 := 2;

  begin
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Scaling_Methods.Main(file,lp.all,bas,scvc);
  end DoblDobl_Scale;

  procedure QuadDobl_Scale is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and calls
  --   the driver to scale in quad double precision.

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Poly_Systems;

    file : file_type;
    lp : Link_to_Poly_Sys;
    scvc : Link_to_Vector;
    bas : natural32 := 2;

  begin
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Scaling_Methods.Main(file,lp.all,bas,scvc);
  end QuadDobl_Scale;

  procedure Multprec_Scale is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and calls
  --   the driver to scale in arbitrary multiprecision.

    use Multprec_Complex_Vectors;
    use Multprec_Complex_Poly_Systems;

    file : file_type;
    lp : Link_to_Poly_Sys;
    scvc : Link_to_Vector;
    bas : natural32 := 2;

  begin
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Scaling_Methods.Main(file,lp.all,bas,scvc);
  end Multprec_Scale;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user to select the level of precision
  --   and then calls the corresponding driver procedure.

    ans : character;

  begin
    new_line;
    put_line("MENU to select the precision : ");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision;");
    put_line("  3. arbitrary multiprecision.");
    put("Type 0, 1, 2, or 3 to select the precision : ");
    Ask_Alternative(ans,"0123");
    case ans is
      when '0' => Standard_Scale;
      when '1' => DoblDobl_Scale;
      when '2' => QuadDobl_Scale;
      when '3' => Multprec_Scale;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mainscal;
