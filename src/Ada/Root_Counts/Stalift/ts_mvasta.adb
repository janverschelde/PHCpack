with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Drivers_for_MixedVol_Algorithm;     use Drivers_for_MixedVol_Algorithm;

procedure ts_mvasta is

-- DESCRIPTION :
--   Create a random coefficient start system, using the MixedVol algorithm.

  procedure Standard_Polyhedral_Homotopies is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system, computes its supports,
  --   and then calls the routine to compute the mixed volume.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    cellfile,startfile : file_type;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    new_line;
    put_line("Reading the name of the file to write the cells ...");
    Read_Name_and_Create_File(cellfile);
    new_line;
    put_line("Reading the name of the file to write the start system ...");
    Read_Name_and_Create_File(startfile);
    new_line;
    Polyhedral_Homotopies(cellfile,startfile,lp.all);
  end Standard_Polyhedral_Homotopies;

  procedure DoblDobl_Polyhedral_Homotopies is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system, computes its supports,
  --   and then calls the routine to compute the mixed volume.

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    cellfile,startfile : file_type;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    new_line;
    put_line("Reading the name of the file to write the cells ...");
    Read_Name_and_Create_File(cellfile);
    new_line;
    put_line("Reading the name of the file to write the start system ...");
    Read_Name_and_Create_File(startfile);
    new_line;
    Polyhedral_Homotopies(cellfile,startfile,lp.all);
  end DoblDobl_Polyhedral_Homotopies;

  procedure QuadDobl_Polyhedral_Homotopies is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system, computes its supports,
  --   and then calls the routine to compute the mixed volume.

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    cellfile,startfile : file_type;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    new_line;
    put_line("Reading the name of the file to write the cells ...");
    Read_Name_and_Create_File(cellfile);
    new_line;
    put_line("Reading the name of the file to write the start system ...");
    Read_Name_and_Create_File(startfile);
    new_line;
    Polyhedral_Homotopies(cellfile,startfile,lp.all);
  end QuadDobl_Polyhedral_Homotopies;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to call the polyhedral homotopies after MixedVol :");
    put_line("  1. use standard double precision for polyhedral homotopies;");
    put_line("  2. use double double precision for polyhedral homotopies;");
    put_line("  3. use quad double precision for polyhedral homotopies;");
    put("Type 1, 2, or 3 to select the precision : ");
    Ask_Alternative(ans,"123");
    case ans is
      when '1' => Standard_Polyhedral_Homotopies;
      when '2' => DoblDobl_Polyhedral_Homotopies;
      when '3' => QuadDobl_Polyhedral_Homotopies;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mvasta;
