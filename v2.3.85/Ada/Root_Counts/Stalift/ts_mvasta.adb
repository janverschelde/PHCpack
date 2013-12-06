with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Drivers_for_MixedVol_Algorithm;     use Drivers_for_MixedVol_Algorithm;

procedure ts_mvasta is

-- DESCRIPTION :
--   Create a random coefficient start system, using the MixedVol algorithm.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system, computes its supports,
  --   and then calls the routine to compute the mixed volume.

    lp : Link_to_Poly_Sys;
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
  end Main;

begin
  Main;
end ts_mvasta;
