with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Drivers_for_m_homogenization;       use Drivers_for_m_homogenization;

procedure ts_drivmhom is

-- DESCRIPTION :
--   Reads a polynomial system and calls the driver.

  file : file_type;
  lp : Link_to_Poly_Sys;
  b : natural64 := 0;

begin
  get(lp);
  declare
    q : Poly_Sys(lp'range);
    qsols : Solution_List;
  begin
    put_line("Reading the output file.");
    Read_Name_and_Create_File(file);
    Driver_for_m_homogenization(file,lp.all,b,q,qsols);
  end;
end ts_drivmhom;
