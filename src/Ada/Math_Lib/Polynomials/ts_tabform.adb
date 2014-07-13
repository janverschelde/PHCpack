with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Tableau_Formats;
with DoblDobl_Tableau_Formats;
with QuadDobl_Tableau_Formats;

procedure ts_tabform is

-- DESCRIPTION :
--   This simple procedure allows the conversion between polynomials in
--   symbolic (PHCpack) format and a more basic tableau format.

  ans : character;

begin
  new_line;
  put_line("MENU to convert polynomials to tableau format : ");
  put_line("  1. with coefficients in standard double precision");
  put_line("  2. with coefficients in double double precision");
  put_line("  3. with coefficients in quad double precision");
  put("Type 1, 2, or 3 to select the precision : ");
  Ask_Alternative(ans,"123");
  case ans is
    when '1' => Standard_Tableau_Formats.Main_Interactive_Driver;
    when '2' => DoblDobl_Tableau_Formats.Main_Interactive_Driver;
    when '3' => QuadDobl_Tableau_Formats.Main_Interactive_Driver;
    when others => null;
  end case;
end ts_tabform;
