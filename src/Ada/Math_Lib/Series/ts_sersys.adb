with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;
with Series_and_Polynomials;
with Series_and_Polynomials_io;
with Standard_Series_Poly_Systems;
with Standard_Series_Poly_SysFun;
with Standard_Series_Jaco_Matrices;

procedure ts_sersys is

-- DESCRIPTION :
--   Tests the methods on systems of series polynomials.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a system of series polynomials.

    ls : Standard_Series_Poly_Systems.Link_to_Poly_Sys;
    ix : integer32 := 0;

  begin
    new_line;
    put("Give the index of the series variable : "); get(ix);
    new_line;
    Series_and_Polynomials.get(ls,ix);
    new_line;
    put_line("The polynomial system : ");
    Series_and_Polynomials.put(ls.all,ix);
  end Main;

begin
  Main;
end ts_sersys;
