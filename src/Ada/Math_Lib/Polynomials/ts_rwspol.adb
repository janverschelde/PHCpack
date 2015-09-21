with text_io;                            use text_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_Strings;      use Standard_Complex_Poly_Strings;

procedure ts_rwspol is

-- DESCRIPTION :
--   Test on reading a polynomial from file
--   and writing it to string.

  function read return Poly is

  -- DESCRIPTION :
  --   Prompts for a polynomial system that consists
  --   of one single polynomial.

    lp : Link_to_Poly_Sys;

  begin
    get(lp);
    return lp(lp'first);
  end read;

  procedure Main is

  -- DESCRIPTION :
  --   Reads a polynomial and writes it to a string.

    p : Poly := read;

  begin
    put_line("The polynomial :"); put(p); new_line;
    declare
      s : constant string := Write(p);
    begin
      put_line("The string : "); put_line(s);
    end;
  end Main;

begin
  Main;
end ts_rwspol;
