with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_Strings;      use Standard_Complex_Poly_Strings;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;    use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Strings;      use DoblDobl_Complex_Poly_Strings;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;    use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Strings;      use QuadDobl_Complex_Poly_Strings;

procedure ts_rwspol is

-- DESCRIPTION :
--   Test on reading a polynomial from file
--   and writing it to string.

  function Standard_Read return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Prompts for a polynomial system that consists
  --   of one single polynomial with coefficients in standard double precision.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    get(lp);
    return lp(lp'first);
  end Standard_Read;

  function DoblDobl_Read return DoblDobl_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Prompts for a polynomial system that consists
  --   of one single polynomial with coefficients in double double precision.

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    get(lp);
    return lp(lp'first);
  end DoblDobl_Read;

  function QuadDobl_Read return QuadDobl_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Prompts for a polynomial system that consists
  --   of one single polynomial with coefficients in quad double precision.

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    get(lp);
    return lp(lp'first);
  end QuadDobl_Read;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Reads a polynomial and writes it to a string.

    p : Standard_Complex_Polynomials.Poly := Standard_Read;

  begin
    put_line("The polynomial :"); put(p); new_line;
    declare
      s : constant string := Write(p);
    begin
      put_line("The string : "); put_line(s);
      put("string size : "); put(natural32(s'last)); new_line;
      put("size limit : "); put(Size_Limit(p)); new_line;
    end;
    Standard_Complex_Polynomials.Clear(p);
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Reads a polynomial and writes it to a string.

    p : DoblDobl_Complex_Polynomials.Poly := DoblDobl_Read;

  begin
    put_line("The polynomial :"); put(p); new_line;
    declare
      s : constant string := Write(p);
    begin
      put_line("The string : "); put_line(s);
      put("string size : "); put(natural32(s'last)); new_line;
      put("size limit : "); put(Size_Limit(p)); new_line;
    end;
    DoblDobl_Complex_Polynomials.Clear(p);
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Reads a polynomial and writes it to a string.

    p : QuadDobl_Complex_Polynomials.Poly := QuadDobl_Read;

  begin
    put_line("The polynomial :"); put(p); new_line;
    declare
      s : constant string := Write(p);
    begin
      put_line("The string : "); put_line(s);
      put("string size : "); put(natural32(s'last)); new_line;
      put("size limit : "); put(Size_Limit(p)); new_line;
    end;
    QuadDobl_Complex_Polynomials.Clear(p);
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision
  --   and then calls the proper test.

    ans : character;

  begin
    new_line;
    put_line("MENU to select the working precision : ");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision;");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_rwspol;
