with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Symbol_Table;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Poly_Systems;
with Standard_Poly_Laur_Convertors;      use Standard_Poly_Laur_Convertors;
with DoblDobl_Poly_Laur_Convertors;      use DoblDobl_Poly_Laur_Convertors;
with QuadDobl_Poly_Laur_Convertors;      use QuadDobl_Poly_Laur_Convertors;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laurentials;
with Random_Polynomial_Systems;          use Random_Polynomial_Systems;
with Regular_Newton_Puiseux;             use Regular_Newton_Puiseux;

procedure ts_puiseux is

-- DESCRIPTION :
--   Development of the Newton-Puiseux algorithm,
--   for regular solution curves, defined by complete intersections,
--   in Noether position, with sufficiently general coefficients.

  procedure Standard_Random_Test is

  -- DESCRIPTION :
  --   Prompts the user for the parameters to generate a system
  --   with random coefficients in standard double precision.
  --   The series are computed for the generated system.

    n,d,m : natural32 := 0;
    c : constant natural32 := 0;
    e : integer32 := 0;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

    use Standard_Complex_Laurentials;

  begin
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give the maximal degree : "); get(d);
    put("Give number of monomials (0 for dense): "); get(m);
    e := integer32(n) - 1;
    Standard_Generate_and_Show(n,d,m,c,e,lp);
    declare
      q : constant Standard_Complex_Laur_Systems.Laur_Sys(lp'range)
        := Polynomial_to_Laurent_System(lp.all);
      nv : constant integer32 := integer32(Number_of_Unknowns(q(q'first)));
    begin
      Standard_Test(q,q'last,nv);
    end;
  end Standard_Random_Test;

  procedure DoblDobl_Random_Test is

  -- DESCRIPTION :
  --   Prompts the user for the parameters to generate a system
  --   with random coefficients in double double precision.
  --   The series are computed for the generated system.

    n,d,m : natural32 := 0;
    c : constant natural32 := 0;
    e : integer32 := 0;
    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

    use DoblDobl_Complex_Laurentials;

  begin
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give the maximal degree : "); get(d);
    put("Give number of monomials (0 for dense): "); get(m);
    e := integer32(n) - 1;
    DoblDobl_Generate_and_Show(n,d,m,c,e,lp);
    declare
      q : constant DoblDobl_Complex_Laur_Systems.Laur_Sys(lp'range)
        := Polynomial_to_Laurent_System(lp.all);
      nv : constant integer32 := integer32(Number_of_Unknowns(q(q'first)));
    begin
      DoblDobl_Test(q,q'last,nv);
    end;
  end DoblDobl_Random_Test;

  procedure QuadDobl_Random_Test is

  -- DESCRIPTION :
  --   Prompts the user for the parameters to generate a system
  --   with random coefficients in quad double precision.
  --   The series are computed for the generated system.

    n,d,m : natural32 := 0;
    c : constant natural32 := 0;
    e : integer32 := 0;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

    use QuadDobl_Complex_Laurentials;

  begin
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give the maximal degree : "); get(d);
    put("Give number of monomials (0 for dense): "); get(m);
    e := integer32(n) - 1;
    QuadDobl_Generate_and_Show(n,d,m,c,e,lp);
    declare
      q : constant QuadDobl_Complex_Laur_Systems.Laur_Sys(lp'range)
        := Polynomial_to_Laurent_System(lp.all);
      nv : constant integer32 := integer32(Number_of_Unknowns(q(q'first)));
    begin
      QuadDobl_Test(q,q'last,nv);
    end;
  end QuadDobl_Random_Test;

  function Prompt_for_Precision return character is

  -- DESCRIPTION :
  --   Displays the menu for the available precision and
  --   returns '0', '1', or '2' for double, double double,
  --   or quad double precision.

    res : character;

  begin
    new_line;
    put_line("MENU to set the precision level :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(res,"012");
    return res;
  end Prompt_for_Precision;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision
  --   and then launches the proper main driver.

    prc : constant character := Prompt_for_Precision;
    ans : character;

  begin
    new_line;
    put("Generate a random polynomial system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      case prc is
        when '0' => Standard_Random_Test;
        when '1' => DoblDobl_Random_Test;
        when '2' => QuadDobl_Random_Test;
        when others => null;
      end case;
    else
      case prc is
        when '0' => Standard_Main;
        when '1' => DoblDobl_Main;
        when '2' => QuadDobl_Main;
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_puiseux;
