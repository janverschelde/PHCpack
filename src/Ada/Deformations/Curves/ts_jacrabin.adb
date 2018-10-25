with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Standard_Deflation_Trees_io;
with DoblDobl_Deflation_Trees_io;
with QuadDobl_Deflation_Trees_io;
with Jacobian_Rabinowitsch_Trick;       use Jacobian_Rabinowitsch_Trick;

procedure ts_jacrabin is

-- DESSCRIPTION :
--   Interactive test on the development of the formulation of the
--   polynomial system along the lines of the trick of Rabinowitsch
--   to move singular solutions to infinity.

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then makes
  --   the system augmented with the Jacobian matrix,
  --   in standard double precision.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    nvar : natural32;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    new_line;
    put_line("Your system :"); put(lp.all);
    new_line;
    put_line("The Rabinowitsch trick applied to the Jacobian :");
    nvar := Standard_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    Standard_Deflation_Trees_io.Add_Multiplier_Symbols(1,nvar);
    Symbol_Table.Enlarge(1);
    declare
      sb : Symbol_Table.Symbol;
      jacrbp : constant Standard_Complex_Poly_Systems.Poly_Sys
             := Jacobian_Rabinowitsch(lp.all);
    begin
      sb := (sb'range => ' ');
      sb(1) := 'y';
      sb(2) := 'r';
      sb(3) := 'b';
      Symbol_Table.add(sb);
      put(jacrbp);
    end;
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then makes
  --   the system augmented with the Jacobian matrix,
  --   in double double precision.

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    nvar : natural32;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    new_line;
    put_line("Your system :"); put(lp.all);
    new_line;
    put_line("The Rabinowitsch trick applied to the Jacobian :");
    nvar := DoblDobl_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    DoblDobl_Deflation_Trees_io.Add_Multiplier_Symbols(1,nvar);
    Symbol_Table.Enlarge(1);
    declare
      sb : Symbol_Table.Symbol;
      jacrbp : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
             := Jacobian_Rabinowitsch(lp.all);
    begin
      sb := (sb'range => ' ');
      sb(1) := 'y';
      sb(2) := 'r';
      sb(3) := 'b';
      Symbol_Table.add(sb);
      put(jacrbp);
    end;
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then makes
  --   the system augmented with the Jacobian matrix,
  --   in quad double precision.

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    nvar : natural32;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    new_line;
    put_line("Your system :"); put(lp.all);
    new_line;
    put_line("The Rabinowitsch trick applied to the Jacobian :");
    nvar := QuadDobl_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    QuadDobl_Deflation_Trees_io.Add_Multiplier_Symbols(1,nvar);
    Symbol_Table.Enlarge(1);
    declare
      sb : Symbol_Table.Symbol;
      jacrbp : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
             := Jacobian_Rabinowitsch(lp.all);
    begin
      sb := (sb'range => ' ');
      sb(1) := 'y';
      sb(2) := 'r';
      sb(3) := 'b';
      Symbol_Table.add(sb);
      put(jacrbp);
    end;
  end QuadDobl_Main;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU for the precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_jacrabin;
