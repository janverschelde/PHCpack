with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;
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
    sols : Standard_Complex_Solutions.Solution_List;
    nvar : natural32;
    ans : character;

    use Standard_Complex_Solutions;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    Standard_System_and_Solutions_io.get(lp,sols);
    new_line;
    put_line("Your system :"); put(lp.all);
    new_line;
    put_line("The Rabinowitsch trick applied to the Jacobian :");
    nvar := Standard_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    Add_Trick_Symbols(nvar);
    declare
      jacrbp : constant Standard_Complex_Poly_Systems.Poly_Sys
             := Jacobian_Rabinowitsch(lp.all);
      jrbsols : constant Solution_List := Jacobian_Rabinowitsch(sols);
    begin
      put(jacrbp);
      new_line;
      put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("The extended solutions : ");
        put(standard_output,Length_Of(jrbsols),
            natural32(Head_Of(jrbsols).n),jrbsols);
      end if;
    end;
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then makes
  --   the system augmented with the Jacobian matrix,
  --   in double double precision.

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    nvar : natural32;
    ans : character;

    use DoblDobl_Complex_Solutions;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    DoblDobl_System_and_Solutions_io.get(lp,sols);
    new_line;
    put_line("Your system :"); put(lp.all);
    new_line;
    put_line("The Rabinowitsch trick applied to the Jacobian :");
    nvar := DoblDobl_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    Add_Trick_Symbols(nvar);
    declare
      jacrbp : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
             := Jacobian_Rabinowitsch(lp.all);
      jrbsols : constant Solution_List := Jacobian_Rabinowitsch(sols);
    begin
      put(jacrbp);
      new_line;
      put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("The extended solutions : ");
        put(standard_output,Length_Of(jrbsols),
            natural32(Head_Of(jrbsols).n),jrbsols);
      end if;
    end;
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then makes
  --   the system augmented with the Jacobian matrix,
  --   in quad double precision.

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    nvar : natural32;
    ans : character;

    use QuadDobl_Complex_Solutions;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    QuadDobl_System_and_Solutions_io.get(lp,sols);
    new_line;
    put_line("Your system :"); put(lp.all);
    new_line;
    put_line("The Rabinowitsch trick applied to the Jacobian :");
    nvar := QuadDobl_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    Add_Trick_Symbols(nvar);
    declare
      jacrbp : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
             := Jacobian_Rabinowitsch(lp.all);
      jrbsols : constant Solution_List := Jacobian_Rabinowitsch(sols);
    begin
      put(jacrbp);
      new_line;
      put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("The extended solutions : ");
        put(standard_output,Length_Of(jrbsols),
            natural32(Head_Of(jrbsols).n),jrbsols);
      end if;
    end;
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the precision and launches the tests.

    ans : constant character := Prompt_for_Precision;

  begin
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
