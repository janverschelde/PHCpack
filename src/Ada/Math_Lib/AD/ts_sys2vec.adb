with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Standard_Polynomial_Vectors;
with Standard_Polynomial_Vectors_io;    use Standard_Polynomial_Vectors_io;
with DoblDobl_Polynomial_Vectors;
with DoblDobl_Polynomial_Vectors_io;    use DoblDobl_Polynomial_Vectors_io;
with QuadDobl_Polynomial_Vectors;
with QuadDobl_Polynomial_Vectors_io;    use QuadDobl_Polynomial_Vectors_io;
with System_Vector_Convertors;          use System_Vector_Convertors;

procedure ts_sys2vec is

-- DESCRIPTION :
--   Development of the conversion from a polynomial system into
--   the polynomial vectors for efficient evaluation and differentiation.

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then converts it
  --   into polynomial vectors, in standard double precision.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    pv : Standard_Polynomial_Vectors.Link_to_Polynomial_Vector;
 
  begin
    put_line("Reading a polynomial system ..."); get(lp);
    pv := Convert(lp);
    put_line("The converted system :"); put(pv);
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then converts it
  --   into polynomial vectors, in double double precision.

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    pv : DoblDobl_Polynomial_Vectors.Link_to_Polynomial_Vector;
 
  begin
    put_line("Reading a polynomial system ..."); get(lp);
    pv := Convert(lp);
    put_line("The converted system :"); put(pv);
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then converts it
  --   into polynomial vectors, in quad double precision.

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    pv : QuadDobl_Polynomial_Vectors.Link_to_Polynomial_Vector;
 
  begin
    put_line("Reading a polynomial system ..."); get(lp);
    pv := Convert(lp);
    put_line("The converted system :"); put(pv);
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a precision and then launches
  --   the corresponding test.

    ans : character;

  begin
    new_line;
    put_line("MENU to select the precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
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
end ts_sys2vec;
