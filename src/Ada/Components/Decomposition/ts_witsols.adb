with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Embeddings_and_Cascades;            use Embeddings_and_Cascades;
with Standard_Witness_Solutions;
with DoblDobl_Witness_Solutions;
with QuadDobl_Witness_Solutions;

procedure ts_witsols is

-- DESCRIPTION :
--   Test on the package to store witness solutions.

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   and the prompts for the top dimension.

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    nq,nv,topdim,lowdim : natural32;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    nq := natural32(lp'last);
    nv := Number_of_Unknowns(lp(lp'first));
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    Standard_Witness_Solutions.Initialize(topdim);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   and the prompts for the top dimension.

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    nq,nv,topdim,lowdim : natural32;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    nq := natural32(lp'last);
    nv := Number_of_Unknowns(lp(lp'first));
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    DoblDobl_Witness_Solutions.Initialize(topdim);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   and the prompts for the top dimension.

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    nq,nv,topdim,lowdim : natural32;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    nq := natural32(lp'last);
    nv := Number_of_Unknowns(lp(lp'first));
    Prompt_for_Top_Dimension(nq,nv,topdim,lowdim);
    QuadDobl_Witness_Solutions.Initialize(topdim);
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision
  --   and then calls the corresponding main test.

    ans : character;

  begin
    new_line;
    put_line("MENU for the precision :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, 2 to select the precision : ");
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
end ts_witsols;
