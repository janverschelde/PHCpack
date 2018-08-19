with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Embeddings_and_Cascades;            use Embeddings_and_Cascades;
with Standard_Witness_Solutions;

procedure ts_witsols is

-- DESCRIPTION :
--   Test on the package to store witness solutions.

  procedure Main is

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
    Standard_Witness_Solutions.Initialize(nv,topdim);
  end Main;

begin
  Main;
end ts_witsols;
