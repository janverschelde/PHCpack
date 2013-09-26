with text_io;                           use text_io;
with Symbol_Table;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Multprec_Complex_Polynomials;      use Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;   use Multprec_COmplex_Polynomials_io;
with Multprec_Scaling;                  use Multprec_Scaling;

procedure ts_scal is

  procedure Main is

    n : natural32 := 0;
    p : Poly;

  begin
    put("Give the number of unknowns : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial, terminate with semi-colon...");
    get(p);
    put_line("Your polynomial : "); put(p);
    Scale(p);
    put_line("The polynomial after scaling : "); put(p);
  end Main;

begin
  Main;
end ts_scal;
