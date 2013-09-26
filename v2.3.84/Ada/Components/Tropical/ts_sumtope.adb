with text_io;                           use text_io;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;

procedure ts_sumtope is

-- DESCRIPTION :
--   Writes the matrix that contains the points that span the
--   Minkowski sum of the Newton polytopes in a given system.

  function Product ( p : Poly_Sys ) return Poly is

  -- DESCRIPTION :
  --   Returns the product of the polynomials in p.

    res : Poly;

  begin
    Copy(p(p'first),res);
    for i in p'first+1..p'last loop
      Mul(res,p(i));
    end loop;
    return res;
  end Product;

  procedure Main is

    p : Link_to_Poly_Sys;
    q : Poly;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(p);
    q := Product(p.all);
    put_line("The product of the polynomials in the system :");
    put(q);
  end Main;

begin
  Main;
end ts_sumtope;
