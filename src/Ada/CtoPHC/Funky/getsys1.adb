with text_io,integer_io;                use text_io,integer_io;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Integer_Arrays_io;               use C_Integer_Arrays_io;
with C_Double_Arrays;                   use C_Double_Arrays;
with C_Double_Arrays_io;                use C_Double_Arrays_io;
with Coefficient_Support_Poly_Systems;  use Coefficient_Support_Poly_Systems;

procedure getsys1 is

-- DESCRIPTION :
--   Calls the PHCpack routines to read a polynomial system.

  function getsys2 ( n,m : integer; moncnt : C_Integer_Array;
                     ns : integer; s : C_Integer_Array;
                     nc : integer; c : C_Double_Array ) return integer;
  pragma Import(C,getsys2,"getsys2");

  procedure Export_to_C ( p : in Poly_Sys ) is

    n : constant natural := p'length;
    m : constant natural := natural(Number_of_Unknowns(p(p'first)));
    moncnt : constant C_Integer_Array := Monomial_Count(p);
    monsum : constant natural := Sum(moncnt);
    sup : constant C_Integer_Array := Support(n,monsum,moncnt,p);
    cff : constant C_Double_Array := Coefficients(monsum,moncnt,p);
    return_of_call : integer;

  begin
    put("The number of equations : "); put(n,1); put_line(".");
    put("The number of variables : "); put(m,1); put_line(".");
    put("#monomials :"); put(moncnt'length,moncnt); put(" with sum : ");
    put(monsum,1); new_line;
    put("The support :"); put(sup'length,sup); new_line;
    put_line("The coefficients "); put(cff'length,cff);
    return_of_call := getsys2(n,m,moncnt,sup'length,sup,cff'length,cff);
    if return_of_call = 0
     then put_line("Call to C terminated normally.");
     else put("Call to C terminated abnormally with exit code ");
          put(return_of_call,1); new_line;
    end if;
  end Export_to_C;

  procedure Main is

    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Reading a polynomial system...");
    get(lp);
    new_line;
    put_line("The polynomial system read:"); 
    put(lp.all);
    new_line;
    put_line("Exporting the polynomial system to C...");
    Export_to_C(lp.all);
  end Main;

begin
  Main;
end getsys1;
