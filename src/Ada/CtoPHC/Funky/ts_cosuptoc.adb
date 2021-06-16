with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Symbol_Table;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Integer_Arrays_io;               use C_Integer_Arrays_io;
with C_Double_Arrays;                   use C_Double_Arrays;
with C_Double_Arrays_io;                use C_Double_Arrays_io;
with Coefficient_Support_Polynomials;   use Coefficient_Support_Polynomials;
with Coefficient_Support_Poly_Systems;  use Coefficient_Support_Poly_Systems;

procedure ts_cosuptoc is

  function cosupoly_to_c ( n,ns : integer32; s : C_Integer_Array;
                             nc : integer32; c : C_Double_Array )
                         return integer32;
  pragma Import(C,cosupoly_to_c,"cosupoly_to_c");

  -- DESCRIPTION :
  --   The function passes the coefficient-support representation of a
  --   polynomial in n variables from Ada to C.

  function cosupsys_to_c ( n,m : integer32; moncnt : C_Integer_Array;
                           ns : integer32; s : C_Integer_Array;
                           nc : integer32; c : C_Double_Array )
                         return integer32;
  pragma Import(C,cosupsys_to_c,"cosupsys_to_c");
 
  -- DESCRIPTION :
  --   The function passes the coefficient-support representation of a
  --   system of n equations in m variables from Ada to C.

  procedure Export_to_C ( n : in natural32; p : in Poly ) is

  -- DESCRIPTION :
  --   Converts the polynomial to a representation for export to C.

    sup : constant C_Integer_Array := Support(p);
    cff : constant C_Double_Array := Coefficients(p);
    return_of_call : integer32;

  begin
    put("The support :"); put(sup'length,sup); new_line;
    put_line("The coefficients : "); put(cff'length,cff);
    return_of_call
      := cosupoly_to_c(integer32(n),sup'length,sup,cff'length,cff);
    if return_of_call = 0 then
      put_line("Call to C terminated normally.");
    else
      put("Call to C terminated abnormally with exit code ");
      put(return_of_call,1); new_line;
    end if;
  end Export_to_C;

  procedure Export_to_C ( p : in Poly_Sys ) is

    n : constant natural32 := natural32(p'length);
    m : constant natural32 := Number_of_Unknowns(p(p'first));
    moncnt : constant C_Integer_Array := Monomial_Count(p);
    monsum : constant natural32 := natural32(Sum(moncnt));
    sup : constant C_Integer_Array := Support(n,monsum,moncnt,p);
    cff : constant C_Double_Array := Coefficients(monsum,moncnt,p);
    return_of_call : integer32;

  begin
    put("The number of equations : "); put(n,1); put_line(".");
    put("The number of variables : "); put(m,1); put_line(".");
    put("#monomials :"); put(moncnt'length,moncnt); put(" with sum : ");
    put(monsum,1); new_line;
    put("The support :"); put(sup'length,sup); new_line;
    put_line("The coefficients "); put(cff'length,cff);
    return_of_call
      := cosupsys_to_c(integer32(n),integer32(m),moncnt,
                       sup'length,sup,cff'length,cff);
    if return_of_call = 0 then
      put_line("Call to C terminated normally.");
    else
      put("Call to C terminated abnormally with exit code ");
      put(return_of_call,1); new_line;
    end if;
  end Export_to_C;

  procedure Main is

    n : natural32 := 0;
    p : Poly;
    ans : character;
    lp : Link_to_Poly_Sys;

  begin
    put_line("Choose one of the following :");
    put_line("  1. export a single complex multivariate polynomial to C; or");
    put_line("  2. export a polynomial system from Ada to C.");
    put("Type 1 or 2 to make your choice : ");
    Ask_Alternative(ans,"12");
    if ans = '1' then
      put("Give the number of variables : "); get(n);
      Symbol_Table.Init(n);
      put("Give a polynomial terminate with semicolon : "); get(p);
      put("Your polynomial is "); put(p); new_line;
      Export_to_C(n,p);
    else
      new_line;
      get(lp);
      Export_to_C(lp.all);
    end if;
  end Main;

begin
  new_line;
  put_line("Passing polynomial and polynomial systems from Ada to C.");
  new_line;
  Main;
end ts_cosuptoc;
