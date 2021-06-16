with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Interfaces.C;                      use Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with C_Integer_io;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_io;
with C_Double_Arrays;                   use C_Double_Arrays;
with Symbol_Table;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Coefficient_Support_Polynomials;   use Coefficient_Support_Polynomials;
with Coefficient_Support_Poly_Systems;  use Coefficient_Support_Poly_Systems;

procedure ts_cosup is

-- DESCRIPTION :
--   This test program develops a basic interface to C for polynomials
--   in several variables with complex coefficients.  Basically we split
--   the representation of a polynomial in two C style arrays: one for
--   the coefficients and one for the supports.

  procedure put ( a : in C_Integer_Array ) is

  -- DESCRIPTION :
  --   Writes the array of C integers to standard output.

  begin
    for i in a'range loop
      put(" "); C_integer_io.put(a(i),1);
    end loop;
  end put;

  procedure put ( a : in C_Double_Array ) is

  -- DESCRIPTION :
  --   Writes the array of C doubles to standard output.

  begin
    for i in a'range loop
      C_double_io.put(a(i)); new_line;
    end loop;
  end put;

  procedure Test_Poly_Creation ( n,m : in natural32 ) is

    numexp : constant size_T := size_T(n*m-1);
    numcff : constant size_T := size_T(2*m-1);
    sup : C_Integer_Array(0..numexp);
    cff : C_Double_Array(0..numcff);
    ind : Interfaces.C.size_T := 0;
    p : Poly;

  begin
    put("Reading "); put(m,1);
    put_line(" integer vectors for the support...");
    for i in 1..m loop
      put("  monomial "); put(i,1); put(" : ");
      for j in 1..n loop
        C_integer_io.get(sup(ind));
        ind := ind+1;
      end loop;
    end loop;
    put("The stacked support vector : ");
    put(sup); new_line;
    ind := 0;
    put("Reading "); put(2*m,1); 
    put_line(" doubles as real and imaginary parts of the coefficients...");
    for i in 1..m loop
      put("  coefficient "); put(i,1); put(" : ");
      C_double_io.get(cff(ind));
      ind := ind + 1;
      C_double_io.get(cff(ind));
      ind := ind + 1;
    end loop;
    put_line("The coefficients : "); put(cff);
    p := Create(n,cff,sup);
    put("The polynomial : "); put_line(p); new_line;
    put_line("Its support : ");
    put(Support(p)); new_line;
    put_line("Its coefficients : ");
    put(Coefficients(p));
  end Test_Poly_Creation;

  procedure Test_Cosup_to_Poly is

    n,m : natural32 := 0;

  begin
    new_line;
    put_line("Turning coefficient and support into a polynomial.");
    new_line;
    put("Give the number of variables : "); get(n);
    put("Give the number of monomials : "); get(m);
    Test_Poly_Creation(n,m);
  end Test_Cosup_to_Poly;

  procedure Test_Cosup_to_Poly_Sys is

    n,m : natural32 := 0;

  begin
    new_line;
    put_line("Turning a coefficient-support representation into a system.");
    new_line;
    put("Give the number of variables : "); get(n);
    put("Give the number of equations : "); get(m);
  end Test_Cosup_to_Poly_Sys;

  procedure Test_Poly_to_Cosup is

    p,q : Poly;
    n : natural32 := 0;

  begin
    new_line;
    put_line("The coefficient-support representation of a polynomial.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Type in a polynomial, terminate by semi-colon :");
    get(p);
    put("Your polynomial : "); put(p); new_line;
    declare
      sup : constant C_Integer_Array := Support(p);
      cff : constant C_Double_Array := Coefficients(p);
    begin
      put("The support vector : "); put(sup); new_line;
      put_line("The coefficients : "); put(cff);
      q := Create(n,cff,sup);
    end;
    put_line("The polynomial created from coefficients and support : ");
    put_line(q);
  end Test_Poly_to_Cosup;

  procedure Test_Poly_Sys_to_Cosup is

    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("The coefficient-support representation of a polynomial system.");
    new_line;
    get(lp);
    declare
      n : constant natural32 := Number_of_Unknowns(lp(lp'first));
      mon : constant C_Integer_Array := Monomial_Count(lp.all);
      moncnt : constant natural32 := natural32(Sum(mon));
      sup : constant C_Integer_Array := Support(n,moncnt,mon,lp.all);
    begin
      put("The number of monomials : "); put(mon); 
      put(" and the sum is "); put(moncnt,1); new_line;
      put("The support : "); put(sup); new_line;
    end;
  end Test_Poly_Sys_to_Cosup;

  procedure Write_Polynomial_System ( x : in C_Double_Array ) is

  -- DESCRIPTION :
  --   Writes the polynomial system from its concatenated 
  --   coefficient support representation.

    n : constant natural32 := Dimension(x);
    m : constant C_Integer_Array := Monomial_Count(x);
    s : constant C_Integer_Array := Support(x);
    c : constant C_Double_Array := Coefficients(x);
    p : constant Poly_Sys := Create(n,m,c,s);

  begin
    put(p);
  end Write_Polynomial_System;

  procedure Test_Concatenation is

    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("The coefficient-support representation of a polynomial system.");
    new_line;
    get(lp);
    declare
      n : constant natural32 := Number_of_Unknowns(lp(lp'first));
      mon : constant C_Integer_Array := Monomial_Count(lp.all);
      moncnt : constant natural32 := natural32(Sum(mon));
      sup : constant C_Integer_Array := Support(n,moncnt,mon,lp.all);
      cff : constant C_Double_Array := Coefficients(moncnt,mon,lp.all);
      cct : constant C_Double_Array := Concat(n,mon,cff,sup);
    begin
      Write_Polynomial_System(cct);
    end;
  end Test_Concatenation;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing a coefficient-support representation of polynomials.");
    new_line;
    put_line("MENU to convert data representations : ");
    put_line("  1. from coefficient-support representation to polynomial;");
    put_line("  2. from polynomial to coefficient-support representation.");
    put_line("  3. from coefficient-support representation to system;");
    put_line("  4. from system to coefficient-support representation.");
    put_line("  5. concatenation of coefficient-support representation.");
    put("Type 1, 2, 3, 4, or 5 to select a conversion : ");
    Ask_Alternative(ans,"12345");
    case ans is
      when '1' => Test_Cosup_to_Poly;
      when '2' => Test_Poly_to_Cosup;
      when '3' => Test_Cosup_to_Poly_Sys;
      when '4' => Test_Poly_Sys_to_Cosup;
      when '5' => Test_Concatenation;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_cosup;
