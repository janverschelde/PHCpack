with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Dense_Series;             use Standard_Dense_Series;
with Standard_Dense_Series_io;          use Standard_Dense_Series_io;
with Standard_Series_Polynomials;
with Series_and_Polynomials;            use Series_and_Polynomials;

procedure ts_serpol is

-- DESCRIPTION :
--   Tests the development of polynomials in several variables,
--   with truncated power series as coefficients.
 
  procedure Write ( s : in Standard_Series_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Very simple output of a polynomial with series coefficients.
 
    cnt : natural32 := 0;

    procedure Visit_Term ( t : in Standard_Series_Polynomials.Term;
                           c : out boolean ) is
   
      cf : Series := t.cf;

    begin
      cnt := cnt + 1;
      put("The coefficient of term "); put(cnt); put_line(" :");
      put(cf);
      put("has order "); put(cf.order,1);
      put(" and degrees : "); put(t.dg.all); new_line;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Series_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
  end Write;

  procedure Main is 

  -- DESCRIPTION :
  --   Prompts the user for the number of variables.
  --   and reads in a regular polynomial in several variables,
  --   for conversion into a series polynomial.

    n : natural32 := 0;
    p,q : Standard_Complex_Polynomials.Poly;
    s : Standard_Series_Polynomials.Poly;
    ans : character;
 
  begin
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.init(n);
    put("Give a polynomial : "); get(p);
    put("> your polynomial : "); put(p); new_line;
    new_line;
    put("Extra output during the conversion ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then s := Polynomial_to_Series_Polynomial(p,true);
     else s := Polynomial_to_Series_Polynomial(p);
    end if;
    new_line;
    put_line("The series polynomial s :");
    Write(s);
    q := Series_Polynomial_to_Polynomial(s,ans = 'y');
    put("s as poly : "); put(q); new_line;
  end Main;

begin
  Main;
end ts_serpol;
