with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Dense_Series;             use Standard_Dense_Series;
with Standard_Dense_Series_io;          use Standard_Dense_Series_io;
with Standard_Dense_Series_Vectors;
with Standard_Random_Series;            use Standard_Random_Series;
with Standard_Series_Polynomials;
with Standard_Series_Poly_Functions;
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

  procedure Test_Conversion is 

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
  end Test_Conversion;

  function Factor ( n,k : integer32; s : Series )
                  return Standard_Series_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns x[k] - s as a polynomial in n variables.
  --   All coefficients of the polynomial on return have the same order,
  --   the same as s.order.

  -- REQUIRED : k is in range 1..n.

    res : Standard_Series_Polynomials.Poly;
    one : constant Series := Create(1.0,s.order);
    trm : Standard_Series_Polynomials.Term;

  begin
    trm.cf := one;
    trm.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    trm.dg(k) := 1;
    res := Standard_Series_Polynomials.Create(trm);
    trm.cf := s;
    trm.dg(k) := 0;
    Standard_Series_Polynomials.Sub(res,trm);
    return res;
  end Factor;

  function Product ( s : Standard_Dense_Series_Vectors.Vector )
                   return Standard_Series_Polynomials.Poly is

  -- DESCRIPION :
  --   Returns the product of the factors x[k] - s[k], for k in s'range,
  --   where s'first = 1.

    dim : constant integer32 := s'last;
    res : Standard_Series_Polynomials.Poly := Factor(dim,1,s(1));
    fac : Standard_Series_Polynomials.Poly;

  begin
    for k in 2..s'last loop
      fac := Factor(dim,k,s(k));
      Standard_Series_Polynomials.Mul(res,fac);
      Standard_Series_Polynomials.Clear(fac);
    end loop;
    return res;
  end Product;

  procedure Test_Evaluation is

  -- DESCRIPTION :
  --   Prompts for the number of variables and the order of the series.
  --   Then as many random series as the number of variables are generated.
  --   The polynomial is of the product of x[k] - s[k], where k ranges
  --   over the number of variables 'x' and series 's'.
  --   So the evaluation at the series should produce zero.

    order,dim : integer32 := 0;

  begin
    new_line;
    put("Give the number of variables in the polynomial : "); get(dim);
    put("Give the order of the power series : "); get(order);
    declare
      rns : constant Standard_Dense_Series_Vectors.Vector(1..dim)
          := Random_Series_Vector(1,dim,order);
      pol : Standard_Series_Polynomials.Poly := Product(rns);
      eva : constant Series := Standard_Series_Poly_Functions.Eval(pol,rns);
    begin
      for i in 1..dim loop
        put("random series "); put(i,1); put_line(" :");
        put(rns(i));
      end loop;
      put_line("The polynomial :"); Write(pol);
      put_line("The value at the polynomial :"); put(eva);
    end;
  end Test_Evaluation;

  
  procedure Main is

  -- DESCRIPTION :
  --   Displays the menu of tests
  --   and prompts the user to make a choice.

    ans : character;

  begin
    new_line;
    put_line("MENU to test polynomials with series coefficients :");
    put_line("  0. test conversion from/to ordinary polynomials;");
    put_line("  1. test evaluation at power series;");
    put("Type 0 or 1 to select a test : ");
    Ask_Alternative(ans,"01");
    case ans is
      when '0' => Test_Conversion;
      when '1' => Test_Evaluation;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_serpol;
