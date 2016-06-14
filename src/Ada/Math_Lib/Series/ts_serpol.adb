with text_io;                           use text_io;
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
with Standard_Series_Polynomials;

procedure ts_serpol is

-- DESCRIPTION :
--   Tests the development of polynomials in several variables,
--   with truncated power series as coefficients.

  function Polynomial_to_Series_Polynomial
             ( p : Standard_Complex_Polynomials.Poly )
             return Standard_Series_Polynomials.Poly is

    res : Standard_Series_Polynomials.Poly
        := Standard_Series_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in Standard_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Adds the information in the term of a multivariate polynomial
    --   as a term of a series polynomial.  The first variable in t
    --   is considered as the variable in the truncated power series.
    --   The other variables are copied to variable of the same index
    --   minus one in the term of the series polynomial.

      rtm : Standard_Series_Polynomials.Term;
      ord : constant integer32 := integer32(t.dg(t.dg'first));
      dim : constant integer32 := t.dg'last-1;
      rcf : Series := Create(0.0,ord);
      cnt : natural32 := 0;

    begin
      rcf.cff(ord) := t.cf;
      rtm.cf := rcf;
      rtm.dg := new Standard_Natural_Vectors.Vector(t.dg'first..dim);
      for i in rtm.dg'range loop
        rtm.dg(i) := t.dg(i+1);
      end loop;
      cnt := cnt + 1;
      put("Adding term "); put(cnt,1); put_line(" with coefficient :");
      put(rtm.cf);
      put("order : "); put(ord,1);
      put(" and degrees : "); put(rtm.dg.all); new_line;
      Standard_Series_Polynomials.Add(res,rtm);
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series_Polynomial;

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
    p : Standard_Complex_Polynomials.Poly;
    s : Standard_Series_Polynomials.Poly;
 
  begin
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.init(n);
    put("Give a polynomial : "); get(p);
    put("> your polynomial : "); put(p); new_line;
    s := Polynomial_to_Series_Polynomial(p);
    put_line("The series polynomial :");
    Write(s);
  end Main;

begin
  Main;
end ts_serpol;
