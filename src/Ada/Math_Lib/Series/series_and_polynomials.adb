with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Dense_Series_io;          use Standard_Dense_Series_io;

package body Series_and_Polynomials is

  function Series_to_Polynomial
             ( s : Standard_Dense_Series.Series )
             return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;
    zero : constant Complex_Number := Create(0.0);

    procedure Add_Term ( k : in integer32; c : in Complex_Number ) is

    -- DESCRIPTION :
    --   Adds c*t^k to the polynomial res.

      t : Standard_Complex_Polynomials.Term;

    begin
      t.cf := c;
      t.dg := new Standard_Natural_Vectors.Vector(1..1);
      t.dg(1) := natural32(k);
      Standard_Complex_Polynomials.Add(res,t);
      Standard_Complex_Polynomials.Clear(t);
    end Add_Term;

  begin
    for i in 0..s.order loop
      if not Equal(s.cff(i),zero)
       then Add_Term(i,s.cff(i));
      end if;
    end loop;
    return res;
  end Series_to_Polynomial;

  function Polynomial_to_Series
             ( p : Standard_Complex_Polynomials.Poly )
             return Standard_Dense_Series.Series is

    res : Standard_Dense_Series.Series;

    procedure Visit_Term ( t : in Standard_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Assigns t.cf to the res.cff(t.dg(1))
    --   and updates the order if needed.

      d : constant integer32 := integer32(t.dg(1));

    begin
      if d > res.order then
        for k in res.order+1..d loop -- fill in with zeroes
          res.cff(k) := Create(0.0);
        end loop;
        res.order := d;
      end if;
      res.cff(d) := t.cf;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    res.order := 0;
    res.cff(0) := Create(0.0);
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series;

  function Polynomial_to_Series_Polynomial
             ( p : Standard_Complex_Polynomials.Poly;
               verbose : boolean := false )
             return Standard_Series_Polynomials.Poly is

    use Standard_Dense_Series;

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
      if verbose then
        put("Adding term "); put(cnt,1); put_line(" with coefficient :");
        put(rtm.cf);
        put("order : "); put(ord,1);
        put(" and degrees : "); put(rtm.dg.all); new_line;
      end if;
      Standard_Series_Polynomials.Add(res,rtm);
      Standard_Series_Polynomials.Clear(rtm);
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series_Polynomial;

  function Series_Polynomial_to_Polynomial
             ( s : Standard_Series_Polynomials.Poly;
               verbose : boolean := false )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Dense_Series;

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in Standard_Series_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Every term in the coefficient of the series in t.cf
    --   with a nonzero coefficient will contribute one monomial
    --   to the result.

      cffs : Series := t.cf;
      zero : constant Complex_Number := Create(0.0);
      rtpc : Complex_Number;
      dim1 : constant integer32 := t.dg'last+1;

    begin
      if verbose then
        put("term with degrees :"); put(t.dg.all);
        put(" has series of order "); put(cffs.order,1); new_line;
        put_line("the series : "); put(cffs);
      end if;
      for k in 0..cffs.order loop
        rtpc := cffs.cff(k);
        if not Equal(rtpc,zero) then
          declare
            rt : Standard_Complex_Polynomials.Term;
          begin
            rt.cf := rtpc;
            rt.dg := new Standard_Natural_Vectors.Vector(1..dim1);
            rt.dg(1) := natural32(k);
            for i in t.dg'range loop
              rt.dg(i+1) := t.dg(i);
            end loop;
            Standard_Complex_Polynomials.Add(res,rt);
            Standard_Complex_Polynomials.Clear(rt);
          end;
        end if;
      end loop;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Series_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
    return res;
  end Series_Polynomial_to_Polynomial;

  function System_to_Series_System
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               verbose : boolean := false )
             return Standard_Series_Poly_Systems.Poly_Sys is

    res : Standard_Series_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Polynomial_to_Series_Polynomial(p(i),verbose);
    end loop;
    return res;
  end System_to_Series_System;

  function Series_System_to_System
             ( s : Standard_Series_Poly_Systems.Poly_Sys;
               verbose : boolean := false )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(s'range);

  begin
    for i in s'range loop
      res(i) := Series_Polynomial_to_Polynomial(s(i),verbose);
    end loop;
    return res;
  end Series_System_to_System;
 
end Series_and_Polynomials;
