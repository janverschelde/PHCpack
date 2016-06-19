with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
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

  function Set_Order ( i : integer32;
                       d : Standard_Complex_Polynomials.Degrees )
                     return integer32 is

  -- DESCRIPTION :
  --   Returns the order of the coefficient series,
  --   which is either zero if i = 0, or otherwise d(i).

  begin
    if i = 0
     then return 0;
     else return integer32(d(i));
    end if;
  end Set_Order;

  function Set_Dimension ( i : integer32;
                           d : Standard_Complex_Polynomials.Degrees )
                         return integer32 is

  -- DESCRIPTION :
  --   Returns the number of variables of the series polynomial
  --   which is either d'last if i = 0, or otherwise d'last-1.

  begin
    if i = 0
     then return d'last;
     else return d'last-1;
    end if;
  end Set_Dimension;

  function Polynomial_to_Series_Polynomial
             ( p : Standard_Complex_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return Standard_Series_Polynomials.Poly is

    use Standard_Dense_Series;

    res : Standard_Series_Polynomials.Poly
        := Standard_Series_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in Standard_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Adds the information in the term of a multivariate polynomial
    --   as a term of a series polynomial.  If idx = 0, then the order of
    --   the series coefficient is zero and the coefficient is copied.
    --   Otherwise, the variable with index idx in t is taken as the
    --   variable in the truncated power series.  The other variables
    --   in t are relocated to fit a polynomial with one variable less.

      rtm : Standard_Series_Polynomials.Term;
      ord : constant integer32 := Set_Order(idx,t.dg);
      dim : constant integer32 := Set_Dimension(idx,t.dg);
      rcf : Series := Create(0.0,ord);
      cnt : natural32 := 0;

    begin
      rcf.cff(ord) := t.cf;
      rtm.cf := rcf;
      rtm.dg := new Standard_Natural_Vectors.Vector(t.dg'first..dim);
      if idx = 0 then
        for i in rtm.dg'range loop
          rtm.dg(i) := t.dg(i);
        end loop;
      else
        for i in 1..(idx-1) loop
          rtm.dg(i) := t.dg(i);
        end loop;
        for i in (idx+1)..t.dg'last loop
          rtm.dg(i-1) := t.dg(i);
        end loop;
      end if;
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
               idx : integer32 := 0; verbose : boolean := false )
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
      dim1 : integer32;

    begin
      if verbose then
        put("term with degrees :"); put(t.dg.all);
        put(" has series of order "); put(cffs.order,1); new_line;
        put_line("the series : "); put(cffs);
      end if;
      if idx = 0 then
        declare
          rt : Standard_Complex_Polynomials.Term;
        begin
          rt.cf := cffs.cff(0);
          rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
          Standard_Complex_Polynomials.Add(res,rt);
          Standard_Complex_Polynomials.Clear(rt);
        end;
      else -- idx > 0
        dim1 := t.dg'last+1;
        for k in 0..cffs.order loop
          rtpc := cffs.cff(k);
          if not Equal(rtpc,zero) then
            declare
              rt : Standard_Complex_Polynomials.Term;
            begin
              rt.cf := rtpc;
              rt.dg := new Standard_Natural_Vectors.Vector(1..dim1);
              for i in 1..(idx-1) loop
                rt.dg(i) := t.dg(i);
              end loop;
              rt.dg(idx) := natural32(k);
              for i in (idx+1)..rt.dg'last loop
                rt.dg(i) := t.dg(i-1);
              end loop;
              Standard_Complex_Polynomials.Add(res,rt);
              Standard_Complex_Polynomials.Clear(rt);
            end;
          end if;
        end loop;
      end if;
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
               idx : integer32 := 0; verbose : boolean := false )
             return Standard_Series_Poly_Systems.Poly_Sys is

    res : Standard_Series_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Polynomial_to_Series_Polynomial(p(i),idx,verbose);
    end loop;
    return res;
  end System_to_Series_System;

  function Series_System_to_System
             ( s : Standard_Series_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(s'range);

  begin
    for i in s'range loop
      res(i) := Series_Polynomial_to_Polynomial(s(i),idx,verbose);
    end loop;
    return res;
  end Series_System_to_System;

  procedure Set_Order ( p : in out Standard_Series_Polynomials.Poly;
                        order : in integer32 ) is

    procedure Change_Order ( t : in out Standard_Series_Polynomials.Term;
                             c : out boolean ) is
    begin
      if t.cf.order < order then
        for k in t.cf.order+1..order loop
          t.cf.cff(k) := Create(0.0);
        end loop;
      end if;
      t.cf.order := order;
      c := true;
    end Change_Order;
    procedure Change_Orders is
      new Standard_Series_Polynomials.Changing_Iterator(Change_Order);

  begin
    Change_Orders(p);
  end Set_Order;

  procedure Set_Order ( p : in out Standard_Series_Poly_Systems.Poly_Sys;
                        order : in integer32 ) is
  begin
    for i in p'range loop
      Set_Order(p(i),order);
    end loop;
  end Set_Order;
 
end Series_and_Polynomials;
