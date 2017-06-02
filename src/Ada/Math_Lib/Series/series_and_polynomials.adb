with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Dense_Series_io;          use Standard_Dense_Series_io;
with DoblDobl_Dense_Series_io;          use DoblDobl_Dense_Series_io;
with QuadDobl_Dense_Series_io;          use QuadDobl_Dense_Series_io;

package body Series_and_Polynomials is

  function Series_to_Polynomial
             ( s : Standard_Dense_Series.Series )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Numbers;

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
      Standard_Complex_Polynomials.Clear(t.dg);
    end Add_Term;

  begin
    for i in 0..s.deg loop
      if not Equal(s.cff(i),zero)
       then Add_Term(i,s.cff(i));
      end if;
    end loop;
    return res;
  end Series_to_Polynomial;

  function Series_to_Polynomial
             ( s : DoblDobl_Dense_Series.Series )
             return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Numbers;

    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Null_Poly;
    zero : constant Complex_Number := Create(integer(0));

    procedure Add_Term ( k : in integer32; c : in Complex_Number ) is

    -- DESCRIPTION :
    --   Adds c*t^k to the polynomial res.

      t : DoblDobl_Complex_Polynomials.Term;

    begin
      t.cf := c;
      t.dg := new Standard_Natural_Vectors.Vector(1..1);
      t.dg(1) := natural32(k);
      DoblDobl_Complex_Polynomials.Add(res,t);
      DoblDobl_Complex_Polynomials.Clear(t.dg);
    end Add_Term;

  begin
    for i in 0..s.deg loop
      if not Equal(s.cff(i),zero)
       then Add_Term(i,s.cff(i));
      end if;
    end loop;
    return res;
  end Series_to_Polynomial;

  function Series_to_Polynomial
             ( s : QuadDobl_Dense_Series.Series )
             return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Numbers;

    res : QuadDobl_Complex_Polynomials.Poly
        := QuadDobl_Complex_Polynomials.Null_Poly;
    zero : constant Complex_Number := Create(integer(0));

    procedure Add_Term ( k : in integer32; c : in Complex_Number ) is

    -- DESCRIPTION :
    --   Adds c*t^k to the polynomial res.

      t : QuadDobl_Complex_Polynomials.Term;

    begin
      t.cf := c;
      t.dg := new Standard_Natural_Vectors.Vector(1..1);
      t.dg(1) := natural32(k);
      QuadDobl_Complex_Polynomials.Add(res,t);
      QuadDobl_Complex_Polynomials.Clear(t.dg);
    end Add_Term;

  begin
    for i in 0..s.deg loop
      if not Equal(s.cff(i),zero)
       then Add_Term(i,s.cff(i));
      end if;
    end loop;
    return res;
  end Series_to_Polynomial;

  function Polynomial_to_Series
             ( p : Standard_Complex_Polynomials.Poly;
               idx : integer32 := 1 )
             return Standard_Dense_Series.Series is

    use Standard_Complex_Numbers;

    res : Standard_Dense_Series.Series;

    procedure Visit_Term ( t : in Standard_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Assigns t.cf to the res.cff(t.dg(idx))
    --   and updates the degree if needed.

      d : constant integer32 := integer32(t.dg(idx));

    begin
      if d > res.deg then
        for k in res.deg+1..d loop -- fill in with zeroes
          res.cff(k) := Create(0.0);
        end loop;
        res.deg := d;
      end if;
      res.cff(d) := t.cf;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    res.deg := 0;
    res.cff(0) := Create(0.0);
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series;

  function Polynomial_to_Series
             ( p : DoblDobl_Complex_Polynomials.Poly;
               idx : integer32 := 1 )
             return DoblDobl_Dense_Series.Series is

    use DoblDobl_Complex_Numbers;

    res : DoblDobl_Dense_Series.Series;

    procedure Visit_Term ( t : in DoblDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Assigns t.cf to the res.cff(t.dg(idx))
    --   and updates the degree if needed.

      d : constant integer32 := integer32(t.dg(idx));

    begin
      if d > res.deg then
        for k in res.deg+1..d loop -- fill in with zeroes
          res.cff(k) := Create(integer(0));
        end loop;
        res.deg := d;
      end if;
      res.cff(d) := t.cf;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new DoblDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    res.deg := 0;
    res.cff(0) := Create(integer(0));
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series;

  function Polynomial_to_Series
             ( p : QuadDobl_Complex_Polynomials.Poly;
               idx : integer32 := 1 )
             return QuadDobl_Dense_Series.Series is

    use QuadDobl_Complex_Numbers;

    res : QuadDobl_Dense_Series.Series;

    procedure Visit_Term ( t : in QuadDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Assigns t.cf to the res.cff(t.dg(idx))
    --   and updates the degree if needed.

      d : constant integer32 := integer32(t.dg(idx));

    begin
      if d > res.deg then
        for k in res.deg+1..d loop -- fill in with zeroes
          res.cff(k) := Create(integer(0));
        end loop;
        res.deg := d;
      end if;
      res.cff(d) := t.cf;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new QuadDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    res.deg := 0;
    res.cff(0) := Create(integer(0));
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series;

  function Series_Vector_to_System
             ( v : Standard_Dense_Series_Vectors.Vector )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(v'range);

  begin
    for k in v'range loop
      res(k) := Series_to_Polynomial(v(k));
    end loop;
    return res;
  end Series_Vector_to_System;

  function Series_Vector_to_System
             ( v : DoblDobl_Dense_Series_Vectors.Vector )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(v'range);

  begin
    for k in v'range loop
      res(k) := Series_to_Polynomial(v(k));
    end loop;
    return res;
  end Series_Vector_to_System;

  function Series_Vector_to_System
             ( v : QuadDobl_Dense_Series_Vectors.Vector )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(v'range);

  begin
    for k in v'range loop
      res(k) := Series_to_Polynomial(v(k));
    end loop;
    return res;
  end Series_Vector_to_System;

  function System_to_Series_Vector
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 1 )
             return Standard_Dense_Series_Vectors.Vector is

    res : Standard_Dense_Series_Vectors.Vector(p'range);

  begin
    for k in p'range loop
      res(k) := Polynomial_to_Series(p(k),idx);
    end loop;
    return res;
  end System_to_Series_Vector;

  function System_to_Series_Vector
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 1 )
             return DoblDobl_Dense_Series_Vectors.Vector is

    res : DoblDobl_Dense_Series_Vectors.Vector(p'range);

  begin
    for k in p'range loop
      res(k) := Polynomial_to_Series(p(k),idx);
    end loop;
    return res;
  end System_to_Series_Vector;

  function System_to_Series_Vector
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 1 )
             return QuadDobl_Dense_Series_Vectors.Vector is

    res : QuadDobl_Dense_Series_Vectors.Vector(p'range);

  begin
    for k in p'range loop
      res(k) := Polynomial_to_Series(p(k),idx);
    end loop;
    return res;
  end System_to_Series_Vector;

  function Series_VecVec_to_System_Array
             ( v : Standard_Dense_Series_VecVecs.VecVec )
             return Standard_Complex_Poly_Systems.Array_of_Poly_Sys is

    res : Standard_Complex_Poly_Systems.Array_of_Poly_Sys(v'range);

    use Standard_Dense_Series_Vectors;

  begin
    for k in v'range loop
      if v(k) /= null then
        declare
          p : constant Standard_Complex_Poly_Systems.Poly_Sys
            := Series_Vector_to_System(v(k).all);
        begin
          res(k) := new Standard_Complex_Poly_Systems.Poly_Sys'(p);
        end;
      end if;
    end loop;
    return res;
  end Series_VecVec_to_System_Array;

  function Series_VecVec_to_System_Array
             ( v : DoblDobl_Dense_Series_VecVecs.VecVec )
             return DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys(v'range);

    use DoblDobl_Dense_Series_Vectors;

  begin
    for k in v'range loop
      if v(k) /= null then
        declare
          p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
            := Series_Vector_to_System(v(k).all);
        begin
          res(k) := new DoblDobl_Complex_Poly_Systems.Poly_Sys'(p);
        end;
      end if;
    end loop;
    return res;
  end Series_VecVec_to_System_Array;

  function Series_VecVec_to_System_Array
             ( v : QuadDobl_Dense_Series_VecVecs.VecVec )
             return QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys(v'range);

    use QuadDobl_Dense_Series_Vectors;

  begin
    for k in v'range loop
      if v(k) /= null then
        declare
          p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
            := Series_Vector_to_System(v(k).all);
        begin
          res(k) := new QuadDobl_Complex_Poly_Systems.Poly_Sys'(p);
        end;
      end if;
    end loop;
    return res;
  end Series_VecVec_to_System_Array;

  function System_Array_to_Series_VecVec
             ( p : Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
               idx : integer32 := 1 )
             return Standard_Dense_Series_VecVecs.VecVec is

    res : Standard_Dense_Series_VecVecs.VecVec(p'range);

    use Standard_Complex_Poly_Systems;

  begin
    for k in p'range loop
      if p(k) /= null then
        declare
          v : constant Standard_Dense_Series_Vectors.Vector
            := System_to_Series_Vector(p(k).all,idx);
        begin
          res(k) := new Standard_Dense_Series_Vectors.Vector'(v);
        end;
      end if;
    end loop;
    return res;
  end System_Array_to_Series_VecVec;

  function System_Array_to_Series_VecVec
             ( p : DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
               idx : integer32 := 1 )
             return DoblDobl_Dense_Series_VecVecs.VecVec is

    res : DoblDobl_Dense_Series_VecVecs.VecVec(p'range);

    use DoblDobl_Complex_Poly_Systems;

  begin
    for k in p'range loop
      if p(k) /= null then
        declare
          v : constant DoblDobl_Dense_Series_Vectors.Vector
            := System_to_Series_Vector(p(k).all,idx);
        begin
          res(k) := new DoblDobl_Dense_Series_Vectors.Vector'(v);
        end;
      end if;
    end loop;
    return res;
  end System_Array_to_Series_VecVec;

  function System_Array_to_Series_VecVec
             ( p : QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
               idx : integer32 := 1 )
             return QuadDobl_Dense_Series_VecVecs.VecVec is

    res : QuadDobl_Dense_Series_VecVecs.VecVec(p'range);

    use QuadDobl_Complex_Poly_Systems;

  begin
    for k in p'range loop
      if p(k) /= null then
        declare
          v : constant QuadDobl_Dense_Series_Vectors.Vector
            := System_to_Series_Vector(p(k).all,idx);
        begin
          res(k) := new QuadDobl_Dense_Series_Vectors.Vector'(v);
        end;
      end if;
    end loop;
    return res;
  end System_Array_to_Series_VecVec;

  function Set_Degree ( i : integer32;
                        d : Standard_Natural_Vectors.Link_to_Vector )
                     return integer32 is

  -- DESCRIPTION :
  --   Returns the degree of the coefficient series,
  --   which is either zero if i = 0, or otherwise d(i).

  begin
    if i = 0 or i > d'last
     then return 0;
     else return integer32(d(i));
    end if;
  end Set_Degree;

  function Set_Dimension ( i : integer32;
                           d : Standard_Natural_Vectors.Link_to_Vector )
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
    --   as a term of a series polynomial.  If idx = 0, then the degree of
    --   the series coefficient is zero and the coefficient is copied.
    --   Otherwise, the variable with index idx in t is taken as the
    --   variable in the truncated power series.  The other variables
    --   in t are relocated to fit a polynomial with one variable less.

      rtm : Standard_Series_Polynomials.Term;
      ord : constant integer32
          := Set_Degree(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      dim : constant integer32
          := Set_Dimension(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
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
        put("degree : "); put(ord,1);
        put(" and degrees : "); put(rtm.dg.all); new_line;
      end if;
      Standard_Series_Polynomials.Add(res,rtm);
      Standard_Series_Polynomials.Clear(rtm.dg);
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series_Polynomial;

  function Polynomial_to_Series_Polynomial
             ( p : DoblDobl_Complex_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return DoblDobl_Series_Polynomials.Poly is

    use DoblDobl_Dense_Series;

    res : DoblDobl_Series_Polynomials.Poly
        := DoblDobl_Series_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in DoblDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Adds the information in the term of a multivariate polynomial
    --   as a term of a series polynomial.  If idx = 0, then the degree of
    --   the series coefficient is zero and the coefficient is copied.
    --   Otherwise, the variable with index idx in t is taken as the
    --   variable in the truncated power series.  The other variables
    --   in t are relocated to fit a polynomial with one variable less.

      rtm : DoblDobl_Series_Polynomials.Term;
      ord : constant integer32
          := Set_Degree(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      dim : constant integer32
          := Set_Dimension(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
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
        put("degree : "); put(ord,1);
        put(" and degrees : "); put(rtm.dg.all); new_line;
      end if;
      DoblDobl_Series_Polynomials.Add(res,rtm);
      DoblDobl_Series_Polynomials.Clear(rtm.dg);
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new DoblDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series_Polynomial;

  function Polynomial_to_Series_Polynomial
             ( p : QuadDobl_Complex_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return QuadDobl_Series_Polynomials.Poly is

    use QuadDobl_Dense_Series;

    res : QuadDobl_Series_Polynomials.Poly
        := QuadDobl_Series_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in QuadDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Adds the information in the term of a multivariate polynomial
    --   as a term of a series polynomial.  If idx = 0, then the degree of
    --   the series coefficient is zero and the coefficient is copied.
    --   Otherwise, the variable with index idx in t is taken as the
    --   variable in the truncated power series.  The other variables
    --   in t are relocated to fit a polynomial with one variable less.

      rtm : QuadDobl_Series_Polynomials.Term;
      ord : constant integer32
          := Set_Degree(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      dim : constant integer32
          := Set_Dimension(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
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
        put("degree : "); put(ord,1);
        put(" and degrees : "); put(rtm.dg.all); new_line;
      end if;
      QuadDobl_Series_Polynomials.Add(res,rtm);
      QuadDobl_Series_Polynomials.Clear(rtm.dg);
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new QuadDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series_Polynomial;

  function Series_Polynomial_to_Polynomial
             ( s : Standard_Series_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Numbers;
    use Standard_Dense_Series;

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in Standard_Series_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Every term in the coefficient of the series in t.cf
    --   with a nonzero coefficient will contribute one monomial
    --   to the result.

      cffs : constant Series := t.cf;
      zero : constant Complex_Number := Create(0.0);
      rtpc : Complex_Number;
      dim1 : integer32;

    begin
      if verbose then
        put("term with degrees :"); put(t.dg.all);
        put(" has series of degree "); put(cffs.deg,1); new_line;
        put_line("the series : "); put(cffs);
      end if;
      if idx = 0 then
        declare
          rt : Standard_Complex_Polynomials.Term;
        begin
          rt.cf := cffs.cff(0);
          rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
          Standard_Complex_Polynomials.Add(res,rt);
          Standard_Complex_Polynomials.Clear(rt.dg);
        end;
      else -- idx > 0
        dim1 := t.dg'last+1;
        for k in 0..cffs.deg loop
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
              if verbose then
                put("the new term has degrees "); put(rt.dg.all); new_line;
                put("and coefficient :"); put(rt.cf); new_line;
              end if;
              Standard_Complex_Polynomials.Add(res,rt);
              Standard_Complex_Polynomials.Clear(rt.dg);
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

  function Series_Polynomial_to_Polynomial
             ( s : DoblDobl_Series_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Dense_Series;

    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in DoblDobl_Series_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Every term in the coefficient of the series in t.cf
    --   with a nonzero coefficient will contribute one monomial
    --   to the result.

      cffs : constant Series := t.cf;
      zero : constant Complex_Number := Create(integer(0));
      rtpc : Complex_Number;
      dim1 : integer32;

    begin
      if verbose then
        put("term with degrees :"); put(t.dg.all);
        put(" has series of degree "); put(cffs.deg,1); new_line;
        put_line("the series : "); put(cffs);
      end if;
      if idx = 0 then
        declare
          rt : DoblDobl_Complex_Polynomials.Term;
        begin
          rt.cf := cffs.cff(0);
          rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
          DoblDobl_Complex_Polynomials.Add(res,rt);
          DoblDobl_Complex_Polynomials.Clear(rt.dg);
        end;
      else -- idx > 0
        dim1 := t.dg'last+1;
        for k in 0..cffs.deg loop
          rtpc := cffs.cff(k);
          if not Equal(rtpc,zero) then
            declare
              rt : DoblDobl_Complex_Polynomials.Term;
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
              if verbose then
                put("the new term has degrees "); put(rt.dg.all); new_line;
                put("and coefficient :"); put(rt.cf); new_line;
              end if;
              DoblDobl_Complex_Polynomials.Add(res,rt);
              DoblDobl_Complex_Polynomials.Clear(rt.dg);
            end;
          end if;
        end loop;
      end if;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new DoblDobl_Series_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
    return res;
  end Series_Polynomial_to_Polynomial;

  function Series_Polynomial_to_Polynomial
             ( s : QuadDobl_Series_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Dense_Series;

    res : QuadDobl_Complex_Polynomials.Poly
        := QuadDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in QuadDobl_Series_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Every term in the coefficient of the series in t.cf
    --   with a nonzero coefficient will contribute one monomial
    --   to the result.

      cffs : constant Series := t.cf;
      zero : constant Complex_Number := Create(integer(0));
      rtpc : Complex_Number;
      dim1 : integer32;

    begin
      if verbose then
        put("term with degrees :"); put(t.dg.all);
        put(" has series of degree "); put(cffs.deg,1); new_line;
        put_line("the series : "); put(cffs);
      end if;
      if idx = 0 then
        declare
          rt : QuadDobl_Complex_Polynomials.Term;
        begin
          rt.cf := cffs.cff(0);
          rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
          QuadDobl_Complex_Polynomials.Add(res,rt);
          QuadDobl_Complex_Polynomials.Clear(rt.dg);
        end;
      else -- idx > 0
        dim1 := t.dg'last+1;
        for k in 0..cffs.deg loop
          rtpc := cffs.cff(k);
          if not Equal(rtpc,zero) then
            declare
              rt : QuadDobl_Complex_Polynomials.Term;
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
              if verbose then
                put("the new term has degrees "); put(rt.dg.all); new_line;
                put("and coefficient :"); put(rt.cf); new_line;
              end if;
              QuadDobl_Complex_Polynomials.Add(res,rt);
              QuadDobl_Complex_Polynomials.Clear(rt.dg);
            end;
          end if;
        end loop;
      end if;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new QuadDobl_Series_Polynomials.Visiting_Iterator(Visit_Term);

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
      if verbose
       then put("converting polynomial "); put(i,1); put_line(" ...");
      end if;
      res(i) := Polynomial_to_Series_Polynomial(p(i),idx,verbose);
    end loop;
    return res;
  end System_to_Series_System;

  function System_to_Series_System
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return DoblDobl_Series_Poly_Systems.Poly_Sys is

    res : DoblDobl_Series_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      if verbose
       then put("converting polynomial "); put(i,1); put_line(" ...");
      end if;
      res(i) := Polynomial_to_Series_Polynomial(p(i),idx,verbose);
    end loop;
    return res;
  end System_to_Series_System;

  function System_to_Series_System
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return QuadDobl_Series_Poly_Systems.Poly_Sys is

    res : QuadDobl_Series_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      if verbose
       then put("converting polynomial "); put(i,1); put_line(" ...");
      end if;
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
      if verbose
       then put("converting series polynomial "); put(i,1); put_line("...");
      end if;
      res(i) := Series_Polynomial_to_Polynomial(s(i),idx,verbose);
    end loop;
    return res;
  end Series_System_to_System;

  function Series_System_to_System
             ( s : DoblDobl_Series_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(s'range);

  begin
    for i in s'range loop
      if verbose
       then put("converting series polynomial "); put(i,1); put_line("...");
      end if;
      res(i) := Series_Polynomial_to_Polynomial(s(i),idx,verbose);
    end loop;
    return res;
  end Series_System_to_System;

  function Series_System_to_System
             ( s : QuadDobl_Series_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(s'range);

  begin
    for i in s'range loop
      if verbose
       then put("converting series polynomial "); put(i,1); put_line("...");
      end if;
      res(i) := Series_Polynomial_to_Polynomial(s(i),idx,verbose);
    end loop;
    return res;
  end Series_System_to_System;

  procedure Set_Degree ( v : in out Standard_Dense_Series_Vectors.Vector;
                         degree : in integer32 ) is
  begin
    for i in v'range loop
      if v(i).deg < degree then
        for k in v(i).deg+1..degree loop
          v(i).cff(k) := Standard_Complex_Numbers.Create(0.0);
        end loop;
      end if;
      v(i).deg := degree;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( v : in out DoblDobl_Dense_Series_Vectors.Vector;
                         degree : in integer32 ) is
  begin
    for i in v'range loop
      if v(i).deg < degree then
        for k in v(i).deg+1..degree loop
          v(i).cff(k) := DoblDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end if;
      v(i).deg := degree;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( v : in out QuadDobl_Dense_Series_Vectors.Vector;
                         degree : in integer32 ) is
  begin
    for i in v'range loop
      if v(i).deg < degree then
        for k in v(i).deg+1..degree loop
          v(i).cff(k) := QuadDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end if;
      v(i).deg := degree;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( m : in out Standard_Dense_Series_Matrices.Matrix;
                         degree : in integer32 ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        if m(i,j).deg < degree then
          for k in m(i,j).deg+1..degree loop
            m(i,j).cff(k) := Standard_Complex_Numbers.Create(0.0);
          end loop;
        end if;
        m(i,j).deg := degree;
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( m : in out DoblDobl_Dense_Series_Matrices.Matrix;
                         degree : in integer32 ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        if m(i,j).deg < degree then
          for k in m(i,j).deg+1..degree loop
            m(i,j).cff(k) := DoblDobl_Complex_Numbers.Create(integer(0));
          end loop;
        end if;
        m(i,j).deg := degree;
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( m : in out QuadDobl_Dense_Series_Matrices.Matrix;
                         degree : in integer32 ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        if m(i,j).deg < degree then
          for k in m(i,j).deg+1..degree loop
            m(i,j).cff(k) := QuadDobl_Complex_Numbers.Create(integer(0));
          end loop;
        end if;
        m(i,j).deg := degree;
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( p : in out Standard_Series_Polynomials.Poly;
                         degree : in integer32 ) is

    procedure Change_Degree ( t : in out Standard_Series_Polynomials.Term;
                              c : out boolean ) is
    begin
      if t.cf.deg < degree then
        for k in t.cf.deg+1..degree loop
          t.cf.cff(k) := Standard_Complex_Numbers.Create(0.0);
        end loop;
      end if;
      t.cf.deg := degree;
      c := true;
    end Change_Degree;
    procedure Change_Degrees is
      new Standard_Series_Polynomials.Changing_Iterator(Change_Degree);

  begin
    Change_Degrees(p);
  end Set_Degree;

  procedure Set_Degree ( p : in out DoblDobl_Series_Polynomials.Poly;
                         degree : in integer32 ) is

    procedure Change_Degree ( t : in out DoblDobl_Series_Polynomials.Term;
                              c : out boolean ) is
    begin
      if t.cf.deg < degree then
        for k in t.cf.deg+1..degree loop
          t.cf.cff(k) := DoblDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end if;
      t.cf.deg := degree;
      c := true;
    end Change_Degree;
    procedure Change_Degrees is
      new DoblDobl_Series_Polynomials.Changing_Iterator(Change_Degree);

  begin
    Change_Degrees(p);
  end Set_Degree;

  procedure Set_Degree ( p : in out QuadDobl_Series_Polynomials.Poly;
                        degree : in integer32 ) is

    procedure Change_Degree ( t : in out QuadDobl_Series_Polynomials.Term;
                              c : out boolean ) is
    begin
      if t.cf.deg < degree then
        for k in t.cf.deg+1..degree loop
          t.cf.cff(k) := QuadDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end if;
      t.cf.deg := degree;
      c := true;
    end Change_Degree;
    procedure Change_Degrees is
      new QuadDobl_Series_Polynomials.Changing_Iterator(Change_Degree);

  begin
    Change_Degrees(p);
  end Set_Degree;

  procedure Set_Degree ( p : in out Standard_Series_Poly_Systems.Poly_Sys;
                         degree : in integer32 ) is
  begin
    for i in p'range loop
      Set_Degree(p(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( p : in out DoblDobl_Series_Poly_Systems.Poly_Sys;
                         degree : in integer32 ) is
  begin
    for i in p'range loop
      Set_Degree(p(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( p : in out QuadDobl_Series_Poly_Systems.Poly_Sys;
                         degree : in integer32 ) is
  begin
    for i in p'range loop
      Set_Degree(p(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( jm : in out Standard_Series_Jaco_Matrices.Jaco_Mat;
                         degree : in integer32 ) is
  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        Set_Degree(jm(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( jm : in out DoblDobl_Series_Jaco_Matrices.Jaco_Mat;
                         degree : in integer32 ) is
  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        Set_Degree(jm(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( jm : in out QuadDobl_Series_Jaco_Matrices.Jaco_Mat;
                         degree : in integer32 ) is
  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        Set_Degree(jm(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Filter ( s : in out Standard_Dense_Series_Vectors.Vector;
                     tol : in double_float ) is
  begin
    for i in s'range loop
      Standard_Dense_Series.Filter(s(i),tol);
    end loop;
  end Filter;

  procedure Filter ( s : in out DoblDobl_Dense_Series_Vectors.Vector;
                     tol : in double_float ) is
  begin
    for i in s'range loop
      DoblDobl_Dense_Series.Filter(s(i),tol);
    end loop;
  end Filter;

  procedure Filter ( s : in out QuadDobl_Dense_Series_Vectors.Vector;
                     tol : in double_float ) is
  begin
    for i in s'range loop
      QuadDobl_Dense_Series.Filter(s(i),tol);
    end loop;
  end Filter;
 
end Series_and_Polynomials;
