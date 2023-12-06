with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with TripDobl_Complex_Numbers;
with TripDobl_Complex_Numbers_io;       use TripDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with PentDobl_Complex_Numbers;
with PentDobl_Complex_Numbers_io;       use PentDobl_Complex_Numbers_io;
with OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers_io;       use OctoDobl_Complex_Numbers_io;
with DecaDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers_io;       use DecaDobl_Complex_Numbers_io;
with HexaDobl_Complex_Numbers;
with HexaDobl_Complex_Numbers_io;       use HexaDobl_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Complex_Series_Functions;
with DoblDobl_Complex_Series_Functions;
with TripDobl_Complex_Series_Functions;
with QuadDobl_Complex_Series_Functions;
with PentDobl_Complex_Series_Functions;
with OctoDobl_Complex_Series_Functions;
with DecaDobl_Complex_Series_Functions;
with HexaDobl_Complex_Series_Functions;
with Standard_Complex_Series_io;        use Standard_Complex_Series_io;
with DoblDobl_Complex_Series_io;        use DoblDobl_Complex_Series_io;
with TripDobl_Complex_Series_io;        use TripDobl_Complex_Series_io;
with QuadDobl_Complex_Series_io;        use QuadDobl_Complex_Series_io;
with PentDobl_Complex_Series_io;        use PentDobl_Complex_Series_io;
with OctoDobl_Complex_Series_io;        use OctoDobl_Complex_Series_io;
with DecaDobl_Complex_Series_io;        use DecaDobl_Complex_Series_io;
with HexaDobl_Complex_Series_io;        use HexaDobl_Complex_Series_io;

package body Complex_Series_and_Polynomials is

  function Series_to_Polynomial
             ( s : Standard_Complex_Series.Series )
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
             ( s : DoblDobl_Complex_Series.Series )
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
             ( s : TripDobl_Complex_Series.Series )
             return TripDobl_Complex_Polynomials.Poly is

    use TripDobl_Complex_Numbers;

    res : TripDobl_Complex_Polynomials.Poly
        := TripDobl_Complex_Polynomials.Null_Poly;
    zero : constant Complex_Number := Create(integer(0));

    procedure Add_Term ( k : in integer32; c : in Complex_Number ) is

    -- DESCRIPTION :
    --   Adds c*t^k to the polynomial res.

      t : TripDobl_Complex_Polynomials.Term;

    begin
      t.cf := c;
      t.dg := new Standard_Natural_Vectors.Vector(1..1);
      t.dg(1) := natural32(k);
      TripDobl_Complex_Polynomials.Add(res,t);
      TripDobl_Complex_Polynomials.Clear(t.dg);
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
             ( s : QuadDobl_Complex_Series.Series )
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

  function Series_to_Polynomial
             ( s : PentDobl_Complex_Series.Series )
             return PentDobl_Complex_Polynomials.Poly is

    use PentDobl_Complex_Numbers;

    res : PentDobl_Complex_Polynomials.Poly
        := PentDobl_Complex_Polynomials.Null_Poly;
    zero : constant Complex_Number := Create(integer(0));

    procedure Add_Term ( k : in integer32; c : in Complex_Number ) is

    -- DESCRIPTION :
    --   Adds c*t^k to the polynomial res.

      t : PentDobl_Complex_Polynomials.Term;

    begin
      t.cf := c;
      t.dg := new Standard_Natural_Vectors.Vector(1..1);
      t.dg(1) := natural32(k);
      PentDobl_Complex_Polynomials.Add(res,t);
      PentDobl_Complex_Polynomials.Clear(t.dg);
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
             ( s : OctoDobl_Complex_Series.Series )
             return OctoDobl_Complex_Polynomials.Poly is

    use OctoDobl_Complex_Numbers;

    res : OctoDobl_Complex_Polynomials.Poly
        := OctoDobl_Complex_Polynomials.Null_Poly;
    zero : constant Complex_Number := Create(integer(0));

    procedure Add_Term ( k : in integer32; c : in Complex_Number ) is

    -- DESCRIPTION :
    --   Adds c*t^k to the polynomial res.

      t : OctoDobl_Complex_Polynomials.Term;

    begin
      t.cf := c;
      t.dg := new Standard_Natural_Vectors.Vector(1..1);
      t.dg(1) := natural32(k);
      OctoDobl_Complex_Polynomials.Add(res,t);
      OctoDobl_Complex_Polynomials.Clear(t.dg);
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
             ( s : DecaDobl_Complex_Series.Series )
             return DecaDobl_Complex_Polynomials.Poly is

    use DecaDobl_Complex_Numbers;

    res : DecaDobl_Complex_Polynomials.Poly
        := DecaDobl_Complex_Polynomials.Null_Poly;
    zero : constant Complex_Number := Create(integer(0));

    procedure Add_Term ( k : in integer32; c : in Complex_Number ) is

    -- DESCRIPTION :
    --   Adds c*t^k to the polynomial res.

      t : DecaDobl_Complex_Polynomials.Term;

    begin
      t.cf := c;
      t.dg := new Standard_Natural_Vectors.Vector(1..1);
      t.dg(1) := natural32(k);
      DecaDobl_Complex_Polynomials.Add(res,t);
      DecaDobl_Complex_Polynomials.Clear(t.dg);
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
             ( s : HexaDobl_Complex_Series.Series )
             return HexaDobl_Complex_Polynomials.Poly is

    use HexaDobl_Complex_Numbers;

    res : HexaDobl_Complex_Polynomials.Poly
        := HexaDobl_Complex_Polynomials.Null_Poly;
    zero : constant Complex_Number := Create(integer(0));

    procedure Add_Term ( k : in integer32; c : in Complex_Number ) is

    -- DESCRIPTION :
    --   Adds c*t^k to the polynomial res.

      t : HexaDobl_Complex_Polynomials.Term;

    begin
      t.cf := c;
      t.dg := new Standard_Natural_Vectors.Vector(1..1);
      t.dg(1) := natural32(k);
      HexaDobl_Complex_Polynomials.Add(res,t);
      HexaDobl_Complex_Polynomials.Clear(t.dg);
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
             return Standard_Complex_Series.Series is

    res : Standard_Complex_Series.Link_to_Series;

    procedure Visit_Term ( t : in Standard_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Assigns t.cf to the res.cff(t.dg(idx))
    --   and updates the degree if needed.

      d : constant integer32 := integer32(t.dg(idx));

    begin
      if d > res.deg
       then Standard_Complex_Series.Set_Degree(res,d);
      end if;
      res.cff(d) := t.cf;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    res := Standard_Complex_Series.Create(0);
    Visit_Terms(p);
    return res.all;
  end Polynomial_to_Series;

  function Polynomial_to_Series
             ( p : DoblDobl_Complex_Polynomials.Poly;
               idx : integer32 := 1 )
             return DoblDobl_Complex_Series.Series is

    res : DoblDobl_Complex_Series.Link_to_Series;

    procedure Visit_Term ( t : in DoblDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Assigns t.cf to the res.cff(t.dg(idx))
    --   and updates the degree if needed.

      d : constant integer32 := integer32(t.dg(idx));

    begin
      if d > res.deg
       then DoblDobl_Complex_Series.Set_Degree(res,d);
      end if;
      res.cff(d) := t.cf;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new DoblDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    res := DoblDobl_Complex_Series.Create(0);
    Visit_Terms(p);
    return res.all;
  end Polynomial_to_Series;

  function Polynomial_to_Series
             ( p : TripDobl_Complex_Polynomials.Poly;
               idx : integer32 := 1 )
             return TripDobl_Complex_Series.Series is

    res : TripDobl_Complex_Series.Link_to_Series;

    procedure Visit_Term ( t : in TripDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Assigns t.cf to the res.cff(t.dg(idx))
    --   and updates the degree if needed.

      d : constant integer32 := integer32(t.dg(idx));

    begin
      if d > res.deg
       then TripDobl_Complex_Series.Set_Degree(res,d);
      end if;
      res.cff(d) := t.cf;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new TripDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    res := TripDobl_Complex_Series.Create(0);
    Visit_Terms(p);
    return res.all;
  end Polynomial_to_Series;

  function Polynomial_to_Series
             ( p : QuadDobl_Complex_Polynomials.Poly;
               idx : integer32 := 1 )
             return QuadDobl_Complex_Series.Series is

    res : QuadDobl_Complex_Series.Link_to_Series;

    procedure Visit_Term ( t : in QuadDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Assigns t.cf to the res.cff(t.dg(idx))
    --   and updates the degree if needed.

      d : constant integer32 := integer32(t.dg(idx));

    begin
      if d > res.deg
       then QuadDobl_Complex_Series.Set_Degree(res,d);
      end if;
      res.cff(d) := t.cf;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new QuadDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    res := QuadDobl_Complex_Series.Create(0);
    Visit_Terms(p);
    return res.all;
  end Polynomial_to_Series;

  function Polynomial_to_Series
             ( p : PentDobl_Complex_Polynomials.Poly;
               idx : integer32 := 1 )
             return PentDobl_Complex_Series.Series is

    res : PentDobl_Complex_Series.Link_to_Series;

    procedure Visit_Term ( t : in PentDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Assigns t.cf to the res.cff(t.dg(idx))
    --   and updates the degree if needed.

      d : constant integer32 := integer32(t.dg(idx));

    begin
      if d > res.deg
       then PentDobl_Complex_Series.Set_Degree(res,d);
      end if;
      res.cff(d) := t.cf;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new PentDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    res := PentDobl_Complex_Series.Create(0);
    Visit_Terms(p);
    return res.all;
  end Polynomial_to_Series;

  function Polynomial_to_Series
             ( p : OctoDobl_Complex_Polynomials.Poly;
               idx : integer32 := 1 )
             return OctoDobl_Complex_Series.Series is

    res : OctoDobl_Complex_Series.Link_to_Series;

    procedure Visit_Term ( t : in OctoDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Assigns t.cf to the res.cff(t.dg(idx))
    --   and updates the degree if needed.

      d : constant integer32 := integer32(t.dg(idx));

    begin
      if d > res.deg
       then OctoDobl_Complex_Series.Set_Degree(res,d);
      end if;
      res.cff(d) := t.cf;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new OctoDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    res := OctoDobl_Complex_Series.Create(0);
    Visit_Terms(p);
    return res.all;
  end Polynomial_to_Series;

  function Polynomial_to_Series
             ( p : DecaDobl_Complex_Polynomials.Poly;
               idx : integer32 := 1 )
             return DecaDobl_Complex_Series.Series is

    res : DecaDobl_Complex_Series.Link_to_Series;

    procedure Visit_Term ( t : in DecaDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Assigns t.cf to the res.cff(t.dg(idx))
    --   and updates the degree if needed.

      d : constant integer32 := integer32(t.dg(idx));

    begin
      if d > res.deg
       then DecaDobl_Complex_Series.Set_Degree(res,d);
      end if;
      res.cff(d) := t.cf;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new DecaDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    res := DecaDobl_Complex_Series.Create(0);
    Visit_Terms(p);
    return res.all;
  end Polynomial_to_Series;

  function Polynomial_to_Series
             ( p : HexaDobl_Complex_Polynomials.Poly;
               idx : integer32 := 1 )
             return HexaDobl_Complex_Series.Series is

    res : HexaDobl_Complex_Series.Link_to_Series;

    procedure Visit_Term ( t : in HexaDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Assigns t.cf to the res.cff(t.dg(idx))
    --   and updates the degree if needed.

      d : constant integer32 := integer32(t.dg(idx));

    begin
      if d > res.deg
       then HexaDobl_Complex_Series.Set_Degree(res,d);
      end if;
      res.cff(d) := t.cf;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new HexaDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    res := HexaDobl_Complex_Series.Create(0);
    Visit_Terms(p);
    return res.all;
  end Polynomial_to_Series;

  function Series_Vector_to_System
             ( v : Standard_Complex_Series_Vectors.Vector )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(v'range);

    use Standard_Complex_Series;

  begin
    for k in v'range loop
      if v(k) /= null
       then res(k) := Series_to_Polynomial(v(k).all);
      end if;
    end loop;
    return res;
  end Series_Vector_to_System;

  function Series_Vector_to_System
             ( v : DoblDobl_Complex_Series_Vectors.Vector )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(v'range);

    use DoblDobl_Complex_Series;

  begin
    for k in v'range loop
      if v(k) /= null
       then res(k) := Series_to_Polynomial(v(k).all);
      end if;
    end loop;
    return res;
  end Series_Vector_to_System;

  function Series_Vector_to_System
             ( v : TripDobl_Complex_Series_Vectors.Vector )
             return TripDobl_Complex_Poly_Systems.Poly_Sys is

    res : TripDobl_Complex_Poly_Systems.Poly_Sys(v'range);

    use TripDobl_Complex_Series;

  begin
    for k in v'range loop
      if v(k) /= null
       then res(k) := Series_to_Polynomial(v(k).all);
      end if;
    end loop;
    return res;
  end Series_Vector_to_System;

  function Series_Vector_to_System
             ( v : QuadDobl_Complex_Series_Vectors.Vector )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(v'range);

    use QuadDobl_Complex_Series;

  begin
    for k in v'range loop
      if v(k) /= null
       then res(k) := Series_to_Polynomial(v(k).all);
      end if;
    end loop;
    return res;
  end Series_Vector_to_System;

  function Series_Vector_to_System
             ( v : PentDobl_Complex_Series_Vectors.Vector )
             return PentDobl_Complex_Poly_Systems.Poly_Sys is

    res : PentDobl_Complex_Poly_Systems.Poly_Sys(v'range);

    use PentDobl_Complex_Series;

  begin
    for k in v'range loop
      if v(k) /= null
       then res(k) := Series_to_Polynomial(v(k).all);
      end if;
    end loop;
    return res;
  end Series_Vector_to_System;

  function Series_Vector_to_System
             ( v : OctoDobl_Complex_Series_Vectors.Vector )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys is

    res : OctoDobl_Complex_Poly_Systems.Poly_Sys(v'range);

    use OctoDobl_Complex_Series;

  begin
    for k in v'range loop
      if v(k) /= null
       then res(k) := Series_to_Polynomial(v(k).all);
      end if;
    end loop;
    return res;
  end Series_Vector_to_System;

  function Series_Vector_to_System
             ( v : DecaDobl_Complex_Series_Vectors.Vector )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys is

    res : DecaDobl_Complex_Poly_Systems.Poly_Sys(v'range);

    use DecaDobl_Complex_Series;

  begin
    for k in v'range loop
      if v(k) /= null
       then res(k) := Series_to_Polynomial(v(k).all);
      end if;
    end loop;
    return res;
  end Series_Vector_to_System;

  function Series_Vector_to_System
             ( v : HexaDobl_Complex_Series_Vectors.Vector )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys is

    res : HexaDobl_Complex_Poly_Systems.Poly_Sys(v'range);

    use HexaDobl_Complex_Series;

  begin
    for k in v'range loop
      if v(k) /= null
       then res(k) := Series_to_Polynomial(v(k).all);
      end if;
    end loop;
    return res;
  end Series_Vector_to_System;

  function System_to_Series_Vector
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 1 )
             return Standard_Complex_Series_Vectors.Vector is

    res : Standard_Complex_Series_Vectors.Vector(p'range);

  begin
    for k in p'range loop
      declare
        s : constant Standard_Complex_Series.Series
          := Polynomial_to_Series(p(k),idx);
      begin
        res(k) := new Standard_Complex_Series.Series'(s);
      end;
    end loop;
    return res;
  end System_to_Series_Vector;

  function System_to_Series_Vector
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 1 )
             return DoblDobl_Complex_Series_Vectors.Vector is

    res : DoblDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    for k in p'range loop
      declare
        s : constant DoblDobl_Complex_Series.Series
          := Polynomial_to_Series(p(k),idx);
      begin
        res(k) := new DoblDobl_Complex_Series.Series'(s);
      end;
    end loop;
    return res;
  end System_to_Series_Vector;

  function System_to_Series_Vector
             ( p : TripDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 1 )
             return TripDobl_Complex_Series_Vectors.Vector is

    res : TripDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    for k in p'range loop
      declare
        s : constant TripDobl_Complex_Series.Series
          := Polynomial_to_Series(p(k),idx);
      begin
        res(k) := new TripDobl_Complex_Series.Series'(s);
      end;
    end loop;
    return res;
  end System_to_Series_Vector;

  function System_to_Series_Vector
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 1 )
             return QuadDobl_Complex_Series_Vectors.Vector is

    res : QuadDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    for k in p'range loop
      declare
        s : constant QuadDobl_Complex_Series.Series
          := Polynomial_to_Series(p(k),idx);
      begin
        res(k) := new QuadDobl_Complex_Series.Series'(s);
      end;
    end loop;
    return res;
  end System_to_Series_Vector;

  function System_to_Series_Vector
             ( p : PentDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 1 )
             return PentDobl_Complex_Series_Vectors.Vector is

    res : PentDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    for k in p'range loop
      declare
        s : constant PentDobl_Complex_Series.Series
          := Polynomial_to_Series(p(k),idx);
      begin
        res(k) := new PentDobl_Complex_Series.Series'(s);
      end;
    end loop;
    return res;
  end System_to_Series_Vector;

  function System_to_Series_Vector
             ( p : OctoDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 1 )
             return OctoDobl_Complex_Series_Vectors.Vector is

    res : OctoDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    for k in p'range loop
      declare
        s : constant OctoDobl_Complex_Series.Series
          := Polynomial_to_Series(p(k),idx);
      begin
        res(k) := new OctoDobl_Complex_Series.Series'(s);
      end;
    end loop;
    return res;
  end System_to_Series_Vector;

  function System_to_Series_Vector
             ( p : DecaDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 1 )
             return DecaDobl_Complex_Series_Vectors.Vector is

    res : DecaDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    for k in p'range loop
      declare
        s : constant DecaDobl_Complex_Series.Series
          := Polynomial_to_Series(p(k),idx);
      begin
        res(k) := new DecaDobl_Complex_Series.Series'(s);
      end;
    end loop;
    return res;
  end System_to_Series_Vector;

  function System_to_Series_Vector
             ( p : HexaDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 1 )
             return HexaDobl_Complex_Series_Vectors.Vector is

    res : HexaDobl_Complex_Series_Vectors.Vector(p'range);

  begin
    for k in p'range loop
      declare
        s : constant HexaDobl_Complex_Series.Series
          := Polynomial_to_Series(p(k),idx);
      begin
        res(k) := new HexaDobl_Complex_Series.Series'(s);
      end;
    end loop;
    return res;
  end System_to_Series_Vector;

  function Series_VecVec_to_System_Array
             ( v : Standard_Complex_Series_VecVecs.VecVec )
             return Standard_Complex_Poly_Systems.Array_of_Poly_Sys is

    res : Standard_Complex_Poly_Systems.Array_of_Poly_Sys(v'range);

    use Standard_Complex_Series_Vectors;

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
             ( v : DoblDobl_Complex_Series_VecVecs.VecVec )
             return DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys(v'range);

    use DoblDobl_Complex_Series_Vectors;

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
             ( v : TripDobl_Complex_Series_VecVecs.VecVec )
             return TripDobl_Complex_Poly_Systems.Array_of_Poly_Sys is

    res : TripDobl_Complex_Poly_Systems.Array_of_Poly_Sys(v'range);

    use TripDobl_Complex_Series_Vectors;

  begin
    for k in v'range loop
      if v(k) /= null then
        declare
          p : constant TripDobl_Complex_Poly_Systems.Poly_Sys
            := Series_Vector_to_System(v(k).all);
        begin
          res(k) := new TripDobl_Complex_Poly_Systems.Poly_Sys'(p);
        end;
      end if;
    end loop;
    return res;
  end Series_VecVec_to_System_Array;

  function Series_VecVec_to_System_Array
             ( v : QuadDobl_Complex_Series_VecVecs.VecVec )
             return QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys(v'range);

    use QuadDobl_Complex_Series_Vectors;

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

  function Series_VecVec_to_System_Array
             ( v : PentDobl_Complex_Series_VecVecs.VecVec )
             return PentDobl_Complex_Poly_Systems.Array_of_Poly_Sys is

    res : PentDobl_Complex_Poly_Systems.Array_of_Poly_Sys(v'range);

    use PentDobl_Complex_Series_Vectors;

  begin
    for k in v'range loop
      if v(k) /= null then
        declare
          p : constant PentDobl_Complex_Poly_Systems.Poly_Sys
            := Series_Vector_to_System(v(k).all);
        begin
          res(k) := new PentDobl_Complex_Poly_Systems.Poly_Sys'(p);
        end;
      end if;
    end loop;
    return res;
  end Series_VecVec_to_System_Array;

  function Series_VecVec_to_System_Array
             ( v : OctoDobl_Complex_Series_VecVecs.VecVec )
             return OctoDobl_Complex_Poly_Systems.Array_of_Poly_Sys is

    res : OctoDobl_Complex_Poly_Systems.Array_of_Poly_Sys(v'range);

    use OctoDobl_Complex_Series_Vectors;

  begin
    for k in v'range loop
      if v(k) /= null then
        declare
          p : constant OctoDobl_Complex_Poly_Systems.Poly_Sys
            := Series_Vector_to_System(v(k).all);
        begin
          res(k) := new OctoDobl_Complex_Poly_Systems.Poly_Sys'(p);
        end;
      end if;
    end loop;
    return res;
  end Series_VecVec_to_System_Array;

  function Series_VecVec_to_System_Array
             ( v : DecaDobl_Complex_Series_VecVecs.VecVec )
             return DecaDobl_Complex_Poly_Systems.Array_of_Poly_Sys is

    res : DecaDobl_Complex_Poly_Systems.Array_of_Poly_Sys(v'range);

    use DecaDobl_Complex_Series_Vectors;

  begin
    for k in v'range loop
      if v(k) /= null then
        declare
          p : constant DecaDobl_Complex_Poly_Systems.Poly_Sys
            := Series_Vector_to_System(v(k).all);
        begin
          res(k) := new DecaDobl_Complex_Poly_Systems.Poly_Sys'(p);
        end;
      end if;
    end loop;
    return res;
  end Series_VecVec_to_System_Array;

  function Series_VecVec_to_System_Array
             ( v : HexaDobl_Complex_Series_VecVecs.VecVec )
             return HexaDobl_Complex_Poly_Systems.Array_of_Poly_Sys is

    res : HexaDobl_Complex_Poly_Systems.Array_of_Poly_Sys(v'range);

    use HexaDobl_Complex_Series_Vectors;

  begin
    for k in v'range loop
      if v(k) /= null then
        declare
          p : constant HexaDobl_Complex_Poly_Systems.Poly_Sys
            := Series_Vector_to_System(v(k).all);
        begin
          res(k) := new HexaDobl_Complex_Poly_Systems.Poly_Sys'(p);
        end;
      end if;
    end loop;
    return res;
  end Series_VecVec_to_System_Array;

  function System_Array_to_Series_VecVec
             ( p : Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
               idx : integer32 := 1 )
             return Standard_Complex_Series_VecVecs.VecVec is

    res : Standard_Complex_Series_VecVecs.VecVec(p'range);

    use Standard_Complex_Poly_Systems;

  begin
    for k in p'range loop
      if p(k) /= null then
        declare
          v : constant Standard_Complex_Series_Vectors.Vector
            := System_to_Series_Vector(p(k).all,idx);
        begin
          res(k) := new Standard_Complex_Series_Vectors.Vector'(v);
        end;
      end if;
    end loop;
    return res;
  end System_Array_to_Series_VecVec;

  function System_Array_to_Series_VecVec
             ( p : DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
               idx : integer32 := 1 )
             return DoblDobl_Complex_Series_VecVecs.VecVec is

    res : DoblDobl_Complex_Series_VecVecs.VecVec(p'range);

    use DoblDobl_Complex_Poly_Systems;

  begin
    for k in p'range loop
      if p(k) /= null then
        declare
          v : constant DoblDobl_Complex_Series_Vectors.Vector
            := System_to_Series_Vector(p(k).all,idx);
        begin
          res(k) := new DoblDobl_Complex_Series_Vectors.Vector'(v);
        end;
      end if;
    end loop;
    return res;
  end System_Array_to_Series_VecVec;

  function System_Array_to_Series_VecVec
             ( p : TripDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
               idx : integer32 := 1 )
             return TripDobl_Complex_Series_VecVecs.VecVec is

    res : TripDobl_Complex_Series_VecVecs.VecVec(p'range);

    use TripDobl_Complex_Poly_Systems;

  begin
    for k in p'range loop
      if p(k) /= null then
        declare
          v : constant TripDobl_Complex_Series_Vectors.Vector
            := System_to_Series_Vector(p(k).all,idx);
        begin
          res(k) := new TripDobl_Complex_Series_Vectors.Vector'(v);
        end;
      end if;
    end loop;
    return res;
  end System_Array_to_Series_VecVec;

  function System_Array_to_Series_VecVec
             ( p : QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
               idx : integer32 := 1 )
             return QuadDobl_Complex_Series_VecVecs.VecVec is

    res : QuadDobl_Complex_Series_VecVecs.VecVec(p'range);

    use QuadDobl_Complex_Poly_Systems;

  begin
    for k in p'range loop
      if p(k) /= null then
        declare
          v : constant QuadDobl_Complex_Series_Vectors.Vector
            := System_to_Series_Vector(p(k).all,idx);
        begin
          res(k) := new QuadDobl_Complex_Series_Vectors.Vector'(v);
        end;
      end if;
    end loop;
    return res;
  end System_Array_to_Series_VecVec;

  function System_Array_to_Series_VecVec
             ( p : PentDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
               idx : integer32 := 1 )
             return PentDobl_Complex_Series_VecVecs.VecVec is

    res : PentDobl_Complex_Series_VecVecs.VecVec(p'range);

    use PentDobl_Complex_Poly_Systems;

  begin
    for k in p'range loop
      if p(k) /= null then
        declare
          v : constant PentDobl_Complex_Series_Vectors.Vector
            := System_to_Series_Vector(p(k).all,idx);
        begin
          res(k) := new PentDobl_Complex_Series_Vectors.Vector'(v);
        end;
      end if;
    end loop;
    return res;
  end System_Array_to_Series_VecVec;

  function System_Array_to_Series_VecVec
             ( p : OctoDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
               idx : integer32 := 1 )
             return OctoDobl_Complex_Series_VecVecs.VecVec is

    res : OctoDobl_Complex_Series_VecVecs.VecVec(p'range);

    use OctoDobl_Complex_Poly_Systems;

  begin
    for k in p'range loop
      if p(k) /= null then
        declare
          v : constant OctoDobl_Complex_Series_Vectors.Vector
            := System_to_Series_Vector(p(k).all,idx);
        begin
          res(k) := new OctoDobl_Complex_Series_Vectors.Vector'(v);
        end;
      end if;
    end loop;
    return res;
  end System_Array_to_Series_VecVec;

  function System_Array_to_Series_VecVec
             ( p : DecaDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
               idx : integer32 := 1 )
             return DecaDobl_Complex_Series_VecVecs.VecVec is

    res : DecaDobl_Complex_Series_VecVecs.VecVec(p'range);

    use DecaDobl_Complex_Poly_Systems;

  begin
    for k in p'range loop
      if p(k) /= null then
        declare
          v : constant DecaDobl_Complex_Series_Vectors.Vector
            := System_to_Series_Vector(p(k).all,idx);
        begin
          res(k) := new DecaDobl_Complex_Series_Vectors.Vector'(v);
        end;
      end if;
    end loop;
    return res;
  end System_Array_to_Series_VecVec;

  function System_Array_to_Series_VecVec
             ( p : HexaDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
               idx : integer32 := 1 )
             return HexaDobl_Complex_Series_VecVecs.VecVec is

    res : HexaDobl_Complex_Series_VecVecs.VecVec(p'range);

    use HexaDobl_Complex_Poly_Systems;

  begin
    for k in p'range loop
      if p(k) /= null then
        declare
          v : constant HexaDobl_Complex_Series_Vectors.Vector
            := System_to_Series_Vector(p(k).all,idx);
        begin
          res(k) := new HexaDobl_Complex_Series_Vectors.Vector'(v);
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
             return Standard_CSeries_Polynomials.Poly is

    use Standard_Complex_Series;

    res : Standard_CSeries_Polynomials.Poly
        := Standard_CSeries_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in Standard_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Adds the information in the term of a multivariate polynomial
    --   as a term of a series polynomial.  If idx = 0, then the degree of
    --   the series coefficient is zero and the coefficient is copied.
    --   Otherwise, the variable with index idx in t is taken as the
    --   variable in the truncated power series.  The other variables
    --   in t are relocated to fit a polynomial with one variable less.

      rtm : Standard_CSeries_Polynomials.Term;
      ord : constant integer32
          := Set_Degree(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      dim : constant integer32
          := Set_Dimension(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      rcf : constant Standard_Complex_Series.Link_to_Series
          := Standard_Complex_Series.Create(0,ord);
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
      Standard_CSeries_Polynomials.Add(res,rtm);
      Standard_CSeries_Polynomials.Clear(rtm.dg);
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
             return DoblDobl_CSeries_Polynomials.Poly is

    use DoblDobl_Complex_Series;

    res : DoblDobl_CSeries_Polynomials.Poly
        := DoblDobl_CSeries_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in DoblDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Adds the information in the term of a multivariate polynomial
    --   as a term of a series polynomial.  If idx = 0, then the degree of
    --   the series coefficient is zero and the coefficient is copied.
    --   Otherwise, the variable with index idx in t is taken as the
    --   variable in the truncated power series.  The other variables
    --   in t are relocated to fit a polynomial with one variable less.

      rtm : DoblDobl_CSeries_Polynomials.Term;
      ord : constant integer32
          := Set_Degree(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      dim : constant integer32
          := Set_Dimension(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      rcf : constant DoblDobl_Complex_Series.Link_to_Series
          := DoblDobl_Complex_Series.Create(0,ord);
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
      DoblDobl_CSeries_Polynomials.Add(res,rtm);
      DoblDobl_CSeries_Polynomials.Clear(rtm.dg);
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new DoblDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series_Polynomial;

  function Polynomial_to_Series_Polynomial
             ( p : TripDobl_Complex_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return TripDobl_CSeries_Polynomials.Poly is

    use TripDobl_Complex_Series;

    res : TripDobl_CSeries_Polynomials.Poly
        := TripDobl_CSeries_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in TripDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Adds the information in the term of a multivariate polynomial
    --   as a term of a series polynomial.  If idx = 0, then the degree of
    --   the series coefficient is zero and the coefficient is copied.
    --   Otherwise, the variable with index idx in t is taken as the
    --   variable in the truncated power series.  The other variables
    --   in t are relocated to fit a polynomial with one variable less.

      rtm : TripDobl_CSeries_Polynomials.Term;
      ord : constant integer32
          := Set_Degree(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      dim : constant integer32
          := Set_Dimension(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      rcf : constant TripDobl_Complex_Series.Link_to_Series
          := TripDobl_Complex_Series.Create(0,ord);
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
      TripDobl_CSeries_Polynomials.Add(res,rtm);
      TripDobl_CSeries_Polynomials.Clear(rtm.dg);
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new TripDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series_Polynomial;

  function Polynomial_to_Series_Polynomial
             ( p : QuadDobl_Complex_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return QuadDobl_CSeries_Polynomials.Poly is

    use QuadDobl_Complex_Series;

    res : QuadDobl_CSeries_Polynomials.Poly
        := QuadDobl_CSeries_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in QuadDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Adds the information in the term of a multivariate polynomial
    --   as a term of a series polynomial.  If idx = 0, then the degree of
    --   the series coefficient is zero and the coefficient is copied.
    --   Otherwise, the variable with index idx in t is taken as the
    --   variable in the truncated power series.  The other variables
    --   in t are relocated to fit a polynomial with one variable less.

      rtm : QuadDobl_CSeries_Polynomials.Term;
      ord : constant integer32
          := Set_Degree(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      dim : constant integer32
          := Set_Dimension(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      rcf : constant QuadDobl_Complex_Series.Link_to_Series
          := QuadDobl_Complex_Series.Create(0,ord);
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
      QuadDobl_CSeries_Polynomials.Add(res,rtm);
      QuadDobl_CSeries_Polynomials.Clear(rtm.dg);
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new QuadDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series_Polynomial;

  function Polynomial_to_Series_Polynomial
             ( p : PentDobl_Complex_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return PentDobl_CSeries_Polynomials.Poly is

    use PentDobl_Complex_Series;

    res : PentDobl_CSeries_Polynomials.Poly
        := PentDobl_CSeries_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in PentDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Adds the information in the term of a multivariate polynomial
    --   as a term of a series polynomial.  If idx = 0, then the degree of
    --   the series coefficient is zero and the coefficient is copied.
    --   Otherwise, the variable with index idx in t is taken as the
    --   variable in the truncated power series.  The other variables
    --   in t are relocated to fit a polynomial with one variable less.

      rtm : PentDobl_CSeries_Polynomials.Term;
      ord : constant integer32
          := Set_Degree(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      dim : constant integer32
          := Set_Dimension(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      rcf : constant PentDobl_Complex_Series.Link_to_Series
          := PentDobl_Complex_Series.Create(0,ord);
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
      PentDobl_CSeries_Polynomials.Add(res,rtm);
      PentDobl_CSeries_Polynomials.Clear(rtm.dg);
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new PentDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series_Polynomial;

  function Polynomial_to_Series_Polynomial
             ( p : OctoDobl_Complex_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return OctoDobl_CSeries_Polynomials.Poly is

    use OctoDobl_Complex_Series;

    res : OctoDobl_CSeries_Polynomials.Poly
        := OctoDobl_CSeries_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in OctoDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Adds the information in the term of a multivariate polynomial
    --   as a term of a series polynomial.  If idx = 0, then the degree of
    --   the series coefficient is zero and the coefficient is copied.
    --   Otherwise, the variable with index idx in t is taken as the
    --   variable in the truncated power series.  The other variables
    --   in t are relocated to fit a polynomial with one variable less.

      rtm : OctoDobl_CSeries_Polynomials.Term;
      ord : constant integer32
          := Set_Degree(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      dim : constant integer32
          := Set_Dimension(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      rcf : constant OctoDobl_Complex_Series.Link_to_Series
          := OctoDobl_Complex_Series.Create(0,ord);
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
      OctoDobl_CSeries_Polynomials.Add(res,rtm);
      OctoDobl_CSeries_Polynomials.Clear(rtm.dg);
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new OctoDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series_Polynomial;

  function Polynomial_to_Series_Polynomial
             ( p : DecaDobl_Complex_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return DecaDobl_CSeries_Polynomials.Poly is

    use DecaDobl_Complex_Series;

    res : DecaDobl_CSeries_Polynomials.Poly
        := DecaDobl_CSeries_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in DecaDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Adds the information in the term of a multivariate polynomial
    --   as a term of a series polynomial.  If idx = 0, then the degree of
    --   the series coefficient is zero and the coefficient is copied.
    --   Otherwise, the variable with index idx in t is taken as the
    --   variable in the truncated power series.  The other variables
    --   in t are relocated to fit a polynomial with one variable less.

      rtm : DecaDobl_CSeries_Polynomials.Term;
      ord : constant integer32
          := Set_Degree(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      dim : constant integer32
          := Set_Dimension(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      rcf : constant DecaDobl_Complex_Series.Link_to_Series
          := DecaDobl_Complex_Series.Create(0,ord);
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
      DecaDobl_CSeries_Polynomials.Add(res,rtm);
      DecaDobl_CSeries_Polynomials.Clear(rtm.dg);
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new DecaDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series_Polynomial;

  function Polynomial_to_Series_Polynomial
             ( p : HexaDobl_Complex_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return HexaDobl_CSeries_Polynomials.Poly is

    use HexaDobl_Complex_Series;

    res : HexaDobl_CSeries_Polynomials.Poly
        := HexaDobl_CSeries_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in HexaDobl_Complex_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Adds the information in the term of a multivariate polynomial
    --   as a term of a series polynomial.  If idx = 0, then the degree of
    --   the series coefficient is zero and the coefficient is copied.
    --   Otherwise, the variable with index idx in t is taken as the
    --   variable in the truncated power series.  The other variables
    --   in t are relocated to fit a polynomial with one variable less.

      rtm : HexaDobl_CSeries_Polynomials.Term;
      ord : constant integer32
          := Set_Degree(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      dim : constant integer32
          := Set_Dimension(idx,Standard_Natural_Vectors.Link_to_Vector(t.dg));
      rcf : constant HexaDobl_Complex_Series.Link_to_Series
          := HexaDobl_Complex_Series.Create(0,ord);
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
      HexaDobl_CSeries_Polynomials.Add(res,rtm);
      HexaDobl_CSeries_Polynomials.Clear(rtm.dg);
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new HexaDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Polynomial_to_Series_Polynomial;

  function Series_Polynomial_to_Polynomial
             ( s : Standard_CSeries_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Numbers;
    use Standard_Complex_Series;

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in Standard_CSeries_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Every term in the coefficient of the series in t.cf
    --   with a nonzero coefficient will contribute one monomial
    --   to the result.

      cffs : constant Link_to_Series := t.cf;
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
      new Standard_CSeries_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
    return res;
  end Series_Polynomial_to_Polynomial;

  function Series_Polynomial_to_Polynomial
             ( s : DoblDobl_CSeries_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Series;

    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in DoblDobl_CSeries_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Every term in the coefficient of the series in t.cf
    --   with a nonzero coefficient will contribute one monomial
    --   to the result.

      cffs : constant Link_to_Series := t.cf;
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
      new DoblDobl_CSeries_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
    return res;
  end Series_Polynomial_to_Polynomial;

  function Series_Polynomial_to_Polynomial
             ( s : TripDobl_CSeries_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return TripDobl_Complex_Polynomials.Poly is

    use TripDobl_Complex_Numbers;
    use TripDobl_Complex_Series;

    res : TripDobl_Complex_Polynomials.Poly
        := TripDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in TripDobl_CSeries_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Every term in the coefficient of the series in t.cf
    --   with a nonzero coefficient will contribute one monomial
    --   to the result.

      cffs : constant Link_to_Series := t.cf;
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
          rt : TripDobl_Complex_Polynomials.Term;
        begin
          rt.cf := cffs.cff(0);
          rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
          TripDobl_Complex_Polynomials.Add(res,rt);
          TripDobl_Complex_Polynomials.Clear(rt.dg);
        end;
      else -- idx > 0
        dim1 := t.dg'last+1;
        for k in 0..cffs.deg loop
          rtpc := cffs.cff(k);
          if not Equal(rtpc,zero) then
            declare
              rt : TripDobl_Complex_Polynomials.Term;
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
              TripDobl_Complex_Polynomials.Add(res,rt);
              TripDobl_Complex_Polynomials.Clear(rt.dg);
            end;
          end if;
        end loop;
      end if;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new TripDobl_CSeries_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
    return res;
  end Series_Polynomial_to_Polynomial;

  function Series_Polynomial_to_Polynomial
             ( s : QuadDobl_CSeries_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Series;

    res : QuadDobl_Complex_Polynomials.Poly
        := QuadDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in QuadDobl_CSeries_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Every term in the coefficient of the series in t.cf
    --   with a nonzero coefficient will contribute one monomial
    --   to the result.

      cffs : constant Link_to_Series := t.cf;
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
      new QuadDobl_CSeries_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
    return res;
  end Series_Polynomial_to_Polynomial;

  function Series_Polynomial_to_Polynomial
             ( s : PentDobl_CSeries_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return PentDobl_Complex_Polynomials.Poly is

    use PentDobl_Complex_Numbers;
    use PentDobl_Complex_Series;

    res : PentDobl_Complex_Polynomials.Poly
        := PentDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in PentDobl_CSeries_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Every term in the coefficient of the series in t.cf
    --   with a nonzero coefficient will contribute one monomial
    --   to the result.

      cffs : constant Link_to_Series := t.cf;
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
          rt : PentDobl_Complex_Polynomials.Term;
        begin
          rt.cf := cffs.cff(0);
          rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
          PentDobl_Complex_Polynomials.Add(res,rt);
          PentDobl_Complex_Polynomials.Clear(rt.dg);
        end;
      else -- idx > 0
        dim1 := t.dg'last+1;
        for k in 0..cffs.deg loop
          rtpc := cffs.cff(k);
          if not Equal(rtpc,zero) then
            declare
              rt : PentDobl_Complex_Polynomials.Term;
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
              PentDobl_Complex_Polynomials.Add(res,rt);
              PentDobl_Complex_Polynomials.Clear(rt.dg);
            end;
          end if;
        end loop;
      end if;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new PentDobl_CSeries_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
    return res;
  end Series_Polynomial_to_Polynomial;

  function Series_Polynomial_to_Polynomial
             ( s : OctoDobl_CSeries_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return OctoDobl_Complex_Polynomials.Poly is

    use OctoDobl_Complex_Numbers;
    use OctoDobl_Complex_Series;

    res : OctoDobl_Complex_Polynomials.Poly
        := OctoDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in OctoDobl_CSeries_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Every term in the coefficient of the series in t.cf
    --   with a nonzero coefficient will contribute one monomial
    --   to the result.

      cffs : constant Link_to_Series := t.cf;
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
          rt : OctoDobl_Complex_Polynomials.Term;
        begin
          rt.cf := cffs.cff(0);
          rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
          OctoDobl_Complex_Polynomials.Add(res,rt);
          OctoDobl_Complex_Polynomials.Clear(rt.dg);
        end;
      else -- idx > 0
        dim1 := t.dg'last+1;
        for k in 0..cffs.deg loop
          rtpc := cffs.cff(k);
          if not Equal(rtpc,zero) then
            declare
              rt : OctoDobl_Complex_Polynomials.Term;
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
              OctoDobl_Complex_Polynomials.Add(res,rt);
              OctoDobl_Complex_Polynomials.Clear(rt.dg);
            end;
          end if;
        end loop;
      end if;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new OctoDobl_CSeries_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
    return res;
  end Series_Polynomial_to_Polynomial;

  function Series_Polynomial_to_Polynomial
             ( s : DecaDobl_CSeries_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return DecaDobl_Complex_Polynomials.Poly is

    use DecaDobl_Complex_Numbers;
    use DecaDobl_Complex_Series;

    res : DecaDobl_Complex_Polynomials.Poly
        := DecaDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in DecaDobl_CSeries_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Every term in the coefficient of the series in t.cf
    --   with a nonzero coefficient will contribute one monomial
    --   to the result.

      cffs : constant Link_to_Series := t.cf;
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
          rt : DecaDobl_Complex_Polynomials.Term;
        begin
          rt.cf := cffs.cff(0);
          rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
          DecaDobl_Complex_Polynomials.Add(res,rt);
          DecaDobl_Complex_Polynomials.Clear(rt.dg);
        end;
      else -- idx > 0
        dim1 := t.dg'last+1;
        for k in 0..cffs.deg loop
          rtpc := cffs.cff(k);
          if not Equal(rtpc,zero) then
            declare
              rt : DecaDobl_Complex_Polynomials.Term;
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
              DecaDobl_Complex_Polynomials.Add(res,rt);
              DecaDobl_Complex_Polynomials.Clear(rt.dg);
            end;
          end if;
        end loop;
      end if;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new DecaDobl_CSeries_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
    return res;
  end Series_Polynomial_to_Polynomial;

  function Series_Polynomial_to_Polynomial
             ( s : HexaDobl_CSeries_Polynomials.Poly;
               idx : integer32 := 0; verbose : boolean := false )
             return HexaDobl_Complex_Polynomials.Poly is

    use HexaDobl_Complex_Numbers;
    use HexaDobl_Complex_Series;

    res : HexaDobl_Complex_Polynomials.Poly
        := HexaDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in HexaDobl_CSeries_Polynomials.Term;
                           c : out boolean ) is

    -- DESCRIPTION :
    --   Every term in the coefficient of the series in t.cf
    --   with a nonzero coefficient will contribute one monomial
    --   to the result.

      cffs : constant Link_to_Series := t.cf;
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
          rt : HexaDobl_Complex_Polynomials.Term;
        begin
          rt.cf := cffs.cff(0);
          rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
          HexaDobl_Complex_Polynomials.Add(res,rt);
          HexaDobl_Complex_Polynomials.Clear(rt.dg);
        end;
      else -- idx > 0
        dim1 := t.dg'last+1;
        for k in 0..cffs.deg loop
          rtpc := cffs.cff(k);
          if not Equal(rtpc,zero) then
            declare
              rt : HexaDobl_Complex_Polynomials.Term;
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
              HexaDobl_Complex_Polynomials.Add(res,rt);
              HexaDobl_Complex_Polynomials.Clear(rt.dg);
            end;
          end if;
        end loop;
      end if;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new HexaDobl_CSeries_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
    return res;
  end Series_Polynomial_to_Polynomial;

  function System_to_Series_System
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return Standard_CSeries_Poly_Systems.Poly_Sys is

    res : Standard_CSeries_Poly_Systems.Poly_Sys(p'range);

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
             return DoblDobl_CSeries_Poly_Systems.Poly_Sys is

    res : DoblDobl_CSeries_Poly_Systems.Poly_Sys(p'range);

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
             ( p : TripDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return TripDobl_CSeries_Poly_Systems.Poly_Sys is

    res : TripDobl_CSeries_Poly_Systems.Poly_Sys(p'range);

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
             return QuadDobl_CSeries_Poly_Systems.Poly_Sys is

    res : QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range);

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
             ( p : PentDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return PentDobl_CSeries_Poly_Systems.Poly_Sys is

    res : PentDobl_CSeries_Poly_Systems.Poly_Sys(p'range);

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
             ( p : OctoDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return OctoDobl_CSeries_Poly_Systems.Poly_Sys is

    res : OctoDobl_CSeries_Poly_Systems.Poly_Sys(p'range);

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
             ( p : DecaDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return DecaDobl_CSeries_Poly_Systems.Poly_Sys is

    res : DecaDobl_CSeries_Poly_Systems.Poly_Sys(p'range);

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
             ( p : HexaDobl_Complex_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return HexaDobl_CSeries_Poly_Systems.Poly_Sys is

    res : HexaDobl_CSeries_Poly_Systems.Poly_Sys(p'range);

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
             ( s : Standard_CSeries_Poly_Systems.Poly_Sys;
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
             ( s : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
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
             ( s : TripDobl_CSeries_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return TripDobl_Complex_Poly_Systems.Poly_Sys is

    res : TripDobl_Complex_Poly_Systems.Poly_Sys(s'range);

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
             ( s : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
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

  function Series_System_to_System
             ( s : PentDobl_CSeries_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return PentDobl_Complex_Poly_Systems.Poly_Sys is

    res : PentDobl_Complex_Poly_Systems.Poly_Sys(s'range);

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
             ( s : OctoDobl_CSeries_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys is

    res : OctoDobl_Complex_Poly_Systems.Poly_Sys(s'range);

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
             ( s : DecaDobl_CSeries_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys is

    res : DecaDobl_Complex_Poly_Systems.Poly_Sys(s'range);

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
             ( s : HexaDobl_CSeries_Poly_Systems.Poly_Sys;
               idx : integer32 := 0; verbose : boolean := false )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys is

    res : HexaDobl_Complex_Poly_Systems.Poly_Sys(s'range);

  begin
    for i in s'range loop
      if verbose
       then put("converting series polynomial "); put(i,1); put_line("...");
      end if;
      res(i) := Series_Polynomial_to_Polynomial(s(i),idx,verbose);
    end loop;
    return res;
  end Series_System_to_System;

  procedure Set_Degree ( v : in out Standard_Complex_Series_Vectors.Vector;
                         degree : in integer32 ) is
  begin
    for i in v'range loop
      Standard_Complex_Series.Set_Degree(v(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( v : in out DoblDobl_Complex_Series_Vectors.Vector;
                         degree : in integer32 ) is
  begin
    for i in v'range loop
      DoblDobl_Complex_Series.Set_Degree(v(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( v : in out TripDobl_Complex_Series_Vectors.Vector;
                         degree : in integer32 ) is
  begin
    for i in v'range loop
      TripDobl_Complex_Series.Set_Degree(v(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( v : in out QuadDobl_Complex_Series_Vectors.Vector;
                         degree : in integer32 ) is
  begin
    for i in v'range loop
      QuadDobl_Complex_Series.Set_Degree(v(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( v : in out PentDobl_Complex_Series_Vectors.Vector;
                         degree : in integer32 ) is
  begin
    for i in v'range loop
      PentDobl_Complex_Series.Set_Degree(v(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( v : in out OctoDobl_Complex_Series_Vectors.Vector;
                         degree : in integer32 ) is
  begin
    for i in v'range loop
      OctoDobl_Complex_Series.Set_Degree(v(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( v : in out DecaDobl_Complex_Series_Vectors.Vector;
                         degree : in integer32 ) is
  begin
    for i in v'range loop
      DecaDobl_Complex_Series.Set_Degree(v(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( v : in out HexaDobl_Complex_Series_Vectors.Vector;
                         degree : in integer32 ) is
  begin
    for i in v'range loop
      HexaDobl_Complex_Series.Set_Degree(v(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( m : in out Standard_Complex_Series_Matrices.Matrix;
                         degree : in integer32 ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        Standard_Complex_Series.Set_Degree(m(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( m : in out DoblDobl_Complex_Series_Matrices.Matrix;
                         degree : in integer32 ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        DoblDobl_Complex_Series.Set_Degree(m(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( m : in out TripDobl_Complex_Series_Matrices.Matrix;
                         degree : in integer32 ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        TripDobl_Complex_Series.Set_Degree(m(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( m : in out QuadDobl_Complex_Series_Matrices.Matrix;
                         degree : in integer32 ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        QuadDobl_Complex_Series.Set_Degree(m(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( m : in out PentDobl_Complex_Series_Matrices.Matrix;
                         degree : in integer32 ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        PentDobl_Complex_Series.Set_Degree(m(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( m : in out OctoDobl_Complex_Series_Matrices.Matrix;
                         degree : in integer32 ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        OctoDobl_Complex_Series.Set_Degree(m(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( m : in out DecaDobl_Complex_Series_Matrices.Matrix;
                         degree : in integer32 ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        DecaDobl_Complex_Series.Set_Degree(m(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( m : in out HexaDobl_Complex_Series_Matrices.Matrix;
                         degree : in integer32 ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        HexaDobl_Complex_Series.Set_Degree(m(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( p : in out Standard_CSeries_Polynomials.Poly;
                         degree : in integer32 ) is

    procedure Change_Degree ( t : in out Standard_CSeries_Polynomials.Term;
                              c : out boolean ) is
    begin
      Standard_Complex_Series.Set_Degree(t.cf,degree);
      c := true;
    end Change_Degree;
    procedure Change_Degrees is
      new Standard_CSeries_Polynomials.Changing_Iterator(Change_Degree);

  begin
    Change_Degrees(p);
  end Set_Degree;

  procedure Set_Degree ( p : in out DoblDobl_CSeries_Polynomials.Poly;
                         degree : in integer32 ) is

    procedure Change_Degree ( t : in out DoblDobl_CSeries_Polynomials.Term;
                              c : out boolean ) is
    begin
      DoblDobl_Complex_Series.Set_Degree(t.cf,degree);
      c := true;
    end Change_Degree;
    procedure Change_Degrees is
      new DoblDobl_CSeries_Polynomials.Changing_Iterator(Change_Degree);

  begin
    Change_Degrees(p);
  end Set_Degree;

  procedure Set_Degree ( p : in out TripDobl_CSeries_Polynomials.Poly;
                         degree : in integer32 ) is

    procedure Change_Degree ( t : in out TripDobl_CSeries_Polynomials.Term;
                              c : out boolean ) is
    begin
      TripDobl_Complex_Series.Set_Degree(t.cf,degree);
      c := true;
    end Change_Degree;
    procedure Change_Degrees is
      new TripDobl_CSeries_Polynomials.Changing_Iterator(Change_Degree);

  begin
    Change_Degrees(p);
  end Set_Degree;

  procedure Set_Degree ( p : in out QuadDobl_CSeries_Polynomials.Poly;
                        degree : in integer32 ) is

    procedure Change_Degree ( t : in out QuadDobl_CSeries_Polynomials.Term;
                              c : out boolean ) is
    begin
      QuadDobl_Complex_Series.Set_Degree(t.cf,degree);
      c := true;
    end Change_Degree;
    procedure Change_Degrees is
      new QuadDobl_CSeries_Polynomials.Changing_Iterator(Change_Degree);

  begin
    Change_Degrees(p);
  end Set_Degree;

  procedure Set_Degree ( p : in out PentDobl_CSeries_Polynomials.Poly;
                        degree : in integer32 ) is

    procedure Change_Degree ( t : in out PentDobl_CSeries_Polynomials.Term;
                              c : out boolean ) is
    begin
      PentDobl_Complex_Series.Set_Degree(t.cf,degree);
      c := true;
    end Change_Degree;
    procedure Change_Degrees is
      new PentDobl_CSeries_Polynomials.Changing_Iterator(Change_Degree);

  begin
    Change_Degrees(p);
  end Set_Degree;

  procedure Set_Degree ( p : in out OctoDobl_CSeries_Polynomials.Poly;
                        degree : in integer32 ) is

    procedure Change_Degree ( t : in out OctoDobl_CSeries_Polynomials.Term;
                              c : out boolean ) is
    begin
      OctoDobl_Complex_Series.Set_Degree(t.cf,degree);
      c := true;
    end Change_Degree;
    procedure Change_Degrees is
      new OctoDobl_CSeries_Polynomials.Changing_Iterator(Change_Degree);

  begin
    Change_Degrees(p);
  end Set_Degree;

  procedure Set_Degree ( p : in out DecaDobl_CSeries_Polynomials.Poly;
                        degree : in integer32 ) is

    procedure Change_Degree ( t : in out DecaDobl_CSeries_Polynomials.Term;
                              c : out boolean ) is
    begin
      DecaDobl_Complex_Series.Set_Degree(t.cf,degree);
      c := true;
    end Change_Degree;
    procedure Change_Degrees is
      new DecaDobl_CSeries_Polynomials.Changing_Iterator(Change_Degree);

  begin
    Change_Degrees(p);
  end Set_Degree;

  procedure Set_Degree ( p : in out HexaDobl_CSeries_Polynomials.Poly;
                        degree : in integer32 ) is

    procedure Change_Degree ( t : in out HexaDobl_CSeries_Polynomials.Term;
                              c : out boolean ) is
    begin
      HexaDobl_Complex_Series.Set_Degree(t.cf,degree);
      c := true;
    end Change_Degree;
    procedure Change_Degrees is
      new HexaDobl_CSeries_Polynomials.Changing_Iterator(Change_Degree);

  begin
    Change_Degrees(p);
  end Set_Degree;

  procedure Set_Degree ( p : in out Standard_CSeries_Poly_Systems.Poly_Sys;
                         degree : in integer32 ) is
  begin
    for i in p'range loop
      Set_Degree(p(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( p : in out DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                         degree : in integer32 ) is
  begin
    for i in p'range loop
      Set_Degree(p(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( p : in out TripDobl_CSeries_Poly_Systems.Poly_Sys;
                         degree : in integer32 ) is
  begin
    for i in p'range loop
      Set_Degree(p(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( p : in out QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                         degree : in integer32 ) is
  begin
    for i in p'range loop
      Set_Degree(p(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( p : in out PentDobl_CSeries_Poly_Systems.Poly_Sys;
                         degree : in integer32 ) is
  begin
    for i in p'range loop
      Set_Degree(p(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( p : in out OctoDobl_CSeries_Poly_Systems.Poly_Sys;
                         degree : in integer32 ) is
  begin
    for i in p'range loop
      Set_Degree(p(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( p : in out DecaDobl_CSeries_Poly_Systems.Poly_Sys;
                         degree : in integer32 ) is
  begin
    for i in p'range loop
      Set_Degree(p(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( p : in out HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                         degree : in integer32 ) is
  begin
    for i in p'range loop
      Set_Degree(p(i),degree);
    end loop;
  end Set_Degree;

  procedure Set_Degree ( jm : in out Standard_CSeries_Jaco_Matrices.Jaco_Mat;
                         degree : in integer32 ) is
  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        Set_Degree(jm(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( jm : in out DoblDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                         degree : in integer32 ) is
  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        Set_Degree(jm(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( jm : in out TripDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                         degree : in integer32 ) is
  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        Set_Degree(jm(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( jm : in out QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                         degree : in integer32 ) is
  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        Set_Degree(jm(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( jm : in out PentDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                         degree : in integer32 ) is
  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        Set_Degree(jm(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( jm : in out OctoDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                         degree : in integer32 ) is
  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        Set_Degree(jm(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( jm : in out DecaDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                         degree : in integer32 ) is
  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        Set_Degree(jm(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Set_Degree ( jm : in out HexaDobl_CSeries_Jaco_Matrices.Jaco_Mat;
                         degree : in integer32 ) is
  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        Set_Degree(jm(i,j),degree);
      end loop;
    end loop;
  end Set_Degree;

  procedure Filter ( s : in out Standard_Complex_Series_Vectors.Vector;
                     tol : in double_float ) is
  begin
    for i in s'range loop
      Standard_Complex_Series_Functions.Filter(s(i),tol);
    end loop;
  end Filter;

  procedure Filter ( s : in out DoblDobl_Complex_Series_Vectors.Vector;
                     tol : in double_float ) is
  begin
    for i in s'range loop
      DoblDobl_Complex_Series_Functions.Filter(s(i),tol);
    end loop;
  end Filter;

  procedure Filter ( s : in out TripDobl_Complex_Series_Vectors.Vector;
                     tol : in double_float ) is
  begin
    for i in s'range loop
      TripDobl_Complex_Series_Functions.Filter(s(i),tol);
    end loop;
  end Filter;

  procedure Filter ( s : in out QuadDobl_Complex_Series_Vectors.Vector;
                     tol : in double_float ) is
  begin
    for i in s'range loop
      QuadDobl_Complex_Series_Functions.Filter(s(i),tol);
    end loop;
  end Filter;

  procedure Filter ( s : in out PentDobl_Complex_Series_Vectors.Vector;
                     tol : in double_float ) is
  begin
    for i in s'range loop
      PentDobl_Complex_Series_Functions.Filter(s(i),tol);
    end loop;
  end Filter;

  procedure Filter ( s : in out OctoDobl_Complex_Series_Vectors.Vector;
                     tol : in double_float ) is
  begin
    for i in s'range loop
      OctoDobl_Complex_Series_Functions.Filter(s(i),tol);
    end loop;
  end Filter;

  procedure Filter ( s : in out DecaDobl_Complex_Series_Vectors.Vector;
                     tol : in double_float ) is
  begin
    for i in s'range loop
      DecaDobl_Complex_Series_Functions.Filter(s(i),tol);
    end loop;
  end Filter;

  procedure Filter ( s : in out HexaDobl_Complex_Series_Vectors.Vector;
                     tol : in double_float ) is
  begin
    for i in s'range loop
      HexaDobl_Complex_Series_Functions.Filter(s(i),tol);
    end loop;
  end Filter;
 
end Complex_Series_and_Polynomials;
