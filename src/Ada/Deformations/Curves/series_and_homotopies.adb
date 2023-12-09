with Standard_Natural_Vectors;
with Standard_Complex_Series;
with Standard_Complex_Series_Functions;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_Functions;
with TripDobl_Complex_Series;
with TripDobl_Complex_Series_Functions;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_Functions;
with PentDobl_Complex_Series;
with PentDobl_Complex_Series_Functions;
with OctoDobl_Complex_Series;
with OctoDobl_Complex_Series_Functions;
with DecaDobl_Complex_Series;
with DecaDobl_Complex_Series_Functions;
with HexaDobl_Complex_Series;
with HexaDobl_Complex_Series_Functions;
with Complex_Series_and_Polynomials;

package body Series_and_Homotopies is

  function Create ( h : in Standard_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return Standard_CSeries_Poly_Systems.Poly_Sys is

    res : constant Standard_CSeries_Poly_Systems.Poly_Sys
        := Complex_Series_and_Polynomials.System_to_Series_System
             (h,idx,verbose);

  begin
    return res;
  end Create;

  function Create ( h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return DoblDobl_CSeries_Poly_Systems.Poly_Sys is

    res : constant DoblDobl_CSeries_Poly_Systems.Poly_Sys
        := Complex_Series_and_Polynomials.System_to_Series_System
             (h,idx,verbose);

  begin
    return res;
  end Create;

  function Create ( h : in TripDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return TripDobl_CSeries_Poly_Systems.Poly_Sys is

    res : constant TripDobl_CSeries_Poly_Systems.Poly_Sys
        := Complex_Series_and_Polynomials.System_to_Series_System
             (h,idx,verbose);

  begin
    return res;
  end Create;

  function Create ( h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return QuadDobl_CSeries_Poly_Systems.Poly_Sys is

    res : constant QuadDobl_CSeries_Poly_Systems.Poly_Sys
        := Complex_Series_and_Polynomials.System_to_Series_System
             (h,idx,verbose);

  begin
    return res;
  end Create;

  function Create ( h : in PentDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return PentDobl_CSeries_Poly_Systems.Poly_Sys is

    res : constant PentDobl_CSeries_Poly_Systems.Poly_Sys
        := Complex_Series_and_Polynomials.System_to_Series_System
             (h,idx,verbose);

  begin
    return res;
  end Create;

  function Create ( h : in OctoDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return OctoDobl_CSeries_Poly_Systems.Poly_Sys is

    res : constant OctoDobl_CSeries_Poly_Systems.Poly_Sys
        := Complex_Series_and_Polynomials.System_to_Series_System
             (h,idx,verbose);

  begin
    return res;
  end Create;

  function Create ( h : in DecaDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return DecaDobl_CSeries_Poly_Systems.Poly_Sys is

    res : constant DecaDobl_CSeries_Poly_Systems.Poly_Sys
        := Complex_Series_and_Polynomials.System_to_Series_System
             (h,idx,verbose);

  begin
    return res;
  end Create;

  function Create ( h : in HexaDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return HexaDobl_CSeries_Poly_Systems.Poly_Sys is

    res : constant HexaDobl_CSeries_Poly_Systems.Poly_Sys
        := Complex_Series_and_Polynomials.System_to_Series_System
             (h,idx,verbose);

  begin
    return res;
  end Create;

  function Eval ( p : Standard_CSeries_Polynomials.Poly;
                  t : double_float )
                return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    procedure Eval_Term ( trm : in Standard_CSeries_Polynomials.Term;
                          cont : out boolean ) is

      rt : Standard_Complex_Polynomials.Term;

    begin
      rt.cf := Standard_Complex_Series_Functions.Eval(trm.cf,t);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      Standard_Complex_Polynomials.Add(res,rt);
      Standard_Complex_Polynomials.Clear(rt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is
      new Standard_CSeries_Polynomials.Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( p : Standard_CSeries_Polynomials.Poly;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    procedure Eval_Term ( trm : in Standard_CSeries_Polynomials.Term;
                          cont : out boolean ) is

      rt : Standard_Complex_Polynomials.Term;

    begin
      rt.cf := Standard_Complex_Series_Functions.Eval(trm.cf,t);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      Standard_Complex_Polynomials.Add(res,rt);
      Standard_Complex_Polynomials.Clear(rt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is
      new Standard_CSeries_Polynomials.Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( p : DoblDobl_CSeries_Polynomials.Poly;
                  t : double_double )
                return DoblDobl_Complex_Polynomials.Poly is

    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Null_Poly;

    procedure Eval_Term ( trm : in DoblDobl_CSeries_Polynomials.Term;
                          cont : out boolean ) is

      rt : DoblDobl_Complex_Polynomials.Term;

    begin
      rt.cf := DoblDobl_Complex_Series_Functions.Eval(trm.cf,t);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      DoblDobl_Complex_Polynomials.Add(res,rt);
      DoblDobl_Complex_Polynomials.Clear(rt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is
      new DoblDobl_CSeries_Polynomials.Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( p : DoblDobl_CSeries_Polynomials.Poly;
                  t : DoblDobl_Complex_Numbers.Complex_Number )
                return DoblDobl_Complex_Polynomials.Poly is

    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Null_Poly;

    procedure Eval_Term ( trm : in DoblDobl_CSeries_Polynomials.Term;
                          cont : out boolean ) is

      rt : DoblDobl_Complex_Polynomials.Term;

    begin
      rt.cf := DoblDobl_Complex_Series_Functions.Eval(trm.cf,t);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      DoblDobl_Complex_Polynomials.Add(res,rt);
      DoblDobl_Complex_Polynomials.Clear(rt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is
      new DoblDobl_CSeries_Polynomials.Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( p : TripDobl_CSeries_Polynomials.Poly;
                  t : triple_double )
                return TripDobl_Complex_Polynomials.Poly is

    res : TripDobl_Complex_Polynomials.Poly
        := TripDobl_Complex_Polynomials.Null_Poly;

    procedure Eval_Term ( trm : in TripDobl_CSeries_Polynomials.Term;
                          cont : out boolean ) is

      rt : TripDobl_Complex_Polynomials.Term;

    begin
      rt.cf := TripDobl_Complex_Series_Functions.Eval(trm.cf,t);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      TripDobl_Complex_Polynomials.Add(res,rt);
      TripDobl_Complex_Polynomials.Clear(rt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is
      new TripDobl_CSeries_Polynomials.Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( p : TripDobl_CSeries_Polynomials.Poly;
                  t : TripDobl_Complex_Numbers.Complex_Number )
                return TripDobl_Complex_Polynomials.Poly is

    res : TripDobl_Complex_Polynomials.Poly
        := TripDobl_Complex_Polynomials.Null_Poly;

    procedure Eval_Term ( trm : in TripDobl_CSeries_Polynomials.Term;
                          cont : out boolean ) is

      rt : TripDobl_Complex_Polynomials.Term;

    begin
      rt.cf := TripDobl_Complex_Series_Functions.Eval(trm.cf,t);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      TripDobl_Complex_Polynomials.Add(res,rt);
      TripDobl_Complex_Polynomials.Clear(rt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is
      new TripDobl_CSeries_Polynomials.Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( p : QuadDobl_CSeries_Polynomials.Poly;
                  t : Quad_double )
                return QuadDobl_Complex_Polynomials.Poly is

    res : QuadDobl_Complex_Polynomials.Poly
        := QuadDobl_Complex_Polynomials.Null_Poly;

    procedure Eval_Term ( trm : in QuadDobl_CSeries_Polynomials.Term;
                          cont : out boolean ) is

      rt : QuadDobl_Complex_Polynomials.Term;

    begin
      rt.cf := QuadDobl_Complex_Series_Functions.Eval(trm.cf,t);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      QuadDobl_Complex_Polynomials.Add(res,rt);
      QuadDobl_Complex_Polynomials.Clear(rt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is
      new QuadDobl_CSeries_Polynomials.Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( p : QuadDobl_CSeries_Polynomials.Poly;
                  t : QuadDobl_Complex_Numbers.Complex_Number )
                return QuadDobl_Complex_Polynomials.Poly is

    res : QuadDobl_Complex_Polynomials.Poly
        := QuadDobl_Complex_Polynomials.Null_Poly;

    procedure Eval_Term ( trm : in QuadDobl_CSeries_Polynomials.Term;
                          cont : out boolean ) is

      rt : QuadDobl_Complex_Polynomials.Term;

    begin
      rt.cf := QuadDobl_Complex_Series_Functions.Eval(trm.cf,t);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      QuadDobl_Complex_Polynomials.Add(res,rt);
      QuadDobl_Complex_Polynomials.Clear(rt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is
      new QuadDobl_CSeries_Polynomials.Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( h : Standard_CSeries_Poly_Systems.Poly_Sys;
                  t : double_float )
                return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(h'range);

  begin
    for i in res'range loop
      res(i) := Eval(h(i),t);
    end loop;
    return res;
  end Eval;

  function Eval ( h : Standard_CSeries_Poly_Systems.Poly_Sys;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(h'range);

  begin
    for i in res'range loop
      res(i) := Eval(h(i),t);
    end loop;
    return res;
  end Eval;

  function Eval ( h : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                  t : double_double )
                return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(h'range);

  begin
    for i in res'range loop
      res(i) := Eval(h(i),t);
    end loop;
    return res;
  end Eval;

  function Eval ( h : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                  t : DoblDobl_Complex_Numbers.Complex_Number )
                return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(h'range);

  begin
    for i in res'range loop
      res(i) := Eval(h(i),t);
    end loop;
    return res;
  end Eval;

  function Eval ( h : TripDobl_CSeries_Poly_Systems.Poly_Sys;
                  t : triple_double )
                return TripDobl_Complex_Poly_Systems.Poly_Sys is

    res : TripDobl_Complex_Poly_Systems.Poly_Sys(h'range);

  begin
    for i in res'range loop
      res(i) := Eval(h(i),t);
    end loop;
    return res;
  end Eval;

  function Eval ( h : TripDobl_CSeries_Poly_Systems.Poly_Sys;
                  t : TripDobl_Complex_Numbers.Complex_Number )
                return TripDobl_Complex_Poly_Systems.Poly_Sys is

    res : TripDobl_Complex_Poly_Systems.Poly_Sys(h'range);

  begin
    for i in res'range loop
      res(i) := Eval(h(i),t);
    end loop;
    return res;
  end Eval;

  function Eval ( h : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                  t : quad_double )
                return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(h'range);

  begin
    for i in res'range loop
      res(i) := Eval(h(i),t);
    end loop;
    return res;
  end Eval;

  function Eval ( h : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                  t : QuadDobl_Complex_Numbers.Complex_Number )
                return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(h'range);

  begin
    for i in res'range loop
      res(i) := Eval(h(i),t);
    end loop;
    return res;
  end Eval;

  function Shift ( p : Standard_CSeries_Polynomials.Poly;
                   c : double_float )
                 return Standard_CSeries_Polynomials.Poly is

    res : Standard_CSeries_Polynomials.Poly
        := Standard_CSeries_Polynomials.Null_Poly;

    procedure Shift_Term ( trm : in Standard_CSeries_Polynomials.Term;
                           cont : out boolean ) is

      rt : Standard_CSeries_Polynomials.Term;

    begin
      rt.cf := Standard_Complex_Series_Functions.Shift(trm.cf,c);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      Standard_CSeries_Polynomials.Add(res,rt);
      Standard_CSeries_Polynomials.Clear(rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new Standard_CSeries_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
    return res;
  end Shift;

  function Shift ( p : Standard_CSeries_Polynomials.Poly;
                   c : Standard_Complex_Numbers.Complex_Number )
                 return Standard_CSeries_Polynomials.Poly is

    res : Standard_CSeries_Polynomials.Poly
        := Standard_CSeries_Polynomials.Null_Poly;

    procedure Shift_Term ( trm : in Standard_CSeries_Polynomials.Term;
                           cont : out boolean ) is

      rt : Standard_CSeries_Polynomials.Term;

    begin
      rt.cf := Standard_Complex_Series_Functions.Shift(trm.cf,c);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      Standard_CSeries_Polynomials.Add(res,rt);
      Standard_CSeries_Polynomials.Clear(rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new Standard_CSeries_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
    return res;
  end Shift;

  function Shift ( p : DoblDobl_CSeries_Polynomials.Poly;
                   c : double_double )
                 return DoblDobl_CSeries_Polynomials.Poly is

    res : DoblDobl_CSeries_Polynomials.Poly
        := DoblDobl_CSeries_Polynomials.Null_Poly;

    procedure Shift_Term ( trm : in DoblDobl_CSeries_Polynomials.Term;
                           cont : out boolean ) is

      rt : DoblDobl_CSeries_Polynomials.Term;

    begin
      rt.cf := DoblDobl_Complex_Series_Functions.Shift(trm.cf,c);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      DoblDobl_CSeries_Polynomials.Add(res,rt);
      DoblDobl_CSeries_Polynomials.Clear(rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new DoblDobl_CSeries_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
    return res;
  end Shift;

  function Shift ( p : DoblDobl_CSeries_Polynomials.Poly;
                   c : DoblDobl_Complex_Numbers.Complex_Number )
                 return DoblDobl_CSeries_Polynomials.Poly is

    res : DoblDobl_CSeries_Polynomials.Poly
        := DoblDobl_CSeries_Polynomials.Null_Poly;

    procedure Shift_Term ( trm : in DoblDobl_CSeries_Polynomials.Term;
                           cont : out boolean ) is

      rt : DoblDobl_CSeries_Polynomials.Term;

    begin
      rt.cf := DoblDobl_Complex_Series_Functions.Shift(trm.cf,c);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      DoblDobl_CSeries_Polynomials.Add(res,rt);
      DoblDobl_CSeries_Polynomials.Clear(rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new DoblDobl_CSeries_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
    return res;
  end Shift;

  function Shift ( p : QuadDobl_CSeries_Polynomials.Poly;
                   c : quad_double )
                 return QuadDobl_CSeries_Polynomials.Poly is

    res : QuadDobl_CSeries_Polynomials.Poly
        := QuadDobl_CSeries_Polynomials.Null_Poly;

    procedure Shift_Term ( trm : in QuadDobl_CSeries_Polynomials.Term;
                           cont : out boolean ) is

      rt : QuadDobl_CSeries_Polynomials.Term;

    begin
      rt.cf := QuadDobl_Complex_Series_Functions.Shift(trm.cf,c);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      QuadDobl_CSeries_Polynomials.Add(res,rt);
      QuadDobl_CSeries_Polynomials.Clear(rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new QuadDobl_CSeries_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
    return res;
  end Shift;

  function Shift ( p : QuadDobl_CSeries_Polynomials.Poly;
                   c : QuadDobl_Complex_Numbers.Complex_Number )
                 return QuadDobl_CSeries_Polynomials.Poly is

    res : QuadDobl_CSeries_Polynomials.Poly
        := QuadDobl_CSeries_Polynomials.Null_Poly;

    procedure Shift_Term ( trm : in QuadDobl_CSeries_Polynomials.Term;
                           cont : out boolean ) is

      rt : QuadDobl_CSeries_Polynomials.Term;

    begin
      rt.cf := QuadDobl_Complex_Series_Functions.Shift(trm.cf,c);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      QuadDobl_CSeries_Polynomials.Add(res,rt);
      QuadDobl_CSeries_Polynomials.Clear(rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new QuadDobl_CSeries_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
    return res;
  end Shift;

  procedure Shift ( p : in out Standard_CSeries_Polynomials.Poly;
                    c : in double_float ) is

    procedure Shift_Term ( trm : in Standard_CSeries_Polynomials.Term;
                           cont : out boolean ) is

      cff : constant Standard_Complex_Series.Link_to_Series := trm.cf;

    begin
      Standard_Complex_Series_Functions.Shift(cff,c);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new Standard_CSeries_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
  end Shift;

  procedure Shift ( p : in out Standard_CSeries_Polynomials.Poly;
                    c : in Standard_Complex_Numbers.Complex_Number ) is

    procedure Shift_Term ( trm : in Standard_CSeries_Polynomials.Term;
                           cont : out boolean ) is

      cff : constant Standard_Complex_Series.Link_to_Series := trm.cf;

    begin
      Standard_Complex_Series_Functions.Shift(cff,c);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new Standard_CSeries_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
  end Shift;

  procedure Shift ( p : in out DoblDobl_CSeries_Polynomials.Poly;
                    c : in double_double ) is

    procedure Shift_Term ( trm : in DoblDobl_CSeries_Polynomials.Term;
                           cont : out boolean ) is

      cff : constant DoblDobl_Complex_Series.Link_to_Series := trm.cf;

    begin
      DoblDobl_Complex_Series_Functions.Shift(cff,c);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new DoblDobl_CSeries_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
  end Shift;

  procedure Shift ( p : in out DoblDobl_CSeries_Polynomials.Poly;
                    c : in DoblDobl_Complex_Numbers.Complex_Number ) is

    procedure Shift_Term ( trm : in DoblDobl_CSeries_Polynomials.Term;
                           cont : out boolean ) is

      cff : DoblDobl_Complex_Series.Link_to_Series := trm.cf;

    begin
      DoblDobl_Complex_Series_Functions.Shift(cff,c);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new DoblDobl_CSeries_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
  end Shift;

  procedure Shift ( p : in out QuadDobl_CSeries_Polynomials.Poly;
                    c : in quad_double ) is

    procedure Shift_Term ( trm : in QuadDobl_CSeries_Polynomials.Term;
                           cont : out boolean ) is

      cff : constant QuadDobl_Complex_Series.Link_to_Series := trm.cf;

    begin
      QuadDobl_Complex_Series_Functions.Shift(cff,c);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new QuadDobl_CSeries_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
  end Shift;

  procedure Shift ( p : in out QuadDobl_CSeries_Polynomials.Poly;
                    c : in QuadDobl_Complex_Numbers.Complex_Number ) is

    procedure Shift_Term ( trm : in QuadDobl_CSeries_Polynomials.Term;
                           cont : out boolean ) is

      cff : constant QuadDobl_Complex_Series.Link_to_Series := trm.cf;

    begin
      QuadDobl_Complex_Series_Functions.Shift(cff,c);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new QuadDobl_CSeries_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
  end Shift;

  function Shift ( p : Standard_CSeries_Poly_Systems.Poly_Sys;
                   c : double_float )
                 return Standard_CSeries_Poly_Systems.Poly_Sys is

    res : Standard_CSeries_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Shift(p(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( p : Standard_CSeries_Poly_Systems.Poly_Sys;
                   c : Standard_Complex_Numbers.Complex_Number )
                 return Standard_CSeries_Poly_Systems.Poly_Sys is

    res : Standard_CSeries_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Shift(p(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( p : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                   c : double_double )
                 return DoblDobl_CSeries_Poly_Systems.Poly_Sys is

    res : DoblDobl_CSeries_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Shift(p(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( p : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                   c : DoblDobl_Complex_Numbers.Complex_Number )
                 return DoblDobl_CSeries_Poly_Systems.Poly_Sys is

    res : DoblDobl_CSeries_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Shift(p(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( p : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                   c : quad_double )
                 return QuadDobl_CSeries_Poly_Systems.Poly_Sys is

    res : QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Shift(p(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( p : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                   c : QuadDobl_Complex_Numbers.Complex_Number )
                 return QuadDobl_CSeries_Poly_Systems.Poly_Sys is

    res : QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Shift(p(i),c);
    end loop;
    return res;
  end Shift;

  procedure Shift ( p : in out Standard_CSeries_Poly_Systems.Poly_Sys;
                    c : in double_float ) is
  begin
    for i in p'range loop
      Shift(p(i),c);
    end loop;
  end Shift;

  procedure Shift ( p : in out Standard_CSeries_Poly_Systems.Poly_Sys;
                    c : in Standard_Complex_Numbers.Complex_Number ) is
  begin
    for i in p'range loop
      Shift(p(i),c);
    end loop;
  end Shift;

  procedure Shift ( p : in out DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                    c : in double_double ) is
  begin
    for i in p'range loop
      Shift(p(i),c);
    end loop;
  end Shift;

  procedure Shift ( p : in out DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                    c : in DoblDobl_Complex_Numbers.Complex_Number ) is
  begin
    for i in p'range loop
      Shift(p(i),c);
    end loop;
  end Shift;

  procedure Shift ( p : in out QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                    c : in quad_double ) is
  begin
    for i in p'range loop
      Shift(p(i),c);
    end loop;
  end Shift;

  procedure Shift ( p : in out QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                    c : in QuadDobl_Complex_Numbers.Complex_Number ) is
  begin
    for i in p'range loop
      Shift(p(i),c);
    end loop;
  end Shift;

end Series_and_Homotopies;
