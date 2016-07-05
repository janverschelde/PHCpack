with Standard_Natural_Vectors;
with Standard_Dense_Series;
with DoblDobl_Dense_Series;
with QuadDobl_Dense_Series;
with Series_and_Polynomials;

package body Series_and_Homotopies is

  function Create ( h : in Standard_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return Standard_Series_Poly_Systems.Poly_Sys is

    res : constant Standard_Series_Poly_Systems.Poly_Sys
        := Series_and_Polynomials.System_to_Series_System(h,idx,verbose);

  begin
    return res;
  end Create;

  function Create ( h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return DoblDobl_Series_Poly_Systems.Poly_Sys is

    res : constant DoblDobl_Series_Poly_Systems.Poly_Sys
        := Series_and_Polynomials.System_to_Series_System(h,idx,verbose);

  begin
    return res;
  end Create;

  function Create ( h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    idx : in integer32; verbose : boolean := false )
                  return QuadDobl_Series_Poly_Systems.Poly_Sys is

    res : constant QuadDobl_Series_Poly_Systems.Poly_Sys
        := Series_and_Polynomials.System_to_Series_System(h,idx,verbose);

  begin
    return res;
  end Create;

  function Eval ( p : Standard_Series_Polynomials.Poly;
                  t : double_float )
                return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    procedure Eval_Term ( trm : in Standard_Series_Polynomials.Term;
                          cont : out boolean ) is

      rt : Standard_Complex_Polynomials.Term;

    begin
      rt.cf := Standard_Dense_Series.Eval(trm.cf,t);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      Standard_Complex_Polynomials.Add(res,rt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is
      new Standard_Series_Polynomials.Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( p : Standard_Series_Polynomials.Poly;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    procedure Eval_Term ( trm : in Standard_Series_Polynomials.Term;
                          cont : out boolean ) is

      rt : Standard_Complex_Polynomials.Term;

    begin
      rt.cf := Standard_Dense_Series.Eval(trm.cf,t);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      Standard_Complex_Polynomials.Add(res,rt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is
      new Standard_Series_Polynomials.Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( p : DoblDobl_Series_Polynomials.Poly;
                  t : double_double )
                return DoblDobl_Complex_Polynomials.Poly is

    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Null_Poly;

    procedure Eval_Term ( trm : in DoblDobl_Series_Polynomials.Term;
                          cont : out boolean ) is

      rt : DoblDobl_Complex_Polynomials.Term;

    begin
      rt.cf := DoblDobl_Dense_Series.Eval(trm.cf,t);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      DoblDobl_Complex_Polynomials.Add(res,rt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is
      new DoblDobl_Series_Polynomials.Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( p : DoblDobl_Series_Polynomials.Poly;
                  t : DoblDobl_Complex_Numbers.Complex_Number )
                return DoblDobl_Complex_Polynomials.Poly is

    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Null_Poly;

    procedure Eval_Term ( trm : in DoblDobl_Series_Polynomials.Term;
                          cont : out boolean ) is

      rt : DoblDobl_Complex_Polynomials.Term;

    begin
      rt.cf := DoblDobl_Dense_Series.Eval(trm.cf,t);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      DoblDobl_Complex_Polynomials.Add(res,rt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is
      new DoblDobl_Series_Polynomials.Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( p : QuadDobl_Series_Polynomials.Poly;
                  t : Quad_double )
                return QuadDobl_Complex_Polynomials.Poly is

    res : QuadDobl_Complex_Polynomials.Poly
        := QuadDobl_Complex_Polynomials.Null_Poly;

    procedure Eval_Term ( trm : in QuadDobl_Series_Polynomials.Term;
                          cont : out boolean ) is

      rt : QuadDobl_Complex_Polynomials.Term;

    begin
      rt.cf := QuadDobl_Dense_Series.Eval(trm.cf,t);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      QuadDobl_Complex_Polynomials.Add(res,rt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is
      new QuadDobl_Series_Polynomials.Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( p : QuadDobl_Series_Polynomials.Poly;
                  t : QuadDobl_Complex_Numbers.Complex_Number )
                return QuadDobl_Complex_Polynomials.Poly is

    res : QuadDobl_Complex_Polynomials.Poly
        := QuadDobl_Complex_Polynomials.Null_Poly;

    procedure Eval_Term ( trm : in QuadDobl_Series_Polynomials.Term;
                          cont : out boolean ) is

      rt : QuadDobl_Complex_Polynomials.Term;

    begin
      rt.cf := QuadDobl_Dense_Series.Eval(trm.cf,t);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      QuadDobl_Complex_Polynomials.Add(res,rt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is
      new QuadDobl_Series_Polynomials.Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( h : Standard_Series_Poly_Systems.Poly_Sys;
                  t : double_float )
                return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(h'range);

  begin
    for i in res'range loop
      res(i) := Eval(h(i),t);
    end loop;
    return res;
  end Eval;

  function Eval ( h : Standard_Series_Poly_Systems.Poly_Sys;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(h'range);

  begin
    for i in res'range loop
      res(i) := Eval(h(i),t);
    end loop;
    return res;
  end Eval;

  function Eval ( h : DoblDobl_Series_Poly_Systems.Poly_Sys;
                  t : double_double )
                return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(h'range);

  begin
    for i in res'range loop
      res(i) := Eval(h(i),t);
    end loop;
    return res;
  end Eval;

  function Eval ( h : DoblDobl_Series_Poly_Systems.Poly_Sys;
                  t : DoblDobl_Complex_Numbers.Complex_Number )
                return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(h'range);

  begin
    for i in res'range loop
      res(i) := Eval(h(i),t);
    end loop;
    return res;
  end Eval;

  function Eval ( h : QuadDobl_Series_Poly_Systems.Poly_Sys;
                  t : quad_double )
                return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(h'range);

  begin
    for i in res'range loop
      res(i) := Eval(h(i),t);
    end loop;
    return res;
  end Eval;

  function Eval ( h : QuadDobl_Series_Poly_Systems.Poly_Sys;
                  t : QuadDobl_Complex_Numbers.Complex_Number )
                return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(h'range);

  begin
    for i in res'range loop
      res(i) := Eval(h(i),t);
    end loop;
    return res;
  end Eval;

  function Shift ( p : Standard_Series_Polynomials.Poly;
                   c : double_float )
                 return Standard_Series_Polynomials.Poly is

    res : Standard_Series_Polynomials.Poly
        := Standard_Series_Polynomials.Null_Poly;

    procedure Shift_Term ( trm : in Standard_Series_Polynomials.Term;
                           cont : out boolean ) is

      rt : Standard_Series_Polynomials.Term;

    begin
      rt.cf := Standard_Dense_Series.Shift(trm.cf,c);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      Standard_Series_Polynomials.Add(res,rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new Standard_Series_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
    return res;
  end Shift;

  function Shift ( p : Standard_Series_Polynomials.Poly;
                   c : Standard_Complex_Numbers.Complex_Number )
                 return Standard_Series_Polynomials.Poly is

    res : Standard_Series_Polynomials.Poly
        := Standard_Series_Polynomials.Null_Poly;

    procedure Shift_Term ( trm : in Standard_Series_Polynomials.Term;
                           cont : out boolean ) is

      rt : Standard_Series_Polynomials.Term;

    begin
      rt.cf := Standard_Dense_Series.Shift(trm.cf,c);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      Standard_Series_Polynomials.Add(res,rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new Standard_Series_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
    return res;
  end Shift;

  function Shift ( p : DoblDobl_Series_Polynomials.Poly;
                   c : double_double )
                 return DoblDobl_Series_Polynomials.Poly is

    res : DoblDobl_Series_Polynomials.Poly
        := DoblDobl_Series_Polynomials.Null_Poly;

    procedure Shift_Term ( trm : in DoblDobl_Series_Polynomials.Term;
                           cont : out boolean ) is

      rt : DoblDobl_Series_Polynomials.Term;

    begin
      rt.cf := DoblDobl_Dense_Series.Shift(trm.cf,c);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      DoblDobl_Series_Polynomials.Add(res,rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new DoblDobl_Series_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
    return res;
  end Shift;

  function Shift ( p : DoblDobl_Series_Polynomials.Poly;
                   c : DoblDobl_Complex_Numbers.Complex_Number )
                 return DoblDobl_Series_Polynomials.Poly is

    res : DoblDobl_Series_Polynomials.Poly
        := DoblDobl_Series_Polynomials.Null_Poly;

    procedure Shift_Term ( trm : in DoblDobl_Series_Polynomials.Term;
                           cont : out boolean ) is

      rt : DoblDobl_Series_Polynomials.Term;

    begin
      rt.cf := DoblDobl_Dense_Series.Shift(trm.cf,c);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      DoblDobl_Series_Polynomials.Add(res,rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new DoblDobl_Series_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
    return res;
  end Shift;

  function Shift ( p : QuadDobl_Series_Polynomials.Poly;
                   c : quad_double )
                 return QuadDobl_Series_Polynomials.Poly is

    res : QuadDobl_Series_Polynomials.Poly
        := QuadDobl_Series_Polynomials.Null_Poly;

    procedure Shift_Term ( trm : in QuadDobl_Series_Polynomials.Term;
                           cont : out boolean ) is

      rt : QuadDobl_Series_Polynomials.Term;

    begin
      rt.cf := QuadDobl_Dense_Series.Shift(trm.cf,c);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      QuadDobl_Series_Polynomials.Add(res,rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new QuadDobl_Series_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
    return res;
  end Shift;

  function Shift ( p : QuadDobl_Series_Polynomials.Poly;
                   c : QuadDobl_Complex_Numbers.Complex_Number )
                 return QuadDobl_Series_Polynomials.Poly is

    res : QuadDobl_Series_Polynomials.Poly
        := QuadDobl_Series_Polynomials.Null_Poly;

    procedure Shift_Term ( trm : in QuadDobl_Series_Polynomials.Term;
                           cont : out boolean ) is

      rt : QuadDobl_Series_Polynomials.Term;

    begin
      rt.cf := QuadDobl_Dense_Series.Shift(trm.cf,c);
      rt.dg := new Standard_Natural_Vectors.Vector(trm.dg'range);
      for k in rt.dg'range loop
        rt.dg(k) := trm.dg(k);
      end loop;
      QuadDobl_Series_Polynomials.Add(res,rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is
      new QuadDobl_Series_Polynomials.Visiting_Iterator(Shift_Term);

  begin
    Shift_Terms(p);
    return res;
  end Shift;

  procedure Shift ( p : in out Standard_Series_Polynomials.Poly;
                    c : in double_float ) is

    res : constant Standard_Series_Polynomials.Poly := Shift(p,c);

  begin
    Standard_Series_Polynomials.Clear(p);
    p := res;
  end Shift;

  procedure Shift ( p : in out Standard_Series_Polynomials.Poly;
                    c : in Standard_Complex_Numbers.Complex_Number ) is

    res : constant Standard_Series_Polynomials.Poly := Shift(p,c);

  begin
    Standard_Series_Polynomials.Clear(p);
    p := res;
  end Shift;

  procedure Shift ( p : in out DoblDobl_Series_Polynomials.Poly;
                    c : in double_double ) is

    res : constant DoblDobl_Series_Polynomials.Poly := Shift(p,c);

  begin
    DoblDobl_Series_Polynomials.Clear(p);
    p := res;
  end Shift;

  procedure Shift ( p : in out DoblDobl_Series_Polynomials.Poly;
                    c : in DoblDobl_Complex_Numbers.Complex_Number ) is

    res : constant DoblDobl_Series_Polynomials.Poly := Shift(p,c);

  begin
    DoblDobl_Series_Polynomials.Clear(p);
    p := res;
  end Shift;

  procedure Shift ( p : in out QuadDobl_Series_Polynomials.Poly;
                    c : in quad_double ) is

    res : constant QuadDobl_Series_Polynomials.Poly := Shift(p,c);

  begin
    QuadDobl_Series_Polynomials.Clear(p);
    p := res;
  end Shift;

  procedure Shift ( p : in out QuadDobl_Series_Polynomials.Poly;
                    c : in QuadDobl_Complex_Numbers.Complex_Number ) is

    res : constant QuadDobl_Series_Polynomials.Poly := Shift(p,c);

  begin
    QuadDobl_Series_Polynomials.Clear(p);
    p := res;
  end Shift;

  function Shift ( p : Standard_Series_Poly_Systems.Poly_Sys;
                   c : double_float )
                 return Standard_Series_Poly_Systems.Poly_Sys is

    res : Standard_Series_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Shift(p(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( p : Standard_Series_Poly_Systems.Poly_Sys;
                   c : Standard_Complex_Numbers.Complex_Number )
                 return Standard_Series_Poly_Systems.Poly_Sys is

    res : Standard_Series_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Shift(p(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( p : DoblDobl_Series_Poly_Systems.Poly_Sys;
                   c : double_double )
                 return DoblDobl_Series_Poly_Systems.Poly_Sys is

    res : DoblDobl_Series_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Shift(p(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( p : DoblDobl_Series_Poly_Systems.Poly_Sys;
                   c : DoblDobl_Complex_Numbers.Complex_Number )
                 return DoblDobl_Series_Poly_Systems.Poly_Sys is

    res : DoblDobl_Series_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Shift(p(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( p : QuadDobl_Series_Poly_Systems.Poly_Sys;
                   c : quad_double )
                 return QuadDobl_Series_Poly_Systems.Poly_Sys is

    res : QuadDobl_Series_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Shift(p(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( p : QuadDobl_Series_Poly_Systems.Poly_Sys;
                   c : QuadDobl_Complex_Numbers.Complex_Number )
                 return QuadDobl_Series_Poly_Systems.Poly_Sys is

    res : QuadDobl_Series_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Shift(p(i),c);
    end loop;
    return res;
  end Shift;

  procedure Shift ( p : in out Standard_Series_Poly_Systems.Poly_Sys;
                    c : in double_float ) is

    res : Standard_Series_Poly_Systems.Poly_Sys(p'range) := Shift(p,c);

  begin
    Standard_Series_Poly_Systems.Clear(p);
    p := res;
  end Shift;

  procedure Shift ( p : in out Standard_Series_Poly_Systems.Poly_Sys;
                    c : in Standard_Complex_Numbers.Complex_Number ) is

    res : Standard_Series_Poly_Systems.Poly_Sys(p'range) := Shift(p,c);

  begin
    Standard_Series_Poly_Systems.Clear(p);
    p := res;
  end Shift;

  procedure Shift ( p : in out DoblDobl_Series_Poly_Systems.Poly_Sys;
                    c : in double_double ) is

    res : DoblDobl_Series_Poly_Systems.Poly_Sys(p'range) := Shift(p,c);

  begin
    DoblDobl_Series_Poly_Systems.Clear(p);
    p := res;
  end Shift;

  procedure Shift ( p : in out DoblDobl_Series_Poly_Systems.Poly_Sys;
                    c : in DoblDobl_Complex_Numbers.Complex_Number ) is

    res : DoblDobl_Series_Poly_Systems.Poly_Sys(p'range) := Shift(p,c);

  begin
    DoblDobl_Series_Poly_Systems.Clear(p);
    p := res;
  end Shift;

  procedure Shift ( p : in out QuadDobl_Series_Poly_Systems.Poly_Sys;
                    c : in quad_double ) is

    res : QuadDobl_Series_Poly_Systems.Poly_Sys(p'range) := Shift(p,c);

  begin
    QuadDobl_Series_Poly_Systems.Clear(p);
    p := res;
  end Shift;

  procedure Shift ( p : in out QuadDobl_Series_Poly_Systems.Poly_Sys;
                    c : in QuadDobl_Complex_Numbers.Complex_Number ) is

    res : QuadDobl_Series_Poly_Systems.Poly_Sys(p'range) := Shift(p,c);

  begin
    QuadDobl_Series_Poly_Systems.Clear(p);
    p := res;
  end Shift;

end Series_and_Homotopies;
