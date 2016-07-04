with Standard_Natural_Vectors;
with Standard_Dense_Series;
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

end Series_and_Homotopies;
