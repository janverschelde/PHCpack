with DecaDobl_Random_Numbers;
with Standard_Natural_Vectors;

package body DecaDobl_Complex_Poly_Randomizers is

  function Complex_Randomize ( p : Poly ) return Poly is

    res : Poly := Null_Poly;

    procedure Randomize_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      rt.cf := DecaDobl_Random_Numbers.Random;
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Add(res,rt);
      Clear(rt);
      cont := true;
    end Randomize_Term;
    procedure Randomize_Terms is new Visiting_Iterator (Randomize_Term);

  begin
    Randomize_Terms(p);
    return res;
  end Complex_Randomize;

  function Complex_Randomize1 ( p : Poly ) return Poly is

    res : Poly := Null_Poly;

    procedure Randomize_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      rt.cf := DecaDobl_Random_Numbers.Random1;
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Add(res,rt);
      Clear(rt);
      cont := true;
    end Randomize_Term;
    procedure Randomize_Terms is new Visiting_Iterator (Randomize_Term);
  begin
    Randomize_Terms(p);
    return res;
  end Complex_Randomize1;

  function Complex_Randomize ( p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in res'range loop
      res(i) := Complex_Randomize(p(i));
    end loop;
    return res;
  end Complex_Randomize;

  function Complex_Randomize1 ( p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in res'range loop
      res(i) := Complex_Randomize1(p(i));
    end loop;
    return res;
  end Complex_Randomize1;

end DecaDobl_Complex_Poly_Randomizers;
