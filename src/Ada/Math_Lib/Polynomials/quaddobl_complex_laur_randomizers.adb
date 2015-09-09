with QuadDobl_Random_Numbers;            use QuadDobl_Random_Numbers;
with Standard_Integer_Vectors;

package body QuadDobl_Complex_Laur_Randomizers is

  function Complex_Randomize ( p : Poly ) return Poly is

    res : Poly := Null_Poly;

    procedure Randomize_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      rt.cf := Random;
      rt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
      Add(res,rt);
      Clear(rt);
      cont := true;
    end Randomize_Term;
    procedure Randomize_Terms is new Visiting_Iterator(Randomize_Term);

  begin
    Randomize_Terms(p);
    return res;
  end Complex_Randomize;

  function Complex_Randomize1 ( p : Poly ) return Poly is

    res : Poly := Null_Poly;

    procedure Randomize_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      rt.cf := Random1;
      rt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
      Add(res,rt);
      Clear(rt);
      cont := true;
    end Randomize_Term;
    procedure Randomize_Terms is new Visiting_Iterator(Randomize_Term);

  begin
    Randomize_Terms(p);
    return res;
  end Complex_Randomize1;

  function Complex_Randomize ( p : Laur_Sys ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for i in res'range loop
      res(i) := Complex_Randomize(p(i));
    end loop;
    return res;
  end Complex_Randomize;

  function Complex_Randomize1 ( p : Laur_Sys ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for i in res'range loop
      res(i) := Complex_Randomize1(p(i));
    end loop;
    return res;
  end Complex_Randomize1;

end QuadDobl_Complex_Laur_Randomizers;
