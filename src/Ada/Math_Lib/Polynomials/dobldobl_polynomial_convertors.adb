with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers_cv;        use DoblDobl_Complex_Numbers_cv;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;
with Multprec_DoblDobl_Convertors;       use Multprec_DoblDobl_Convertors;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;

package body DoblDobl_Polynomial_Convertors is

  function Standard_Polynomial_to_Double_Double
             ( p : Standard_Complex_Polynomials.Poly )
             return Double_Double_Polynomials.Poly is
  
    res : Double_Double_Polynomials.Poly
        := Double_Double_Polynomials.Null_Poly;       

    use Standard_Complex_Polynomials;

    procedure Term_to_Double_Double_Term ( t : in Term; c : out boolean ) is

      rt : Double_Double_Polynomials.Term;
      cf : constant double_float := Standard_Complex_Numbers.REAL_PART(t.cf);
      dd : constant double_double := create(cf);

    begin
      rt.cf := dd;
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Double_Double_Polynomials.Add(res,rt);
      Double_Double_Polynomials.Clear(rt);
      c := true;
    end Term_to_Double_Double_Term;
    procedure P2DD is new Visiting_Iterator(Term_to_Double_Double_Term);

  begin
    P2DD(p);
    return res;
  end Standard_Polynomial_to_Double_Double;

  function Standard_Poly_Sys_to_Double_Double
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Double_Double_Poly_Systems.Poly_Sys is

    res : Double_Double_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Standard_Polynomial_to_Double_Double(p(i));
    end loop;
    return res;
  end Standard_Poly_Sys_to_Double_Double;

  function Multprec_Polynomial_to_Double_Double
             ( p : Multprec_Complex_Polynomials.Poly )
             return Double_Double_Polynomials.Poly is
  
    res : Double_Double_Polynomials.Poly
        := Double_Double_Polynomials.Null_Poly;       

    use Multprec_Complex_Polynomials;

    procedure Term_to_Double_Double_Term ( t : in Term; c : out boolean ) is

      rt : Double_Double_Polynomials.Term;
      cf : Floating_Number := Multprec_Complex_Numbers.REAL_PART(t.cf);
      dd : constant double_double := to_double_double(cf);

    begin
      rt.cf := dd;
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Double_Double_Polynomials.Add(res,rt);
      Double_Double_Polynomials.Clear(rt);
      Multprec_Floating_Numbers.Clear(cf);
      c := true;
    end Term_to_Double_Double_Term;
    procedure MP2DD is new Visiting_Iterator(Term_to_Double_Double_Term);

  begin
    MP2DD(p);
    return res;
  end Multprec_Polynomial_to_Double_Double;

  function Multprec_Poly_Sys_to_Double_Double
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return Double_Double_Poly_Systems.Poly_Sys is

    res : Double_Double_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Multprec_Polynomial_to_Double_Double(p(i));
    end loop;
    return res;
  end Multprec_Poly_Sys_to_Double_Double;

  function Standard_Polynomial_to_DoblDobl_Complex
             ( p : Standard_Complex_Polynomials.Poly )
             return DoblDobl_Complex_Polynomials.Poly is

    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Null_Poly;

    use Standard_Complex_Polynomials;

    procedure Term_to_DoblDobl_Complex_Term ( t : in Term; c : out boolean ) is

      rt : DoblDobl_Complex_Polynomials.Term;

    begin
      rt.cf := Standard_to_DoblDobl_Complex(t.cf);
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      DoblDobl_Complex_Polynomials.Add(res,rt);
      DoblDobl_Complex_Polynomials.Clear(rt);
      c := true;
    end Term_to_DoblDobl_Complex_Term;
    procedure P2DC is new Visiting_Iterator(Term_to_DoblDobl_Complex_Term);

  begin
    P2DC(p);
    return res;
  end Standard_Polynomial_to_DoblDobl_Complex;

  function Standard_Poly_Sys_to_DoblDobl_Complex
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Standard_Polynomial_to_DoblDobl_Complex(p(i));
    end loop;
    return res;
  end Standard_Poly_Sys_to_DoblDobl_Complex;

  function Multprec_Polynomial_to_DoblDobl_Complex
             ( p : Multprec_Complex_Polynomials.Poly )
             return DoblDobl_Complex_Polynomials.Poly is

    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Null_Poly;

    use Multprec_Complex_Polynomials;

    procedure Term_to_DoblDobl_Complex_Term ( t : in Term; c : out boolean ) is

      rt : DoblDobl_Complex_Polynomials.Term;

    begin
      rt.cf := Multprec_to_DoblDobl_Complex(t.cf);
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      DoblDobl_Complex_Polynomials.Add(res,rt);
      DoblDobl_Complex_Polynomials.Clear(rt);
      c := true;
    end Term_to_DoblDobl_Complex_Term;
    procedure MP2DC is new Visiting_Iterator(Term_to_DoblDobl_Complex_Term);

  begin
    MP2DC(p);
    return res;
  end Multprec_Polynomial_to_DoblDobl_Complex;

  function Multprec_Poly_Sys_to_DoblDobl_Complex
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Multprec_Polynomial_to_DoblDobl_Complex(p(i));
    end loop;
    return res;
  end Multprec_Poly_Sys_to_DoblDobl_Complex;

  function Standard_Laurential_to_DoblDobl_Complex
             ( p : Standard_Complex_Laurentials.Poly )
             return DoblDobl_Complex_Laurentials.Poly is

    res : DoblDobl_Complex_Laurentials.Poly
        := DoblDobl_Complex_Laurentials.Null_Poly;

    use Standard_Complex_Laurentials;

    procedure Term_to_DoblDobl_Complex_Term ( t : in Term; c : out boolean ) is

      rt : DoblDobl_Complex_Laurentials.Term;

    begin
      rt.cf := Standard_to_DoblDobl_Complex(t.cf);
      rt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
      DoblDobl_Complex_Laurentials.Add(res,rt);
      DoblDobl_Complex_Laurentials.Clear(rt);
      c := true;
    end Term_to_DoblDobl_Complex_Term;
    procedure L2DC is new Visiting_Iterator(Term_to_DoblDobl_Complex_Term);

  begin
    L2DC(p);
    return res;
  end Standard_Laurential_to_DoblDobl_Complex;

  function Standard_Laur_Sys_to_DoblDobl_Complex
             ( p : Standard_Complex_Laur_Systems.Laur_Sys )
             return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Standard_Laurential_to_DoblDobl_Complex(p(i));
    end loop;
    return res;
  end Standard_Laur_Sys_to_DoblDobl_Complex;

  function Multprec_Laurential_to_DoblDobl_Complex
             ( p : Multprec_Complex_Laurentials.Poly )
             return DoblDobl_Complex_Laurentials.Poly is

    res : DoblDobl_Complex_Laurentials.Poly
        := DoblDobl_Complex_Laurentials.Null_Poly;

    use Multprec_Complex_Laurentials;

    procedure Term_to_DoblDobl_Complex_Term ( t : in Term; c : out boolean ) is

      rt : DoblDobl_Complex_Laurentials.Term;

    begin
      rt.cf := Multprec_to_DoblDobl_Complex(t.cf);
      rt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
      DoblDobl_Complex_Laurentials.Add(res,rt);
      DoblDobl_Complex_Laurentials.Clear(rt);
      c := true;
    end Term_to_DoblDobl_Complex_Term;
    procedure ML2DC is new Visiting_Iterator(Term_to_DoblDobl_Complex_Term);

  begin
    ML2DC(p);
    return res;
  end Multprec_Laurential_to_DoblDobl_Complex;

  function Multprec_Laur_Sys_to_DoblDobl_Complex
             ( p : Multprec_Complex_Laur_Systems.Laur_Sys )
             return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Multprec_Laurential_to_DoblDobl_Complex(p(i));
    end loop;
    return res;
  end Multprec_Laur_Sys_to_DoblDobl_Complex;

  function Double_Double_to_Standard_Polynomial
             ( p : Double_Double_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    use Double_Double_Polynomials;

    procedure Term_to_Standard_Term ( t : in Term; c : out boolean ) is

      rt : Standard_Complex_Polynomials.Term;
      re : constant double_float := to_double(t.cf);

    begin
      rt.cf := Standard_Complex_Numbers.Create(re);
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Standard_Complex_Polynomials.Add(res,rt);
      Standard_Complex_Polynomials.Clear(rt);
      c := true;
    end Term_to_Standard_Term;
    procedure DD2P is new Visiting_Iterator(Term_to_Standard_Term);

  begin
    DD2P(p);
    return res;
  end Double_Double_to_Standard_Polynomial;

  function Double_Double_to_Standard_Poly_Sys
             ( p : Double_Double_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Double_Double_to_Standard_Polynomial(p(i));
    end loop;
    return res;
  end Double_Double_to_Standard_Poly_Sys;

  function Double_Double_to_Multprec_Polynomial
             ( p : Double_Double_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly is

    res : Multprec_Complex_Polynomials.Poly
        := Multprec_Complex_Polynomials.Null_Poly;

    use Double_Double_Polynomials;

    procedure Term_to_Multprec_Term ( t : in Term; c : out boolean ) is

      rt : Multprec_Complex_Polynomials.Term;
      re : Floating_Number := to_floating_number(t.cf);

    begin
      rt.cf := Multprec_Complex_Numbers.Create(re);
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Multprec_Complex_Polynomials.Add(res,rt);
      Multprec_Complex_Polynomials.Clear(rt);
      Multprec_Floating_Numbers.Clear(re);
      c := true;
    end Term_to_Multprec_Term;
    procedure DD2MP is new Visiting_Iterator(Term_to_Multprec_Term);

  begin
    DD2MP(p);
    return res;
  end Double_Double_to_Multprec_Polynomial;

  function Double_Double_to_Multprec_Poly_Sys
             ( p : Double_Double_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys is

    res : Multprec_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Double_Double_to_Multprec_Polynomial(p(i));
    end loop;
    return res;
  end Double_Double_to_Multprec_Poly_Sys;

  function DoblDobl_Complex_to_Standard_Polynomial
             ( p : DoblDobl_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    use DoblDobl_Complex_Polynomials;

    procedure Term_to_Standard_Term ( t : in Term; c : out boolean ) is

      rt : Standard_Complex_Polynomials.Term;

    begin
      rt.cf := DoblDobl_Complex_to_Standard(t.cf);
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Standard_Complex_Polynomials.Add(res,rt);
      Standard_Complex_Polynomials.Clear(rt);
      c := true;
    end Term_to_Standard_Term;
    procedure DC2P is new Visiting_Iterator(Term_to_Standard_Term);

  begin
    DC2P(p);
    return res;
  end DoblDobl_Complex_to_Standard_Polynomial;

  function DoblDobl_Complex_to_Standard_Poly_Sys
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := DoblDobl_Complex_to_Standard_Polynomial(p(i));
    end loop;
    return res;
  end DoblDobl_Complex_to_Standard_Poly_Sys;

  function DoblDobl_Complex_to_Multprec_Polynomial
             ( p : DoblDobl_Complex_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly is

    res : Multprec_Complex_Polynomials.Poly
        := Multprec_Complex_Polynomials.Null_Poly;

    use DoblDobl_Complex_Polynomials;

    procedure Term_to_Multprec_Term ( t : in Term; c : out boolean ) is

      rt : Multprec_Complex_Polynomials.Term;

    begin
      rt.cf := DoblDobl_Complex_to_Multprec(t.cf);
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Multprec_Complex_Polynomials.Add(res,rt);
      Multprec_Complex_Polynomials.Clear(rt);
      c := true;
    end Term_to_Multprec_Term;
    procedure DC2MP is new Visiting_Iterator(Term_to_Multprec_Term);

  begin
    DC2MP(p);
    return res;
  end DoblDobl_Complex_to_Multprec_Polynomial;

  function DoblDobl_Complex_to_Multprec_Poly_Sys
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys is

    res : Multprec_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := DoblDobl_Complex_to_Multprec_Polynomial(p(i));
    end loop;
    return res;
  end DoblDobl_Complex_to_Multprec_Poly_Sys;

  function DoblDobl_Complex_to_Standard_Laurential
             ( p : DoblDobl_Complex_Laurentials.Poly )
             return Standard_Complex_Laurentials.Poly is

    res : Standard_Complex_Laurentials.Poly
        := Standard_Complex_Laurentials.Null_Poly;

    use DoblDobl_Complex_Laurentials;

    procedure Term_to_Standard_Term ( t : in Term; c : out boolean ) is

      rt : Standard_Complex_Laurentials.Term;

    begin
      rt.cf := DoblDobl_Complex_to_Standard(t.cf);
      rt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
      Standard_Complex_Laurentials.Add(res,rt);
      Standard_Complex_Laurentials.Clear(rt);
      c := true;
    end Term_to_Standard_Term;
    procedure DC2L is new Visiting_Iterator(Term_to_Standard_Term);

  begin
    DC2L(p);
    return res;
  end DoblDobl_Complex_to_Standard_Laurential;

  function DoblDobl_Complex_to_Standard_Laur_Sys
             ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
             return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := DoblDobl_Complex_to_Standard_Laurential(p(i));
    end loop;
    return res;
  end DoblDobl_Complex_to_Standard_Laur_Sys;

  function DoblDobl_Complex_to_Multprec_Laurential
             ( p : DoblDobl_Complex_Laurentials.Poly )
             return Multprec_Complex_Laurentials.Poly is

    res : Multprec_Complex_Laurentials.Poly
        := Multprec_Complex_Laurentials.Null_Poly;

    use DoblDobl_Complex_Laurentials;

    procedure Term_to_Multprec_Term ( t : in Term; c : out boolean ) is

      rt : Multprec_Complex_Laurentials.Term;

    begin
      rt.cf := DoblDobl_Complex_to_Multprec(t.cf);
      rt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
      Multprec_Complex_Laurentials.Add(res,rt);
      Multprec_Complex_Laurentials.Clear(rt);
      c := true;
    end Term_to_Multprec_Term;
    procedure DC2ML is new Visiting_Iterator(Term_to_Multprec_Term);

  begin
    DC2ML(p);
    return res;
  end DoblDobl_Complex_to_Multprec_Laurential;

  function DoblDobl_Complex_to_Multprec_Laur_Sys
             ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
             return Multprec_Complex_Laur_Systems.Laur_Sys is

    res : Multprec_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := DoblDobl_Complex_to_Multprec_Laurential(p(i));
    end loop;
    return res;
  end DoblDobl_Complex_to_Multprec_Laur_Sys;

end DoblDobl_Polynomial_Convertors;
