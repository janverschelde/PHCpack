with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with HexaDobl_Complex_Numbers_cv;        use HexaDobl_Complex_Numbers_cv;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;
with Multprec_HexaDobl_Convertors;       use Multprec_HexaDobl_Convertors;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;

package body HexaDobl_Polynomial_Convertors is

  function Standard_Polynomial_to_Hexa_Double
             ( p : Standard_Complex_Polynomials.Poly )
             return Hexa_Double_Polynomials.Poly is
  
    res : Hexa_Double_Polynomials.Poly
        := Hexa_Double_Polynomials.Null_Poly;       

    use Standard_Complex_Polynomials;

    procedure Term_to_Hexa_Double_Term ( t : in Term; c : out boolean ) is

      rt : Hexa_Double_Polynomials.Term;
      cf : constant double_float := Standard_Complex_Numbers.REAL_PART(t.cf);
      dd : constant hexa_double := create(cf);

    begin
      rt.cf := dd;
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Hexa_Double_Polynomials.Add(res,rt);
      Hexa_Double_Polynomials.Clear(rt);
      c := true;
    end Term_to_Hexa_Double_Term;
    procedure P2DD is new Visiting_Iterator(Term_to_Hexa_Double_Term);

  begin
    P2DD(p);
    return res;
  end Standard_Polynomial_to_Hexa_Double;

  function Standard_Poly_Sys_to_Hexa_Double
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Hexa_Double_Poly_Systems.Poly_Sys is

    res : Hexa_Double_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Standard_Polynomial_to_Hexa_Double(p(i));
    end loop;
    return res;
  end Standard_Poly_Sys_to_Hexa_Double;

  function Multprec_Polynomial_to_Hexa_Double
             ( p : Multprec_Complex_Polynomials.Poly )
             return Hexa_Double_Polynomials.Poly is
  
    res : Hexa_Double_Polynomials.Poly
        := Hexa_Double_Polynomials.Null_Poly;       

    use Multprec_Complex_Polynomials;

    procedure Term_to_Hexa_Double_Term ( t : in Term; c : out boolean ) is

      rt : Hexa_Double_Polynomials.Term;
      cf : Floating_Number := Multprec_Complex_Numbers.REAL_PART(t.cf);
      dd : constant Hexa_Double := to_Hexa_Double(cf);

    begin
      rt.cf := dd;
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Hexa_Double_Polynomials.Add(res,rt);
      Hexa_Double_Polynomials.Clear(rt);
      Multprec_Floating_Numbers.Clear(cf);
      c := true;
    end Term_to_Hexa_Double_Term;
    procedure MP2DD is new Visiting_Iterator(Term_to_Hexa_Double_Term);

  begin
    MP2DD(p);
    return res;
  end Multprec_Polynomial_to_Hexa_Double;

  function Multprec_Poly_Sys_to_Hexa_Double
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return Hexa_Double_Poly_Systems.Poly_Sys is

    res : Hexa_Double_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Multprec_Polynomial_to_Hexa_Double(p(i));
    end loop;
    return res;
  end Multprec_Poly_Sys_to_Hexa_Double;

  function Standard_Polynomial_to_HexaDobl_Complex
             ( p : Standard_Complex_Polynomials.Poly )
             return HexaDobl_Complex_Polynomials.Poly is

    res : HexaDobl_Complex_Polynomials.Poly
        := HexaDobl_Complex_Polynomials.Null_Poly;

    use Standard_Complex_Polynomials;

    procedure Term_to_HexaDobl_Complex_Term ( t : in Term; c : out boolean ) is

      rt : HexaDobl_Complex_Polynomials.Term;

    begin
      rt.cf := Standard_to_HexaDobl_Complex(t.cf);
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      HexaDobl_Complex_Polynomials.Add(res,rt);
      HexaDobl_Complex_Polynomials.Clear(rt);
      c := true;
    end Term_to_HexaDobl_Complex_Term;
    procedure P2DC is new Visiting_Iterator(Term_to_HexaDobl_Complex_Term);

  begin
    P2DC(p);
    return res;
  end Standard_Polynomial_to_HexaDobl_Complex;

  function Standard_Poly_Sys_to_HexaDobl_Complex
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys is

    res : HexaDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Standard_Polynomial_to_HexaDobl_Complex(p(i));
    end loop;
    return res;
  end Standard_Poly_Sys_to_HexaDobl_Complex;

  function Multprec_Polynomial_to_HexaDobl_Complex
             ( p : Multprec_Complex_Polynomials.Poly )
             return HexaDobl_Complex_Polynomials.Poly is

    res : HexaDobl_Complex_Polynomials.Poly
        := HexaDobl_Complex_Polynomials.Null_Poly;

    use Multprec_Complex_Polynomials;

    procedure Term_to_HexaDobl_Complex_Term ( t : in Term; c : out boolean ) is

      rt : HexaDobl_Complex_Polynomials.Term;

    begin
      rt.cf := Multprec_to_HexaDobl_Complex(t.cf);
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      HexaDobl_Complex_Polynomials.Add(res,rt);
      HexaDobl_Complex_Polynomials.Clear(rt);
      c := true;
    end Term_to_HexaDobl_Complex_Term;
    procedure MP2DC is new Visiting_Iterator(Term_to_HexaDobl_Complex_Term);

  begin
    MP2DC(p);
    return res;
  end Multprec_Polynomial_to_HexaDobl_Complex;

  function Multprec_Poly_Sys_to_HexaDobl_Complex
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys is

    res : HexaDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Multprec_Polynomial_to_HexaDobl_Complex(p(i));
    end loop;
    return res;
  end Multprec_Poly_Sys_to_HexaDobl_Complex;

  function Standard_Laurential_to_HexaDobl_Complex
             ( p : Standard_Complex_Laurentials.Poly )
             return HexaDobl_Complex_Laurentials.Poly is

    res : HexaDobl_Complex_Laurentials.Poly
        := HexaDobl_Complex_Laurentials.Null_Poly;

    use Standard_Complex_Laurentials;

    procedure Term_to_HexaDobl_Complex_Term ( t : in Term; c : out boolean ) is

      rt : HexaDobl_Complex_Laurentials.Term;

    begin
      rt.cf := Standard_to_HexaDobl_Complex(t.cf);
      rt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
      HexaDobl_Complex_Laurentials.Add(res,rt);
      HexaDobl_Complex_Laurentials.Clear(rt);
      c := true;
    end Term_to_HexaDobl_Complex_Term;
    procedure L2DC is new Visiting_Iterator(Term_to_HexaDobl_Complex_Term);

  begin
    L2DC(p);
    return res;
  end Standard_Laurential_to_HexaDobl_Complex;

  function Standard_Laur_Sys_to_HexaDobl_Complex
             ( p : Standard_Complex_Laur_Systems.Laur_Sys )
             return HexaDobl_Complex_Laur_Systems.Laur_Sys is

    res : HexaDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Standard_Laurential_to_HexaDobl_Complex(p(i));
    end loop;
    return res;
  end Standard_Laur_Sys_to_HexaDobl_Complex;

  function Multprec_Laurential_to_HexaDobl_Complex
             ( p : Multprec_Complex_Laurentials.Poly )
             return HexaDobl_Complex_Laurentials.Poly is

    res : HexaDobl_Complex_Laurentials.Poly
        := HexaDobl_Complex_Laurentials.Null_Poly;

    use Multprec_Complex_Laurentials;

    procedure Term_to_HexaDobl_Complex_Term ( t : in Term; c : out boolean ) is

      rt : HexaDobl_Complex_Laurentials.Term;

    begin
      rt.cf := Multprec_to_HexaDobl_Complex(t.cf);
      rt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
      HexaDobl_Complex_Laurentials.Add(res,rt);
      HexaDobl_Complex_Laurentials.Clear(rt);
      c := true;
    end Term_to_HexaDobl_Complex_Term;
    procedure ML2DC is new Visiting_Iterator(Term_to_HexaDobl_Complex_Term);

  begin
    ML2DC(p);
    return res;
  end Multprec_Laurential_to_HexaDobl_Complex;

  function Multprec_Laur_Sys_to_HexaDobl_Complex
             ( p : Multprec_Complex_Laur_Systems.Laur_Sys )
             return HexaDobl_Complex_Laur_Systems.Laur_Sys is

    res : HexaDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Multprec_Laurential_to_HexaDobl_Complex(p(i));
    end loop;
    return res;
  end Multprec_Laur_Sys_to_HexaDobl_Complex;

  function Hexa_Double_to_Standard_Polynomial
             ( p : Hexa_Double_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    use Hexa_Double_Polynomials;

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
  end Hexa_Double_to_Standard_Polynomial;

  function Hexa_Double_to_Standard_Poly_Sys
             ( p : Hexa_Double_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Hexa_Double_to_Standard_Polynomial(p(i));
    end loop;
    return res;
  end Hexa_Double_to_Standard_Poly_Sys;

  function Hexa_Double_to_Multprec_Polynomial
             ( p : Hexa_Double_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly is

    res : Multprec_Complex_Polynomials.Poly
        := Multprec_Complex_Polynomials.Null_Poly;

    use Hexa_Double_Polynomials;

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
  end Hexa_Double_to_Multprec_Polynomial;

  function Hexa_Double_to_Multprec_Poly_Sys
             ( p : Hexa_Double_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys is

    res : Multprec_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Hexa_Double_to_Multprec_Polynomial(p(i));
    end loop;
    return res;
  end Hexa_Double_to_Multprec_Poly_Sys;

  function HexaDobl_Complex_to_Standard_Polynomial
             ( p : HexaDobl_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    use HexaDobl_Complex_Polynomials;

    procedure Term_to_Standard_Term ( t : in Term; c : out boolean ) is

      rt : Standard_Complex_Polynomials.Term;

    begin
      rt.cf := HexaDobl_Complex_to_Standard(t.cf);
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Standard_Complex_Polynomials.Add(res,rt);
      Standard_Complex_Polynomials.Clear(rt);
      c := true;
    end Term_to_Standard_Term;
    procedure DC2P is new Visiting_Iterator(Term_to_Standard_Term);

  begin
    DC2P(p);
    return res;
  end HexaDobl_Complex_to_Standard_Polynomial;

  function HexaDobl_Complex_to_Standard_Poly_Sys
             ( p : HexaDobl_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := HexaDobl_Complex_to_Standard_Polynomial(p(i));
    end loop;
    return res;
  end HexaDobl_Complex_to_Standard_Poly_Sys;

  function HexaDobl_Complex_to_Multprec_Polynomial
             ( p : HexaDobl_Complex_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly is

    res : Multprec_Complex_Polynomials.Poly
        := Multprec_Complex_Polynomials.Null_Poly;

    use HexaDobl_Complex_Polynomials;

    procedure Term_to_Multprec_Term ( t : in Term; c : out boolean ) is

      rt : Multprec_Complex_Polynomials.Term;

    begin
      rt.cf := HexaDobl_Complex_to_Multprec(t.cf);
      rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Multprec_Complex_Polynomials.Add(res,rt);
      Multprec_Complex_Polynomials.Clear(rt);
      c := true;
    end Term_to_Multprec_Term;
    procedure DC2MP is new Visiting_Iterator(Term_to_Multprec_Term);

  begin
    DC2MP(p);
    return res;
  end HexaDobl_Complex_to_Multprec_Polynomial;

  function HexaDobl_Complex_to_Multprec_Poly_Sys
             ( p : HexaDobl_Complex_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys is

    res : Multprec_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := HexaDobl_Complex_to_Multprec_Polynomial(p(i));
    end loop;
    return res;
  end HexaDobl_Complex_to_Multprec_Poly_Sys;

  function HexaDobl_Complex_to_Standard_Laurential
             ( p : HexaDobl_Complex_Laurentials.Poly )
             return Standard_Complex_Laurentials.Poly is

    res : Standard_Complex_Laurentials.Poly
        := Standard_Complex_Laurentials.Null_Poly;

    use HexaDobl_Complex_Laurentials;

    procedure Term_to_Standard_Term ( t : in Term; c : out boolean ) is

      rt : Standard_Complex_Laurentials.Term;

    begin
      rt.cf := HexaDobl_Complex_to_Standard(t.cf);
      rt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
      Standard_Complex_Laurentials.Add(res,rt);
      Standard_Complex_Laurentials.Clear(rt);
      c := true;
    end Term_to_Standard_Term;
    procedure DC2L is new Visiting_Iterator(Term_to_Standard_Term);

  begin
    DC2L(p);
    return res;
  end HexaDobl_Complex_to_Standard_Laurential;

  function HexaDobl_Complex_to_Standard_Laur_Sys
             ( p : HexaDobl_Complex_Laur_Systems.Laur_Sys )
             return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := HexaDobl_Complex_to_Standard_Laurential(p(i));
    end loop;
    return res;
  end HexaDobl_Complex_to_Standard_Laur_Sys;

  function HexaDobl_Complex_to_Multprec_Laurential
             ( p : HexaDobl_Complex_Laurentials.Poly )
             return Multprec_Complex_Laurentials.Poly is

    res : Multprec_Complex_Laurentials.Poly
        := Multprec_Complex_Laurentials.Null_Poly;

    use HexaDobl_Complex_Laurentials;

    procedure Term_to_Multprec_Term ( t : in Term; c : out boolean ) is

      rt : Multprec_Complex_Laurentials.Term;

    begin
      rt.cf := HexaDobl_Complex_to_Multprec(t.cf);
      rt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
      Multprec_Complex_Laurentials.Add(res,rt);
      Multprec_Complex_Laurentials.Clear(rt);
      c := true;
    end Term_to_Multprec_Term;
    procedure DC2ML is new Visiting_Iterator(Term_to_Multprec_Term);

  begin
    DC2ML(p);
    return res;
  end HexaDobl_Complex_to_Multprec_Laurential;

  function HexaDobl_Complex_to_Multprec_Laur_Sys
             ( p : HexaDobl_Complex_Laur_Systems.Laur_Sys )
             return Multprec_Complex_Laur_Systems.Laur_Sys is

    res : Multprec_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := HexaDobl_Complex_to_Multprec_Laurential(p(i));
    end loop;
    return res;
  end HexaDobl_Complex_to_Multprec_Laur_Sys;

end HexaDobl_Polynomial_Convertors;
