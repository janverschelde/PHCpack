with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Complex_Numbers;
with Standard_Natural_Vectors;

package body Coefficient_Supported_Polynomials is

  function Create_Standard_Polynomial
             ( e : Standard_Natural_VecVecs.VecVec )
             return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly;
    t : Standard_Complex_Polynomials.Term;

  begin
    for i in e'range loop
      t.cf := Standard_Complex_Numbers.Create(1.0);
      t.dg := new Standard_Natural_Vectors.Vector'(e(i).all);
      Standard_Complex_Polynomials.Add(res,t);
      Standard_Complex_Polynomials.Clear(t);
    end loop;
    return res;
  end Create_Standard_Polynomial;

  function Create_Standard_Polynomial
             ( c : Standard_Complex_Vectors.Vector;
               e : Standard_Natural_VecVecs.VecVec )
             return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly;
    t : Standard_Complex_Polynomials.Term;

  begin
    for i in e'range loop
      t.cf := c(i);
      t.dg := new Standard_Natural_Vectors.Vector'(e(i).all);
      Standard_Complex_Polynomials.Add(res,t);
      Standard_Complex_Polynomials.Clear(t);
    end loop;
    return res;
  end Create_Standard_Polynomial;

  procedure Coefficients_and_Supports
              ( p : in Standard_Complex_Polynomials.Poly;
                c : out Standard_Complex_Vectors.Vector;
                e : out Standard_Natural_VecVecs.VecVec ) is

    ind : integer32 := 0;
    use Standard_Complex_Polynomials;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      c(ind) := t.cf;
      e(ind) := new Standard_Natural_Vectors.Vector(t.dg'range);
      for i in t.dg'range loop
        e(ind)(i) := t.dg(i);
      end loop;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
  end Coefficients_and_Supports;

  procedure Coefficients_and_Supports
              ( p : in DoblDobl_Complex_Polynomials.Poly;
                c : out DoblDobl_Complex_Vectors.Vector;
                e : out Standard_Natural_VecVecs.VecVec ) is

    ind : integer32 := 0;
    use DoblDobl_Complex_Polynomials;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      c(ind) := t.cf;
      e(ind) := new Standard_Natural_Vectors.Vector(t.dg'range);
      for i in t.dg'range loop
        e(ind)(i) := t.dg(i);
      end loop;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
  end Coefficients_and_Supports;

  procedure Coefficients_and_Supports
              ( p : in QuadDobl_Complex_Polynomials.Poly;
                c : out QuadDobl_Complex_Vectors.Vector;
                e : out Standard_Natural_VecVecs.VecVec ) is

    ind : integer32 := 0;
    use QuadDobl_Complex_Polynomials;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      c(ind) := t.cf;
      e(ind) := new Standard_Natural_Vectors.Vector(t.dg'range);
      for i in t.dg'range loop
        e(ind)(i) := t.dg(i);
      end loop;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
  end Coefficients_and_Supports;

  procedure Coefficients_and_Supports
              ( p : in Multprec_Complex_Polynomials.Poly;
                c : out Multprec_Complex_Vectors.Vector;
                e : out Standard_Natural_VecVecs.VecVec ) is

    ind : integer32 := 0;
    use Multprec_Complex_Polynomials;

    procedure Visit_Term ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      Multprec_Complex_Numbers.Copy(t.cf,c(ind));
      e(ind) := new Standard_Natural_Vectors.Vector(t.dg'range);
      for i in t.dg'range loop
        e(ind)(i) := t.dg(i);
      end loop;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
  end Coefficients_and_Supports;

  procedure Coefficients_and_Supports
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                c : out Standard_Complex_VecVecs.VecVec;
                e : out Standard_Natural_VecVecs.Array_of_VecVecs ) is
  begin
    for i in p'range loop
      declare
        cnt : constant integer32
            := integer32(Standard_Complex_Polynomials.Number_of_Terms(p(i)));
        cff : Standard_Complex_Vectors.Vector(1..cnt);
        exp : Standard_Natural_VecVecs.VecVec(1..cnt);
      begin
        Coefficients_and_Supports(p(i),cff,exp);
        c(i) := new Standard_Complex_Vectors.Vector'(cff);
        e(i) := new Standard_Natural_VecVecs.VecVec'(exp);
      end;
    end loop;
  end Coefficients_and_Supports;

  function Create_DoblDobl_Polynomial
             ( e : Standard_Natural_VecVecs.VecVec )
             return DoblDobl_Complex_Polynomials.Poly is

    res : DoblDobl_Complex_Polynomials.Poly;
    t : DoblDobl_Complex_Polynomials.Term;

  begin
    for i in e'range loop
      t.cf := DoblDobl_Complex_Numbers.Create(integer(1));
      t.dg := new Standard_Natural_Vectors.Vector'(e(i).all);
      DoblDobl_Complex_Polynomials.Add(res,t);
      DoblDobl_Complex_Polynomials.Clear(t);
    end loop;
    return res;
  end Create_DoblDobl_Polynomial;

  function Create_DoblDobl_Polynomial
             ( c : DoblDobl_Complex_Vectors.Vector;
               e : Standard_Natural_VecVecs.VecVec )
             return DoblDobl_Complex_Polynomials.Poly is

    res : DoblDobl_Complex_Polynomials.Poly;
    t : DoblDobl_Complex_Polynomials.Term;

  begin
    for i in e'range loop
      t.cf := c(i);
      t.dg := new Standard_Natural_Vectors.Vector'(e(i).all);
      DoblDobl_Complex_Polynomials.Add(res,t);
      DoblDobl_Complex_Polynomials.Clear(t);
    end loop;
    return res;
  end Create_DoblDobl_Polynomial;

  function Create_QuadDobl_Polynomial
             ( e : Standard_Natural_VecVecs.VecVec )
             return QuadDobl_Complex_Polynomials.Poly is

    res : QuadDobl_Complex_Polynomials.Poly;
    t : QuadDobl_Complex_Polynomials.Term;

  begin
    for i in e'range loop
      t.cf := QuadDobl_Complex_Numbers.Create(integer(1));
      t.dg := new Standard_Natural_Vectors.Vector'(e(i).all);
      QuadDobl_Complex_Polynomials.Add(res,t);
      QuadDobl_Complex_Polynomials.Clear(t);
    end loop;
    return res;
  end Create_QuadDobl_Polynomial;

  function Create_QuadDobl_Polynomial
             ( c : QuadDobl_Complex_Vectors.Vector;
               e : Standard_Natural_VecVecs.VecVec )
             return QuadDobl_Complex_Polynomials.Poly is

    res : QuadDobl_Complex_Polynomials.Poly;
    t : QuadDobl_Complex_Polynomials.Term;

  begin
    for i in e'range loop
      t.cf := c(i);
      t.dg := new Standard_Natural_Vectors.Vector'(e(i).all);
      QuadDobl_Complex_Polynomials.Add(res,t);
      QuadDobl_Complex_Polynomials.Clear(t);
    end loop;
    return res;
  end Create_QuadDobl_Polynomial;

  function Create_Multprec_Polynomial
             ( e : Standard_Natural_VecVecs.VecVec )
             return Multprec_Complex_Polynomials.Poly is

    res : Multprec_Complex_Polynomials.Poly;
    t : Multprec_Complex_Polynomials.Term;

  begin
    for i in e'range loop
      t.cf := Multprec_Complex_Numbers.Create(integer(1));
      t.dg := new Standard_Natural_Vectors.Vector'(e(i).all);
      Multprec_Complex_Polynomials.Add(res,t);
      Multprec_Complex_Polynomials.Clear(t);
    end loop;
    return res;
  end Create_Multprec_Polynomial;

  function Create_Multprec_Polynomial
             ( c : Multprec_Complex_Vectors.Vector;
               e : Standard_Natural_VecVecs.VecVec )
             return Multprec_Complex_Polynomials.Poly is

    res : Multprec_Complex_Polynomials.Poly;
    t : Multprec_Complex_Polynomials.Term;

  begin
    for i in e'range loop
      Multprec_Complex_Numbers.Copy(c(i),t.cf);
      t.dg := new Standard_Natural_Vectors.Vector'(e(i).all);
      Multprec_Complex_Polynomials.Add(res,t);
      Multprec_Complex_Polynomials.Clear(t);
    end loop;
    return res;
  end Create_Multprec_Polynomial;

end Coefficient_Supported_Polynomials;
