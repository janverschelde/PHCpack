with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;

package body QuadDobl_Laur_Poly_Convertors is

  function Negative ( d : QuadDobl_Complex_Laurentials.Degrees ) 
                    return boolean is

  -- DESCRIPTION :
  --   Returns true if there is at least one negative element in d.

  begin
    for i in d'range loop
      if d(i) < 0
       then return true;
      end if;
    end loop;
    return false;
  end Negative;

  function Is_Genuine_Laurent
              ( p : QuadDobl_Complex_Laurentials.Poly ) return boolean is

    use QuadDobl_Complex_Laurentials;

    d : Degrees := Minimal_Degrees(p);
    res : constant boolean := Negative(d);

  begin
    Clear(d);
    return res;
  end Is_Genuine_Laurent;

  function Is_Genuine_Laurent ( p : Laur_Sys ) return boolean is
  begin
    for i in p'range loop
      if Is_Genuine_Laurent(p(i))
       then return true;
      end if;
    end loop;
    return false;
  end Is_Genuine_Laurent;

  function Positive_Laurent_Polynomial
              ( p : QuadDobl_Complex_Laurentials.Poly )
              return QuadDobl_Complex_Polynomials.Poly is

    res : QuadDobl_Complex_Polynomials.Poly
        := QuadDobl_Complex_Polynomials.Null_Poly;

    procedure Convert ( t : in QuadDobl_Complex_Laurentials.Term;
                        c : out boolean ) is

      pt : QuadDobl_Complex_Polynomials.Term;

    begin
      pt.cf := t.cf;
      pt.dg := new Standard_Natural_Vectors.Vector(t.dg'range);
      for i in pt.dg'range loop
        pt.dg(i) := natural32(t.dg(i));
      end loop;
      QuadDobl_Complex_Polynomials.Add(res,pt);
      c := true;
    end Convert;
    procedure Convert_Terms is 
      new QuadDobl_Complex_Laurentials.Visiting_Iterator(Convert);

  begin
    Convert_Terms(p);
    return res;
  end Positive_Laurent_Polynomial;

  function Positive_Laurent_Polynomial_System
             ( p : Laur_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Positive_Laurent_Polynomial(p(i));
    end loop;
    return res;
  end Positive_Laurent_Polynomial_System;

  function Laurent_Polynomial_to_Polynomial
             ( p : QuadDobl_Complex_Laurentials.Poly )
             return QuadDobl_Complex_Polynomials.Poly is

    res : QuadDobl_Complex_Polynomials.Poly;
    tt : QuadDobl_Complex_Laurentials.Term;

  begin
    Laurent_Polynomial_to_Polynomial(p,tt,res);
    QuadDobl_Complex_Laurentials.Clear(tt);
    return res;
  end Laurent_Polynomial_to_Polynomial;

  procedure Laurent_Polynomial_to_Polynomial
             ( L : in QuadDobl_Complex_Laurentials.Poly;
               t : out QuadDobl_Complex_Laurentials.Term; 
               p : out QuadDobl_Complex_Polynomials.Poly ) is

    min : QuadDobl_Complex_Laurentials.Degrees 
        := QuadDobl_Complex_Laurentials.Minimal_Degrees(L);
    tt : QuadDobl_Complex_Laurentials.Term;
    one : constant quad_double := create(1.0);

  begin
    for i in min'range loop
      if min(i) < 0
       then min(i) := -min(i);  -- only multiply if negative!
       else min(i) := 0;
      end if;
    end loop;
    tt.cf := QuadDobl_Complex_Numbers.Create(one);
    tt.dg := min;
    p := Laurent_Polynomial_to_Polynomial(L,tt); t := tt;
  end Laurent_Polynomial_to_Polynomial;

  function Laurent_Polynomial_to_Polynomial
            ( L : QuadDobl_Complex_Laurentials.Poly;
              t : QuadDobl_Complex_Laurentials.Term )
            return QuadDobl_Complex_Polynomials.Poly is

    res : QuadDobl_Complex_Polynomials.Poly;
    use QuadDobl_Complex_Laurentials;

    procedure Laurent_Term_to_Term ( tt : in Term; cont : out boolean ) is

      rt : QuadDobl_Complex_Polynomials.Term;

    begin
      rt.cf := tt.cf;
      rt.dg := new Standard_Natural_Vectors.Vector(tt.dg'range);
      for i in tt.dg'range loop
        rt.dg(i) := natural32(tt.dg(i) + t.dg(i));
      end loop;
      QuadDobl_Complex_Polynomials.Add(res,rt);
      QuadDobl_Complex_Polynomials.Clear(rt);
      cont := true;
    end Laurent_Term_to_Term;
    procedure LP2P is new Visiting_Iterator(Laurent_Term_to_Term);

  begin
    LP2P(L);
    return res;
  end Laurent_Polynomial_to_Polynomial;

  function Laurent_to_Polynomial_System ( p : Laur_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Laurent_Polynomial_to_Polynomial(p(i));
    end loop;
    return res;
  end Laurent_to_Polynomial_System;

end QuadDobl_Laur_Poly_Convertors;
