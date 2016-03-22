with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;

package body Standard_Laur_Poly_Convertors is

  function Negative ( d : Standard_Complex_Laurentials.Degrees ) 
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
              ( p : Standard_Complex_Laurentials.Poly ) return boolean is

    use Standard_Complex_Laurentials;

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
              ( p : Standard_Complex_Laurentials.Poly )
              return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    procedure Convert ( t : in Standard_Complex_Laurentials.Term;
                        c : out boolean ) is

      pt : Standard_Complex_Polynomials.Term;

    begin
      pt.cf := t.cf;
      pt.dg := new Standard_Natural_Vectors.Vector(t.dg'range);
      for i in pt.dg'range loop
        pt.dg(i) := natural32(t.dg(i));
      end loop;
      Standard_Complex_Polynomials.Add(res,pt);
      c := true;
    end Convert;
    procedure Convert_Terms is 
      new Standard_Complex_Laurentials.Visiting_Iterator(Convert);

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
             ( p : Standard_Complex_Laurentials.Poly )
             return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly;
    tt : Standard_Complex_Laurentials.Term;

  begin
    Laurent_Polynomial_to_Polynomial(p,tt,res);
    Standard_Complex_Laurentials.Clear(tt);
    return res;
  end Laurent_Polynomial_to_Polynomial;

  procedure Laurent_Polynomial_to_Polynomial
             ( L : in Standard_Complex_Laurentials.Poly;
               t : out Standard_Complex_Laurentials.Term; 
               p : out Standard_Complex_Polynomials.Poly ) is

    min : Standard_Complex_Laurentials.Degrees 
        := Standard_Complex_Laurentials.Minimal_Degrees(L);
    tt : Standard_Complex_Laurentials.Term;

  begin
    for i in min'range loop
      if min(i) < 0
       then min(i) := -min(i);  -- only multiply if negative!
       else min(i) := 0;
      end if;
    end loop;
    tt.cf := Create(1.0);
    tt.dg := min;
    p := Laurent_Polynomial_to_Polynomial(l,tt); t := tt;
  end Laurent_Polynomial_to_Polynomial;

  function Laurent_Polynomial_to_Polynomial
            ( L : Standard_Complex_Laurentials.Poly;
              t : Standard_Complex_Laurentials.Term )
            return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly;
    use Standard_Complex_Laurentials;

    procedure Laurent_Term_to_Term ( tt : in Term; cont : out boolean ) is

      rt : Standard_Complex_Polynomials.Term;

    begin
      rt.cf := tt.cf;
      rt.dg := new Standard_Natural_Vectors.Vector(tt.dg'range);
      for i in tt.dg'range loop
        rt.dg(i) := natural32(tt.dg(i) + t.dg(i));
      end loop;
      Standard_Complex_Polynomials.Add(res,rt);
      Standard_Complex_Polynomials.Clear(rt);
      cont := true;
    end Laurent_Term_to_Term;
    procedure LP2P is new Visiting_Iterator(Laurent_Term_to_Term);

  begin
    LP2P(l);
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

end Standard_Laur_Poly_Convertors;
