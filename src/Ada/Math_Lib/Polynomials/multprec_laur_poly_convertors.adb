with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Standard_Natural_Vectors;

package body Multprec_Laur_Poly_Convertors is

  function Negative ( d : Multprec_Complex_Laurentials.Degrees ) 
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
              ( p : Multprec_Complex_Laurentials.Poly ) return boolean is

    use Multprec_Complex_Laurentials;

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
              ( p : Multprec_Complex_Laurentials.Poly )
              return Multprec_Complex_Polynomials.Poly is

    res : Multprec_Complex_Polynomials.Poly
        := Multprec_Complex_Polynomials.Null_Poly;

    procedure Convert ( t : in Multprec_Complex_Laurentials.Term;
                        c : out boolean ) is

      pt : Multprec_Complex_Polynomials.Term;

    begin
      Copy(t.cf,pt.cf);
      pt.dg := new Standard_Natural_Vectors.Vector(t.dg'range);
      for i in pt.dg'range loop
        pt.dg(i) := natural32(t.dg(i));
      end loop;
      Multprec_Complex_Polynomials.Add(res,pt);
      c := true;
    end Convert;
    procedure Convert_Terms is 
      new Multprec_Complex_Laurentials.Visiting_Iterator(Convert);

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
             ( p : Multprec_Complex_Laurentials.Poly )
             return Multprec_Complex_Polynomials.Poly is

    res : Multprec_Complex_Polynomials.Poly;
    tt : Multprec_Complex_Laurentials.Term;

  begin
    Laurent_Polynomial_to_Polynomial(p,tt,res);
    Multprec_Complex_Laurentials.Clear(tt);
    return res;
  end Laurent_Polynomial_to_Polynomial;

  procedure Laurent_Polynomial_to_Polynomial
             ( l : in Multprec_Complex_Laurentials.Poly;
               t : out Multprec_Complex_Laurentials.Term; 
               p : out Multprec_Complex_Polynomials.Poly ) is

    min : constant Multprec_Complex_Laurentials.Degrees 
        := Multprec_Complex_Laurentials.Minimal_Degrees(l);
    tt : Multprec_Complex_Laurentials.Term;

  begin
    for i in min'range loop
      if min(i) < 0
       then min(i) := -min(i);   -- only multiply if negative!
       else min(i) := 0;
      end if;
    end loop;
    tt.cf := Create(integer32(1));
    tt.dg := min;
    p := Laurent_Polynomial_to_Polynomial(l,tt); t := tt;
  end Laurent_Polynomial_to_Polynomial;

  function Laurent_Polynomial_to_Polynomial
            ( l : Multprec_Complex_Laurentials.Poly;
              t : Multprec_Complex_Laurentials.Term )
            return Multprec_Complex_Polynomials.Poly is

    res : Multprec_Complex_Polynomials.Poly;
    use Multprec_Complex_Laurentials;

    procedure Laurent_Term_to_Term ( tt : in Term; cont : out boolean ) is

      rt : Multprec_Complex_Polynomials.Term;

    begin
      Copy(tt.cf,rt.cf);
      rt.dg := new Standard_Natural_Vectors.Vector(tt.dg'range);
      for i in tt.dg'range loop
        rt.dg(i) := natural32(tt.dg(i)) + natural32(t.dg(i));
      end loop;
      Multprec_Complex_Polynomials.Add(res,rt);
      Multprec_Complex_Polynomials.Clear(rt);
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

end Multprec_Laur_Poly_Convertors;
