with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors; 

package body QuadDobl_Poly_Laur_Convertors is

  function Polynomial_to_Laurent_Polynomial
             ( p : QuadDobl_Complex_Polynomials.Poly )
             return QuadDobl_Complex_Laurentials.Poly is

    res : QuadDobl_Complex_Laurentials.Poly
        := QuadDobl_Complex_Laurentials.Null_Poly;

    use QuadDobl_Complex_Polynomials;

    procedure Term_to_Laurent_Term ( t : in Term; cont : out boolean ) is

      rt : QuadDobl_Complex_Laurentials.Term;

    begin
      rt.cf := t.cf;
      rt.dg := new Standard_Integer_Vectors.Vector(t.dg'range);
      for i in t.dg'range loop
        rt.dg(i) := integer32(t.dg(i));
      end loop;
      QuadDobl_Complex_Laurentials.Add(res,rt);
      QuadDobl_Complex_Laurentials.Clear(rt);
      cont := true;
    end Term_to_Laurent_Term;
    procedure P2LP is new Visiting_Iterator(Term_to_Laurent_Term);

  begin
    P2LP(p);
    return res;
  end Polynomial_to_Laurent_Polynomial;

  function Polynomial_to_Laurent_System ( p : Poly_Sys ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Polynomial_to_Laurent_Polynomial(p(i));
    end loop;
    return res;
  end Polynomial_to_Laurent_System;

end QuadDobl_Poly_Laur_Convertors;
