with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Standard_Natural_Vectors;

package body Multprec_Complex_to_Real_Poly is

  function Convert_Complex_to_Real
             ( p : Multprec_Complex_Polynomials.Poly )
             return Multprec_Floating_Polynomials.Poly is

    res : Multprec_Floating_Polynomials.Poly
        := Multprec_Floating_Polynomials.Null_Poly;

    procedure Convert ( t : in Multprec_Complex_Polynomials.Term;
                        continue : out boolean ) is

      ct : Multprec_Floating_Polynomials.Term;

    begin
      ct.cf := Real_Part(t.cf);
      ct.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Multprec_Floating_Polynomials.Add(res,ct);
      Multprec_Floating_Polynomials.Clear(ct);
      continue := true;
    end Convert;
    procedure Convert is
      new Multprec_Complex_Polynomials.Visiting_Iterator(Convert);

  begin
    Convert(p);
    return res;
  end Convert_Complex_to_Real;

  function Convert_Complex_to_Real
             ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
             return Multprec_Floating_Poly_Systems.Poly_Sys is

    res : Multprec_Floating_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Convert_Complex_to_Real(p(i));
    end loop;
    return res;
  end Convert_Complex_to_Real;

  function Convert_Real_to_Complex
             ( p : Multprec_Floating_Polynomials.Poly )
             return Multprec_Complex_Polynomials.Poly is

    res : Multprec_Complex_Polynomials.Poly
        := Multprec_Complex_Polynomials.Null_Poly;

    procedure Convert ( t : in Multprec_Floating_Polynomials.Term;
                        continue : out boolean ) is

      ct : Multprec_Complex_Polynomials.Term;

    begin
      ct.cf := Create(t.cf);
      ct.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Multprec_Complex_Polynomials.Add(res,ct);
      Multprec_Complex_Polynomials.Clear(ct);
      continue := true;
    end Convert;
    procedure Convert is
      new Multprec_Floating_Polynomials.Visiting_Iterator(Convert);

  begin
    Convert(p);
    return res;
  end Convert_Real_to_Complex;

  function Convert_Real_to_Complex
             ( p : Multprec_Floating_Poly_Systems.Poly_Sys )
             return Multprec_Complex_Poly_Systems.Poly_Sys is

    res : Multprec_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Convert_Real_to_Complex(p(i));
    end loop;
    return res;
  end Convert_Real_to_Complex;

end Multprec_Complex_to_Real_Poly;
