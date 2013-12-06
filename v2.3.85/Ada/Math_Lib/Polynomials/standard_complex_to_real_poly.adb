with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;

package body Standard_Complex_to_Real_Poly is

  function Is_Real ( p : Standard_Complex_Polynomials.Poly ) return boolean is

    res : boolean := true;

    procedure Is_Real_Term ( t : in Standard_Complex_Polynomials.Term;
                             continue : out boolean ) is
    begin
      if Imag_Part(t.cf) = 0.0
       then continue := true;
       else continue := false;
            res := false;
      end if;
    end Is_Real_Term;
    procedure Real_Terms is
      new Standard_Complex_Polynomials.Visiting_Iterator(Is_Real_Term);

  begin
    Real_Terms(p);
    return res;
  end Is_Real;

  function Is_Real
             ( p : Standard_Complex_Poly_Systems.Poly_Sys ) return boolean is
  begin
    for i in p'range loop
      if not Is_Real(p(i))
       then return false;
      end if;
    end loop;
    return true;
  end Is_Real;

  function Convert_Complex_to_Real
             ( p : Standard_Complex_Polynomials.Poly )
             return Standard_Floating_Polynomials.Poly is

    res : Standard_Floating_Polynomials.Poly
        := Standard_Floating_Polynomials.Null_Poly;

    procedure Convert ( t : in Standard_Complex_Polynomials.Term;
                        continue : out boolean ) is

      ct : Standard_Floating_Polynomials.Term;

    begin
      ct.cf := Real_Part(t.cf);
      ct.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Standard_Floating_Polynomials.Add(res,ct);
      Standard_Floating_Polynomials.Clear(ct);
      continue := true;
    end Convert;
    procedure Convert is
      new Standard_Complex_Polynomials.Visiting_Iterator(Convert);

  begin
    Convert(p);
    return res;
  end Convert_Complex_to_Real;

  function Convert_Complex_to_Real
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Floating_Poly_Systems.Poly_Sys is

    res : Standard_Floating_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Convert_Complex_to_Real(p(i));
    end loop;
    return res;
  end Convert_Complex_to_Real;

  function Convert_Complex_to_Real
             ( p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys )
             return Standard_Floating_Poly_Systems.Link_to_Poly_Sys is

    res : Standard_Floating_Poly_Systems.Link_to_Poly_Sys;

    use Standard_Complex_Poly_Systems;

  begin
    if p /= null then
      declare
        s : Standard_Floating_Poly_Systems.Poly_Sys(p'range);
      begin
        s := Standard_Complex_to_Real_Poly.Convert_Complex_to_Real(p.all);
        res := new Standard_Floating_Poly_Systems.Poly_Sys'(s);
      end;
    end if;
    return res;
  end Convert_Complex_to_Real;

  function Convert_Real_to_Complex
             ( p : Standard_Floating_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    procedure Convert ( t : in Standard_Floating_Polynomials.Term;
                        continue : out boolean ) is

      ct : Standard_Complex_Polynomials.Term;

    begin
      ct.cf := Create(t.cf);
      ct.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Standard_Complex_Polynomials.Add(res,ct);
      Standard_Complex_Polynomials.Clear(ct);
      continue := true;
    end Convert;
    procedure Convert is
      new Standard_Floating_Polynomials.Visiting_Iterator(Convert);

  begin
    Convert(p);
    return res;
  end Convert_Real_to_Complex;

  function Convert_Real_to_Complex
             ( p : Standard_Floating_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Convert_Real_to_Complex(p(i));
    end loop;
    return res;
  end Convert_Real_to_Complex;

  function Convert_Real_to_Complex
             ( p : Standard_Floating_Poly_Systems.Link_to_Poly_Sys )
             return Standard_Complex_Poly_Systems.Link_to_Poly_Sys is

    res : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

    use Standard_Floating_Poly_Systems;

  begin
    if p /= null then
      declare
        s : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
      begin
        s := Standard_Complex_to_Real_Poly.Convert_Real_to_Complex(p.all);
        res := new Standard_Complex_Poly_Systems.Poly_Sys'(s);
      end;
    end if;
    return res;
  end Convert_Real_to_Complex;

end Standard_Complex_to_Real_Poly;
