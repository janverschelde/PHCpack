with Double_Double_Numbers;          use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with Standard_Natural_Vectors;

package body DoblDobl_Complex_to_Real_Poly is

  function Is_Real ( p : DoblDobl_Complex_Polynomials.Poly ) return boolean is

    res : boolean := true;
    zero : constant double_double := create(0.0);

    procedure Is_Real_Term ( t : in DoblDobl_Complex_Polynomials.Term;
                             continue : out boolean ) is
    begin
      if Imag_Part(t.cf) = zero
       then continue := true;
       else continue := false; res := false;
      end if;
    end Is_Real_Term;
    procedure Real_Terms is
      new DoblDobl_Complex_Polynomials.Visiting_Iterator(Is_Real_Term);

  begin
    Real_Terms(p);
    return res;
  end Is_Real;

  function Is_Real
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys ) return boolean is
  begin
    for i in p'range loop
      if not Is_Real(p(i))
       then return false;
      end if;
    end loop;
    return true;
  end Is_Real;

  function Convert_Complex_to_Real
             ( p : DoblDobl_Complex_Polynomials.Poly )
             return Double_Double_Polynomials.Poly is

    res : Double_Double_Polynomials.Poly
        := Double_Double_Polynomials.Null_Poly;

    procedure Convert ( t : in DoblDobl_Complex_Polynomials.Term;
                        continue : out boolean ) is

      ct : Double_Double_Polynomials.Term;

    begin
      ct.cf := Real_Part(t.cf);
      ct.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      Double_Double_Polynomials.Add(res,ct);
      Double_Double_Polynomials.Clear(ct);
      continue := true;
    end Convert;
    procedure Convert is
      new DoblDobl_Complex_Polynomials.Visiting_Iterator(Convert);

  begin
    Convert(p);
    return res;
  end Convert_Complex_to_Real;

  function Convert_Complex_to_Real
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return Double_Double_Poly_Systems.Poly_Sys is

    res : Double_Double_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Convert_Complex_to_Real(p(i));
    end loop;
    return res;
  end Convert_Complex_to_Real;

  function Convert_Complex_to_Real
             ( p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys )
             return Double_Double_Poly_Systems.Link_to_Poly_Sys is

    res : Double_Double_Poly_Systems.Link_to_Poly_Sys;

    use DoblDobl_Complex_Poly_Systems;

  begin
    if p /= null then
      declare
        s : Double_Double_Poly_Systems.Poly_Sys(p'range);
      begin
        s := DoblDobl_Complex_to_Real_Poly.Convert_Complex_to_Real(p.all);
        res := new Double_Double_Poly_Systems.Poly_Sys'(s);
      end;
    end if;
    return res;
  end Convert_Complex_to_Real;

  function Convert_Real_to_Complex
             ( p : Double_Double_Polynomials.Poly )
             return DoblDobl_Complex_Polynomials.Poly is

    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Null_Poly;

    procedure Convert ( t : in Double_Double_Polynomials.Term;
                        continue : out boolean ) is

      ct : DoblDobl_Complex_Polynomials.Term;

    begin
      ct.cf := Create(t.cf);
      ct.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
      DoblDobl_Complex_Polynomials.Add(res,ct);
      DoblDobl_Complex_Polynomials.Clear(ct);
      continue := true;
    end Convert;
    procedure Convert is
      new Double_Double_Polynomials.Visiting_Iterator(Convert);

  begin
    Convert(p);
    return res;
  end Convert_Real_to_Complex;

  function Convert_Real_to_Complex
             ( p : Double_Double_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Convert_Real_to_Complex(p(i));
    end loop;
    return res;
  end Convert_Real_to_Complex;

  function Convert_Real_to_Complex
             ( p : Double_Double_Poly_Systems.Link_to_Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

    use Double_Double_Poly_Systems;

  begin
    if p /= null then
      declare
        s : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
      begin
        s := DoblDobl_Complex_to_Real_Poly.Convert_Real_to_Complex(p.all);
        res := new DoblDobl_Complex_Poly_Systems.Poly_Sys'(s);
      end;
    end if;
    return res;
  end Convert_Real_to_Complex;

end DoblDobl_Complex_to_Real_Poly;
