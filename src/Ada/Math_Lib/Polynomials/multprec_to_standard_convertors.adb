with Standard_Natural_Vectors;
with Multprec_Complex_Number_Tools;      use Multprec_Complex_Number_Tools;

package body Multprec_to_Standard_Convertors is

-- AUXILIARIES :

  function Convert ( d : Multprec_Complex_Polynomials.Degrees )
                   return Standard_Complex_Polynomials.Degrees is

    res : Standard_Complex_Polynomials.Degrees;

  begin
    res := new Standard_Natural_Vectors.Vector(d'range);
    for i in res'range loop
      res(i) := d(i);
    end loop;
    return res;
  end Convert;

  function Convert ( t : Multprec_Complex_Polynomials.Term ) 
                   return Standard_Complex_Polynomials.Term is

    res : Standard_Complex_Polynomials.Term;

  begin
    res.cf := Round(t.cf);
    res.dg := Convert(t.dg);
    return res;
  end Convert;

-- TARGET FUNCTIONS :

  function Convert ( p : Multprec_Complex_Polynomials.Poly )
                   return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    procedure Convert_Term ( t : in Multprec_Complex_Polynomials.Term;
                             continue : out boolean ) is

      ct : Standard_Complex_Polynomials.Term := Convert(t);

    begin
      Standard_Complex_Polynomials.Add(res,ct);
      Standard_Complex_Polynomials.Clear(ct);
      continue := true;
    end Convert_Term;
    procedure Convert_Terms is 
      new Multprec_Complex_Polynomials.Visiting_Iterator(Convert_Term);

  begin
    Convert_Terms(p);
    return res;
  end Convert;

  function Convert ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
                   return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Convert(p(i));
    end loop;
    return res;
  end Convert;

end Multprec_to_Standard_Convertors;
