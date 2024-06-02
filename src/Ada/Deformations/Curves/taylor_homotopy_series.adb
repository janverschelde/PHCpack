with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Series;

package body Taylor_Homotopy_Series is

  function Make ( m : Double_Taylor_Homotopies.Taylor_Monomial )
                return Standard_CSeries_Polynomials.Term is

    res : Standard_CSeries_Polynomials.Term;
    srs : Standard_Complex_Series.Series(m.deg);

  begin
    if m.pwr = 0.0 then
      srs.cff(0) := m.cst;
      for i in 1..m.deg loop
        srs.cff(i) := Create(0.0);
      end loop;
    else
      srs.cff(0) := m.cff(0)*m.cst;
      for i in 1..m.deg loop
        srs.cff(i) := m.cff(i)*m.cst;
      end loop;
    end if;
    res.cf := new Standard_Complex_Series.Series'(srs);
    res.dg := new Standard_Natural_Vectors.Vector(m.exp'range);
    for i in m.exp'range loop
      res.dg(i) := natural32(m.exp(i));
    end loop;
    return res;
  end Make;

  function Make ( m : Double_Taylor_Homotopies.Link_to_Taylor_Monomial )
                return Standard_CSeries_Polynomials.Term is

    res : Standard_CSeries_Polynomials.Term;

    use Double_Taylor_Homotopies;

  begin
    if m /= null
     then res := Make(m.all);
    end if;
    return res;
  end Make;

  function Make ( v : Double_Taylor_Homotopies.Taylor_Monomial_Vector )
                return Standard_CSeries_Polynomials.Poly is

    res : Standard_CSeries_Polynomials.Poly;

  begin
    for i in v'range loop
      declare
        t : Standard_CSeries_Polynomials.Term := Make(v(i));
      begin
        Standard_CSeries_Polynomials.Add(res,t);
        Standard_CSeries_Polynomials.Clear(t);
      end;
    end loop;
    return res;
  end Make;

  function Make ( v : Double_Taylor_Homotopies.Link_to_Taylor_Monomial_Vector )
                return Standard_CSeries_Polynomials.Poly is

    res : Standard_CSeries_Polynomials.Poly;

    use Double_Taylor_Homotopies;

  begin
    if v /= null
     then res := Make(v.all);
    end if;
    return res;
  end Make;

  function Make ( h : Double_Taylor_Homotopies.Taylor_Homotopy )
                return Standard_CSeries_Poly_Systems.Poly_Sys is

    res : Standard_CSeries_Poly_Systems.Poly_Sys(h'range);

  begin
    for i in h'range loop
      res(i) := Make(h(i));
    end loop;
    return res;
  end Make;

end Taylor_Homotopy_Series;
