with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Series;

package body Standard_Polynomial_CSeries is

-- CONSTRUCTORS :

  function Create ( p : Standard_CSeries_Polynomials.Poly )
                  return Standard_Polynomial_CSeries.Poly is

    maxdeg : constant integer32 := 16;
    res : Standard_Polynomial_CSeries.Poly(maxdeg);

    procedure Add_Term ( t : in Standard_CSeries_Polynomials.Term;
                         cont : out boolean ) is

      cf : constant Standard_Complex_Series.Link_to_Series := t.cf;
   
    begin
      for k in 0..cf.deg loop
        declare
          rt : Standard_Complex_Polynomials.Term;
        begin
          rt.cf := cf.cff(k);
          rt.dg := Standard_Complex_Polynomials.Degrees(t.dg);
          Standard_Complex_Polynomials.Add(res.cff(k),rt);
        end;
      end loop;
      cont := true;
    end Add_Term;
    procedure Add_Terms is
      new Standard_CSeries_Polynomials.Visiting_Iterator(Add_Term);

  begin
    res.cff := (res.cff'range => Standard_Complex_Polynomials.Null_Poly);
    Add_Terms(p);
    return res;
  end Create;

  function Convert ( p : Standard_Complex_Polynomials.Poly;
                     d : integer32 )
                   return Standard_CSeries_Polynomials.Poly is

  -- DESCRIPTION :
  --   Each coefficient c of the polynomial p is stored as c*t^d
  --   as series coefficient in the polynomial on return.

    res : Standard_CSeries_Polynomials.Poly
        := Standard_CSeries_Polynomials.Null_Poly;

    procedure Add_Term ( t : in Standard_Complex_Polynomials.Term;
                         cont : out boolean ) is

      rt : Standard_CSeries_Polynomials.Term;

    begin
      rt.cf := Standard_Complex_Series.Create(0,d);
      for k in 0..d-1 loop
        rt.cf.cff(k) := Create(0.0);
      end loop;
      rt.cf.cff(d) := t.cf;
      rt.dg := Standard_CSeries_Polynomials.Degrees(t.dg);
      Standard_CSeries_Polynomials.Add(res,rt);
      cont := true;
    end Add_Term;
    procedure Add_Terms is new
      Standard_Complex_Polynomials.Visiting_Iterator(Add_Term);

  begin
    Add_Terms(p);
    return res;
  end Convert;

  function Create ( p : Standard_Polynomial_CSeries.Poly )
                  return Standard_CSeries_Polynomials.Poly is

    res : Standard_CSeries_Polynomials.Poly := Convert(p.cff(0),0);

  begin
    for k in 1..p.deg loop
      declare
        tpk : Standard_CSeries_Polynomials.Poly := Convert(p.cff(k),k);
      begin
        Standard_CSeries_Polynomials.Add(res,tpk);
        Standard_CSeries_Polynomials.Clear(tpk);
      end;
    end loop;
    return res;
  end Create;

-- DESTRUCTOR :

  procedure Clear ( p : in out Standard_Polynomial_CSeries.Poly ) is
  begin
    for k in 0..p.deg loop
      Standard_Complex_Polynomials.Clear(p.cff(k));
    end loop;
  end Clear;

end Standard_Polynomial_CSeries;
