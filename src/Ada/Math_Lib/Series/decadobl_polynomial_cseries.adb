with Deca_Double_Numbers;               use Deca_Double_Numbers;
with DecaDobl_Complex_Numbers;          use DecaDobl_Complex_Numbers;
with DecaDobl_Complex_Polynomials;
with DecaDobl_Complex_Series;

package body DecaDobl_Polynomial_CSeries is

-- CONSTRUCTORS :

  function Create ( p : DecaDobl_CSeries_Polynomials.Poly )
                  return DecaDobl_Polynomial_CSeries.Poly is

    maxdeg : constant integer32 := 16;
    res : DecaDobl_Polynomial_CSeries.Poly(maxdeg);

    procedure Add_Term ( t : in DecaDobl_CSeries_Polynomials.Term;
                         cont : out boolean ) is

      cf : constant DecaDobl_Complex_Series.Link_to_Series := t.cf;
   
    begin
      for k in 0..cf.deg loop
        declare
          rt : DecaDobl_Complex_Polynomials.Term;
        begin
          rt.cf := cf.cff(k);
          rt.dg := DecaDobl_Complex_Polynomials.Degrees(t.dg);
          DecaDobl_Complex_Polynomials.Add(res.cff(k),rt);
        end;
      end loop;
      cont := true;
    end Add_Term;
    procedure Add_Terms is
      new DecaDobl_CSeries_Polynomials.Visiting_Iterator(Add_Term);

  begin
    res.cff := (res.cff'range => DecaDobl_Complex_Polynomials.Null_Poly);
    Add_Terms(p);
    return res;
  end Create;

  function Convert ( p : DecaDobl_Complex_Polynomials.Poly;
                     d : integer32 )
                   return DecaDobl_CSeries_Polynomials.Poly is

  -- DESCRIPTION :
  --   Each coefficient c of the polynomial p is stored as c*t^d
  --   as series coefficient in the polynomial on return.

    res : DecaDobl_CSeries_Polynomials.Poly
        := DecaDobl_CSeries_Polynomials.Null_Poly;
    da_zero : constant deca_double := create(0.0);
    zero : constant Complex_Number := create(da_zero);

    procedure Add_Term ( t : in DecaDobl_Complex_Polynomials.Term;
                         cont : out boolean ) is

      rt : DecaDobl_CSeries_Polynomials.Term;

    begin
      rt.cf := DecaDobl_Complex_Series.Create(0,d);
      for k in 0..d-1 loop
        rt.cf.cff(k) := zero;
      end loop;
      rt.cf.cff(d) := t.cf;
      rt.dg := DecaDobl_CSeries_Polynomials.Degrees(t.dg);
      DecaDobl_CSeries_Polynomials.Add(res,rt);
      cont := true;
    end Add_Term;
    procedure Add_Terms is new
      DecaDobl_Complex_Polynomials.Visiting_Iterator(Add_Term);

  begin
    Add_Terms(p);
    return res;
  end Convert;

  function Create ( p : DecaDobl_Polynomial_CSeries.Poly )
                  return DecaDobl_CSeries_Polynomials.Poly is

    res : DecaDobl_CSeries_Polynomials.Poly := Convert(p.cff(0),0);

  begin
    for k in 1..p.deg loop
      declare
        tpk : DecaDobl_CSeries_Polynomials.Poly := Convert(p.cff(k),k);
      begin
        DecaDobl_CSeries_Polynomials.Add(res,tpk);
        DecaDobl_CSeries_Polynomials.Clear(tpk);
      end;
    end loop;
    return res;
  end Create;

-- DESTRUCTOR :

  procedure Clear ( p : in out DecaDobl_Polynomial_CSeries.Poly ) is
  begin
    for k in 0..p.deg loop
      DecaDobl_Complex_Polynomials.Clear(p.cff(k));
    end loop;
  end Clear;

end DecaDobl_Polynomial_CSeries;
