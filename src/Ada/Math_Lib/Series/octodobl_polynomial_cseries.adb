with Octo_Double_Numbers;               use Octo_Double_Numbers;
with OctoDobl_Complex_Numbers;          use OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Polynomials;
with OctoDobl_Complex_Series;

package body OctoDobl_Polynomial_CSeries is

-- CONSTRUCTORS :

  function Create ( p : OctoDobl_CSeries_Polynomials.Poly )
                  return OctoDobl_Polynomial_CSeries.Poly is

    maxdeg : constant integer32 := 16;
    res : OctoDobl_Polynomial_CSeries.Poly(maxdeg);

    procedure Add_Term ( t : in OctoDobl_CSeries_Polynomials.Term;
                         cont : out boolean ) is

      cf : constant OctoDobl_Complex_Series.Link_to_Series := t.cf;
   
    begin
      for k in 0..cf.deg loop
        declare
          rt : OctoDobl_Complex_Polynomials.Term;
        begin
          rt.cf := cf.cff(k);
          rt.dg := OctoDobl_Complex_Polynomials.Degrees(t.dg);
          OctoDobl_Complex_Polynomials.Add(res.cff(k),rt);
        end;
      end loop;
      cont := true;
    end Add_Term;
    procedure Add_Terms is
      new OctoDobl_CSeries_Polynomials.Visiting_Iterator(Add_Term);

  begin
    res.cff := (res.cff'range => OctoDobl_Complex_Polynomials.Null_Poly);
    Add_Terms(p);
    return res;
  end Create;

  function Convert ( p : OctoDobl_Complex_Polynomials.Poly;
                     d : integer32 )
                   return OctoDobl_CSeries_Polynomials.Poly is

  -- DESCRIPTION :
  --   Each coefficient c of the polynomial p is stored as c*t^d
  --   as series coefficient in the polynomial on return.

    res : OctoDobl_CSeries_Polynomials.Poly
        := OctoDobl_CSeries_Polynomials.Null_Poly;
    od_zero : constant octo_double := create(0.0);
    zero : constant Complex_Number := create(od_zero);

    procedure Add_Term ( t : in OctoDobl_Complex_Polynomials.Term;
                         cont : out boolean ) is

      rt : OctoDobl_CSeries_Polynomials.Term;

    begin
      rt.cf := OctoDobl_Complex_Series.Create(0,d);
      for k in 0..d-1 loop
        rt.cf.cff(k) := zero;
      end loop;
      rt.cf.cff(d) := t.cf;
      rt.dg := OctoDobl_CSeries_Polynomials.Degrees(t.dg);
      OctoDobl_CSeries_Polynomials.Add(res,rt);
      cont := true;
    end Add_Term;
    procedure Add_Terms is new
      OctoDobl_Complex_Polynomials.Visiting_Iterator(Add_Term);

  begin
    Add_Terms(p);
    return res;
  end Convert;

  function Create ( p : OctoDobl_Polynomial_CSeries.Poly )
                  return OctoDobl_CSeries_Polynomials.Poly is

    res : OctoDobl_CSeries_Polynomials.Poly := Convert(p.cff(0),0);

  begin
    for k in 1..p.deg loop
      declare
        tpk : OctoDobl_CSeries_Polynomials.Poly := Convert(p.cff(k),k);
      begin
        OctoDobl_CSeries_Polynomials.Add(res,tpk);
        OctoDobl_CSeries_Polynomials.Clear(tpk);
      end;
    end loop;
    return res;
  end Create;

-- DESTRUCTOR :

  procedure Clear ( p : in out OctoDobl_Polynomial_CSeries.Poly ) is
  begin
    for k in 0..p.deg loop
      OctoDobl_Complex_Polynomials.Clear(p.cff(k));
    end loop;
  end Clear;

end OctoDobl_Polynomial_CSeries;
