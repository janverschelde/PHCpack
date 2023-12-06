with Hexa_Double_Numbers;               use Hexa_Double_Numbers;
with HexaDobl_Complex_Numbers;          use HexaDobl_Complex_Numbers;
with HexaDobl_Complex_Polynomials;
with HexaDobl_Complex_Series;

package body HexaDobl_Polynomial_CSeries is

-- CONSTRUCTORS :

  function Create ( p : HexaDobl_CSeries_Polynomials.Poly )
                  return HexaDobl_Polynomial_CSeries.Poly is

    maxdeg : constant integer32 := 16;
    res : HexaDobl_Polynomial_CSeries.Poly(maxdeg);

    procedure Add_Term ( t : in HexaDobl_CSeries_Polynomials.Term;
                         cont : out boolean ) is

      cf : constant HexaDobl_Complex_Series.Link_to_Series := t.cf;
   
    begin
      for k in 0..cf.deg loop
        declare
          rt : HexaDobl_Complex_Polynomials.Term;
        begin
          rt.cf := cf.cff(k);
          rt.dg := HexaDobl_Complex_Polynomials.Degrees(t.dg);
          HexaDobl_Complex_Polynomials.Add(res.cff(k),rt);
        end;
      end loop;
      cont := true;
    end Add_Term;
    procedure Add_Terms is
      new HexaDobl_CSeries_Polynomials.Visiting_Iterator(Add_Term);

  begin
    res.cff := (res.cff'range => HexaDobl_Complex_Polynomials.Null_Poly);
    Add_Terms(p);
    return res;
  end Create;

  function Convert ( p : HexaDobl_Complex_Polynomials.Poly;
                     d : integer32 )
                   return HexaDobl_CSeries_Polynomials.Poly is

  -- DESCRIPTION :
  --   Each coefficient c of the polynomial p is stored as c*t^d
  --   as series coefficient in the polynomial on return.

    res : HexaDobl_CSeries_Polynomials.Poly
        := HexaDobl_CSeries_Polynomials.Null_Poly;
    da_zero : constant hexa_double := create(0.0);
    zero : constant Complex_Number := create(da_zero);

    procedure Add_Term ( t : in HexaDobl_Complex_Polynomials.Term;
                         cont : out boolean ) is

      rt : HexaDobl_CSeries_Polynomials.Term;

    begin
      rt.cf := HexaDobl_Complex_Series.Create(0,d);
      for k in 0..d-1 loop
        rt.cf.cff(k) := zero;
      end loop;
      rt.cf.cff(d) := t.cf;
      rt.dg := HexaDobl_CSeries_Polynomials.Degrees(t.dg);
      HexaDobl_CSeries_Polynomials.Add(res,rt);
      cont := true;
    end Add_Term;
    procedure Add_Terms is new
      HexaDobl_Complex_Polynomials.Visiting_Iterator(Add_Term);

  begin
    Add_Terms(p);
    return res;
  end Convert;

  function Create ( p : HexaDobl_Polynomial_CSeries.Poly )
                  return HexaDobl_CSeries_Polynomials.Poly is

    res : HexaDobl_CSeries_Polynomials.Poly := Convert(p.cff(0),0);

  begin
    for k in 1..p.deg loop
      declare
        tpk : HexaDobl_CSeries_Polynomials.Poly := Convert(p.cff(k),k);
      begin
        HexaDobl_CSeries_Polynomials.Add(res,tpk);
        HexaDobl_CSeries_Polynomials.Clear(tpk);
      end;
    end loop;
    return res;
  end Create;

-- DESTRUCTOR :

  procedure Clear ( p : in out HexaDobl_Polynomial_CSeries.Poly ) is
  begin
    for k in 0..p.deg loop
      HexaDobl_Complex_Polynomials.Clear(p.cff(k));
    end loop;
  end Clear;

end HexaDobl_Polynomial_CSeries;
