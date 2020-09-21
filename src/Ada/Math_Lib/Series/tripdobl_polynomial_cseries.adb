with Triple_Double_Numbers;             use Triple_Double_Numbers;
with TripDobl_Complex_Numbers;          use TripDobl_Complex_Numbers;
with TripDobl_Complex_Polynomials;
with TripDobl_Complex_Series;

package body TripDobl_Polynomial_CSeries is

-- CONSTRUCTORS :

  function Create ( p : TripDobl_CSeries_Polynomials.Poly )
                  return TripDobl_Polynomial_CSeries.Poly is

    maxdeg : constant integer32 := 16;
    res : TripDobl_Polynomial_CSeries.Poly(maxdeg);

    procedure Add_Term ( t : in TripDobl_CSeries_Polynomials.Term;
                         cont : out boolean ) is

      cf : constant TripDobl_Complex_Series.Link_to_Series := t.cf;
   
    begin
      for k in 0..cf.deg loop
        declare
          rt : TripDobl_Complex_Polynomials.Term;
        begin
          rt.cf := cf.cff(k);
          rt.dg := TripDobl_Complex_Polynomials.Degrees(t.dg);
          TripDobl_Complex_Polynomials.Add(res.cff(k),rt);
        end;
      end loop;
      cont := true;
    end Add_Term;
    procedure Add_Terms is
      new TripDobl_CSeries_Polynomials.Visiting_Iterator(Add_Term);

  begin
    res.cff := (res.cff'range => TripDobl_Complex_Polynomials.Null_Poly);
    Add_Terms(p);
    return res;
  end Create;

  function Convert ( p : TripDobl_Complex_Polynomials.Poly;
                     d : integer32 )
                   return TripDobl_CSeries_Polynomials.Poly is

  -- DESCRIPTION :
  --   Each coefficient c of the polynomial p is stored as c*t^d
  --   as series coefficient in the polynomial on return.

    res : TripDobl_CSeries_Polynomials.Poly
        := TripDobl_CSeries_Polynomials.Null_Poly;
    td_zero : constant triple_double := create(0.0);
    zero : constant Complex_Number := create(td_zero);

    procedure Add_Term ( t : in TripDobl_Complex_Polynomials.Term;
                         cont : out boolean ) is

      rt : TripDobl_CSeries_Polynomials.Term;

    begin
      rt.cf := TripDobl_Complex_Series.Create(0,d);
      for k in 0..d-1 loop
        rt.cf.cff(k) := zero;
      end loop;
      rt.cf.cff(d) := t.cf;
      rt.dg := TripDobl_CSeries_Polynomials.Degrees(t.dg);
      TripDobl_CSeries_Polynomials.Add(res,rt);
      cont := true;
    end Add_Term;
    procedure Add_Terms is new
      TripDobl_Complex_Polynomials.Visiting_Iterator(Add_Term);

  begin
    Add_Terms(p);
    return res;
  end Convert;

  function Create ( p : TripDobl_Polynomial_CSeries.Poly )
                  return TripDobl_CSeries_Polynomials.Poly is

    res : TripDobl_CSeries_Polynomials.Poly := Convert(p.cff(0),0);

  begin
    for k in 1..p.deg loop
      declare
        tpk : TripDobl_CSeries_Polynomials.Poly := Convert(p.cff(k),k);
      begin
        TripDobl_CSeries_Polynomials.Add(res,tpk);
        TripDobl_CSeries_Polynomials.Clear(tpk);
      end;
    end loop;
    return res;
  end Create;

-- DESTRUCTOR :

  procedure Clear ( p : in out TripDobl_Polynomial_CSeries.Poly ) is
  begin
    for k in 0..p.deg loop
      TripDobl_Complex_Polynomials.Clear(p.cff(k));
    end loop;
  end Clear;

end TripDobl_Polynomial_CSeries;
