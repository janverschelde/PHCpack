with Penta_Double_Numbers;              use Penta_Double_Numbers;
with PentDobl_Complex_Numbers;          use PentDobl_Complex_Numbers;
with PentDobl_Complex_Polynomials;
with PentDobl_Complex_Series;

package body PentDobl_Polynomial_CSeries is

-- CONSTRUCTORS :

  function Create ( p : PentDobl_CSeries_Polynomials.Poly )
                  return PentDobl_Polynomial_CSeries.Poly is

    maxdeg : constant integer32 := 16;
    res : PentDobl_Polynomial_CSeries.Poly(maxdeg);

    procedure Add_Term ( t : in PentDobl_CSeries_Polynomials.Term;
                         cont : out boolean ) is

      cf : constant PentDobl_Complex_Series.Link_to_Series := t.cf;
   
    begin
      for k in 0..cf.deg loop
        declare
          rt : PentDobl_Complex_Polynomials.Term;
        begin
          rt.cf := cf.cff(k);
          rt.dg := PentDobl_Complex_Polynomials.Degrees(t.dg);
          PentDobl_Complex_Polynomials.Add(res.cff(k),rt);
        end;
      end loop;
      cont := true;
    end Add_Term;
    procedure Add_Terms is
      new PentDobl_CSeries_Polynomials.Visiting_Iterator(Add_Term);

  begin
    res.cff := (res.cff'range => PentDobl_Complex_Polynomials.Null_Poly);
    Add_Terms(p);
    return res;
  end Create;

  function Convert ( p : PentDobl_Complex_Polynomials.Poly;
                     d : integer32 )
                   return PentDobl_CSeries_Polynomials.Poly is

  -- DESCRIPTION :
  --   Each coefficient c of the polynomial p is stored as c*t^d
  --   as series coefficient in the polynomial on return.

    res : PentDobl_CSeries_Polynomials.Poly
        := PentDobl_CSeries_Polynomials.Null_Poly;
    pd_zero : constant penta_double := create(0.0);
    zero : constant Complex_Number := create(pd_zero);

    procedure Add_Term ( t : in PentDobl_Complex_Polynomials.Term;
                         cont : out boolean ) is

      rt : PentDobl_CSeries_Polynomials.Term;

    begin
      rt.cf := PentDobl_Complex_Series.Create(0,d);
      for k in 0..d-1 loop
        rt.cf.cff(k) := zero;
      end loop;
      rt.cf.cff(d) := t.cf;
      rt.dg := PentDobl_CSeries_Polynomials.Degrees(t.dg);
      PentDobl_CSeries_Polynomials.Add(res,rt);
      cont := true;
    end Add_Term;
    procedure Add_Terms is new
      PentDobl_Complex_Polynomials.Visiting_Iterator(Add_Term);

  begin
    Add_Terms(p);
    return res;
  end Convert;

  function Create ( p : PentDobl_Polynomial_CSeries.Poly )
                  return PentDobl_CSeries_Polynomials.Poly is

    res : PentDobl_CSeries_Polynomials.Poly := Convert(p.cff(0),0);

  begin
    for k in 1..p.deg loop
      declare
        tpk : PentDobl_CSeries_Polynomials.Poly := Convert(p.cff(k),k);
      begin
        PentDobl_CSeries_Polynomials.Add(res,tpk);
        PentDobl_CSeries_Polynomials.Clear(tpk);
      end;
    end loop;
    return res;
  end Create;

-- DESTRUCTOR :

  procedure Clear ( p : in out PentDobl_Polynomial_CSeries.Poly ) is
  begin
    for k in 0..p.deg loop
      PentDobl_Complex_Polynomials.Clear(p.cff(k));
    end loop;
  end Clear;

end PentDobl_Polynomial_CSeries;
