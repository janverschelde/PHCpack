with Quad_Double_Numbers;               use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Series;

package body QuadDobl_Polynomial_CSeries is

-- CONSTRUCTORS :

  function Create ( p : QuadDobl_CSeries_Polynomials.Poly )
                  return QuadDobl_Polynomial_CSeries.Poly is

    maxdeg : constant integer32 := 16;
    res : QuadDobl_Polynomial_CSeries.Poly(maxdeg);

    procedure Add_Term ( t : in QuadDobl_CSeries_Polynomials.Term;
                         cont : out boolean ) is

      cf : constant QuadDobl_Complex_Series.Link_to_Series := t.cf;
   
    begin
      for k in 0..cf.deg loop
        declare
          rt : QuadDobl_Complex_Polynomials.Term;
        begin
          rt.cf := cf.cff(k);
          rt.dg := QuadDobl_Complex_Polynomials.Degrees(t.dg);
          QuadDobl_Complex_Polynomials.Add(res.cff(k),rt);
        end;
      end loop;
      cont := true;
    end Add_Term;
    procedure Add_Terms is
      new QuadDobl_CSeries_Polynomials.Visiting_Iterator(Add_Term);

  begin
    res.cff := (res.cff'range => QuadDobl_Complex_Polynomials.Null_Poly);
    Add_Terms(p);
    return res;
  end Create;

  function Convert ( p : QuadDobl_Complex_Polynomials.Poly;
                     d : integer32 )
                   return QuadDobl_CSeries_Polynomials.Poly is

  -- DESCRIPTION :
  --   Each coefficient c of the polynomial p is stored as c*t^d
  --   as series coefficient in the polynomial on return.

    res : QuadDobl_CSeries_Polynomials.Poly
        := QuadDobl_CSeries_Polynomials.Null_Poly;
    qd_zero : constant quad_double := create(0.0);
    zero : constant Complex_Number := create(qd_zero);

    procedure Add_Term ( t : in QuadDobl_Complex_Polynomials.Term;
                         cont : out boolean ) is

      rt : QuadDobl_CSeries_Polynomials.Term;

    begin
      rt.cf := QuadDobl_Complex_Series.Create(0,d);
      for k in 0..d-1 loop
        rt.cf.cff(k) := zero;
      end loop;
      rt.cf.cff(d) := t.cf;
      rt.dg := QuadDobl_CSeries_Polynomials.Degrees(t.dg);
      QuadDobl_CSeries_Polynomials.Add(res,rt);
      cont := true;
    end Add_Term;
    procedure Add_Terms is new
      QuadDobl_Complex_Polynomials.Visiting_Iterator(Add_Term);

  begin
    Add_Terms(p);
    return res;
  end Convert;

  function Create ( p : QuadDobl_Polynomial_CSeries.Poly )
                  return QuadDobl_CSeries_Polynomials.Poly is

    res : QuadDobl_CSeries_Polynomials.Poly := Convert(p.cff(0),0);

  begin
    for k in 1..p.deg loop
      declare
        tpk : QuadDobl_CSeries_Polynomials.Poly := Convert(p.cff(k),k);
      begin
        QuadDobl_CSeries_Polynomials.Add(res,tpk);
        QuadDobl_CSeries_Polynomials.Clear(tpk);
      end;
    end loop;
    return res;
  end Create;

-- DESTRUCTOR :

  procedure Clear ( p : in out QuadDobl_Polynomial_CSeries.Poly ) is
  begin
    for k in 0..p.deg loop
      QuadDobl_Complex_Polynomials.Clear(p.cff(k));
    end loop;
  end Clear;

end QuadDobl_Polynomial_CSeries;
