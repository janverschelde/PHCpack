with Quad_Double_Numbers;               use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Polynomials;

package body QuadDobl_Polynomial_Series is

-- CONSTRUCTORS :

  function Create ( p : QuadDobl_Series_Polynomials.Poly )
                  return QuadDobl_Polynomial_Series.Poly is

    res : QuadDobl_Polynomial_Series.Poly;

    procedure Add_Term ( t : in QuadDobl_Series_Polynomials.Term;
                         cont : out boolean ) is

      cf : constant QuadDobl_Dense_Series.Series := t.cf;
   
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
      if cf.deg > res.deg
       then res.deg := cf.deg;
      end if;
      cont := true;
    end Add_Term;
    procedure Add_Terms is
      new QuadDobl_Series_Polynomials.Visiting_Iterator(Add_Term);

  begin
    res.deg := -1;
    res.cff := (res.cff'range => QuadDobl_Complex_Polynomials.Null_Poly);
    Add_Terms(p);
    return res;
  end Create;

  function Convert ( p : QuadDobl_Complex_Polynomials.Poly;
                     d : integer32 )
                   return QuadDobl_Series_Polynomials.Poly is

  -- DESCRIPTION :
  --   Each coefficient c of the polynomial p is stored as c*t^d
  --   as series coefficient in the polynomial on return.

    res : QuadDobl_Series_Polynomials.Poly
        := QuadDobl_Series_Polynomials.Null_Poly;
    qd_zero : constant quad_double := create(0.0);
    zero : constant Complex_Number := create(qd_zero);

    procedure Add_Term ( t : in QuadDobl_Complex_Polynomials.Term;
                         cont : out boolean ) is

      rt : QuadDobl_Series_Polynomials.Term;

    begin
      rt.cf.deg := d;
      for k in 0..d-1 loop
        rt.cf.cff(k) := zero;
      end loop;
      rt.cf.cff(d) := t.cf;
      rt.dg := QuadDobl_Series_Polynomials.Degrees(t.dg);
      QuadDobl_Series_Polynomials.Add(res,rt);
      cont := true;
    end Add_Term;
    procedure Add_Terms is new
      QuadDobl_Complex_Polynomials.Visiting_Iterator(Add_Term);

  begin
    Add_Terms(p);
    return res;
  end Convert;

  function Create ( p : QuadDobl_Polynomial_Series.Poly )
                  return QuadDobl_Series_Polynomials.Poly is

    res : QuadDobl_Series_Polynomials.Poly := Convert(p.cff(0),0);

  begin
    for k in 1..p.deg loop
      declare
        tpk : QuadDobl_Series_Polynomials.Poly := Convert(p.cff(k),k);
      begin
        QuadDobl_Series_Polynomials.Add(res,tpk);
        QuadDobl_Series_Polynomials.Clear(tpk);
      end;
    end loop;
    return res;
  end Create;

-- DESTRUCTOR :

  procedure Clear ( p : in out QuadDobl_Polynomial_Series.Poly ) is
  begin
    for k in 0..p.deg loop
      QuadDobl_Complex_Polynomials.Clear(p.cff(k));
    end loop;
    p.deg := -1;
  end Clear;

end QuadDobl_Polynomial_Series;
