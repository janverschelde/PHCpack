with Double_Double_Numbers;             use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Polynomials;

package body DoblDobl_Polynomial_Series is

-- CONSTRUCTORS :

  function Create ( p : DoblDobl_Series_Polynomials.Poly )
                  return DoblDobl_Polynomial_Series.Poly is

    res : DoblDobl_Polynomial_Series.Poly;

    procedure Add_Term ( t : in DoblDobl_Series_Polynomials.Term;
                         cont : out boolean ) is

      cf : constant DoblDobl_Dense_Series.Series := t.cf;
   
    begin
      for k in 0..cf.deg loop
        declare
          rt : DoblDobl_Complex_Polynomials.Term;
        begin
          rt.cf := cf.cff(k);
          rt.dg := DoblDobl_Complex_Polynomials.Degrees(t.dg);
          DoblDobl_Complex_Polynomials.Add(res.cff(k),rt);
        end;
      end loop;
      if cf.deg > res.deg
       then res.deg := cf.deg;
      end if;
      cont := true;
    end Add_Term;
    procedure Add_Terms is
      new DoblDobl_Series_Polynomials.Visiting_Iterator(Add_Term);

  begin
    res.deg := -1;
    res.cff := (res.cff'range => DoblDobl_Complex_Polynomials.Null_Poly);
    Add_Terms(p);
    return res;
  end Create;

  function Convert ( p : DoblDobl_Complex_Polynomials.Poly;
                     d : integer32 )
                   return DoblDobl_Series_Polynomials.Poly is

  -- DESCRIPTION :
  --   Each coefficient c of the polynomial p is stored as c*t^d
  --   as series coefficient in the polynomial on return.

    res : DoblDobl_Series_Polynomials.Poly
        := DoblDobl_Series_Polynomials.Null_Poly;
    dd_zero : constant double_double := create(0.0);
    zero : constant Complex_Number := create(dd_zero);

    procedure Add_Term ( t : in DoblDobl_Complex_Polynomials.Term;
                         cont : out boolean ) is

      rt : DoblDobl_Series_Polynomials.Term;

    begin
      rt.cf.deg := d;
      for k in 0..d-1 loop
        rt.cf.cff(k) := zero;
      end loop;
      rt.cf.cff(d) := t.cf;
      rt.dg := DoblDobl_Series_Polynomials.Degrees(t.dg);
      DoblDobl_Series_Polynomials.Add(res,rt);
      cont := true;
    end Add_Term;
    procedure Add_Terms is new
      DoblDobl_Complex_Polynomials.Visiting_Iterator(Add_Term);

  begin
    Add_Terms(p);
    return res;
  end Convert;

  function Create ( p : DoblDobl_Polynomial_Series.Poly )
                  return DoblDobl_Series_Polynomials.Poly is

    res : DoblDobl_Series_Polynomials.Poly := Convert(p.cff(0),0);

  begin
    for k in 1..p.deg loop
      declare
        tpk : DoblDobl_Series_Polynomials.Poly := Convert(p.cff(k),k);
      begin
        DoblDobl_Series_Polynomials.Add(res,tpk);
        DoblDobl_Series_Polynomials.Clear(tpk);
      end;
    end loop;
    return res;
  end Create;

-- DESTRUCTOR :

  procedure Clear ( p : in out DoblDobl_Polynomial_Series.Poly ) is
  begin
    for k in 0..p.deg loop
      DoblDobl_Complex_Polynomials.Clear(p.cff(k));
    end loop;
    p.deg := -1;
  end Clear;

end DoblDobl_Polynomial_Series;
