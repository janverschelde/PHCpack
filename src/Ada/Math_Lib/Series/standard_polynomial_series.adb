with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Polynomials;

package body Standard_Polynomial_Series is

-- CONSTRUCTORS :

  function Create ( p : Standard_Series_Polynomials.Poly )
                  return Standard_Polynomial_Series.Poly is

    res : Standard_Polynomial_Series.Poly;

    procedure Add_Term ( t : in Standard_Series_Polynomials.Term;
                         cont : out boolean ) is

      cf : constant Standard_Dense_Series.Series := t.cf;
   
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
      if cf.deg > res.deg
       then res.deg := cf.deg;
      end if;
      cont := true;
    end Add_Term;
    procedure Add_Terms is
      new Standard_Series_Polynomials.Visiting_Iterator(Add_Term);

  begin
    res.deg := -1;
    res.cff := (res.cff'range => Standard_Complex_Polynomials.Null_Poly);
    Add_Terms(p);
    return res;
  end Create;

  function Convert ( p : Standard_Complex_Polynomials.Poly;
                     d : integer32 )
                   return Standard_Series_Polynomials.Poly is

  -- DESCRIPTION :
  --   Each coefficient c of the polynomial p is stored as c*t^d
  --   as series coefficient in the polynomial on return.

    res : Standard_Series_Polynomials.Poly
        := Standard_Series_Polynomials.Null_Poly;

    procedure Add_Term ( t : in Standard_Complex_Polynomials.Term;
                         cont : out boolean ) is

      rt : Standard_Series_Polynomials.Term;

    begin
      rt.cf.deg := d;
      for k in 0..d-1 loop
        rt.cf.cff(k) := Create(0.0);
      end loop;
      rt.cf.cff(d) := t.cf;
      rt.dg := Standard_Series_Polynomials.Degrees(t.dg);
      Standard_Series_Polynomials.Add(res,rt);
      cont := true;
    end Add_Term;
    procedure Add_Terms is new
      Standard_Complex_Polynomials.Visiting_Iterator(Add_Term);

  begin
    Add_Terms(p);
    return res;
  end Convert;

  function Create ( p : Standard_Polynomial_Series.Poly )
                  return Standard_Series_Polynomials.Poly is

    res : Standard_Series_Polynomials.Poly := Convert(p.cff(0),0);

  begin
    for k in 1..p.deg loop
      declare
        tpk : Standard_Series_Polynomials.Poly := Convert(p.cff(k),k);
      begin
        Standard_Series_Polynomials.Add(res,tpk);
        Standard_Series_Polynomials.Clear(tpk);
      end;
    end loop;
    return res;
  end Create;

-- DESTRUCTOR :

  procedure Clear ( p : in out Standard_Polynomial_Series.Poly ) is
  begin
    for k in 0..p.deg loop
      Standard_Complex_Polynomials.Clear(p.cff(k));
    end loop;
    p.deg := -1;
  end Clear;

end Standard_Polynomial_Series;
