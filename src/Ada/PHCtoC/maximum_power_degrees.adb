with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;

package body Maximum_Power_Degrees is

  function Maximum_Power
             ( t : Standard_Complex_Polynomials.Term ) return integer32 is

    use Standard_Complex_Polynomials;

    res : integer32 := -1;

  begin
    if t.dg /= null then
      res := integer32(t.dg(t.dg'first));
      for k in t.dg'first+1..t.dg'last loop
        if integer32(t.dg(k)) > res
         then res := integer32(t.dg(k));
        end if;
      end loop;
    end if;
    return res;
  end Maximum_Power;

  function Maximum_Power
             ( t : DoblDobl_Complex_Polynomials.Term ) return integer32 is

    use DoblDobl_Complex_Polynomials;

    res : integer32 := -1;

  begin
    if t.dg /= null then
      res := integer32(t.dg(t.dg'first));
      for k in t.dg'first+1..t.dg'last loop
        if integer32(t.dg(k)) > res
         then res := integer32(t.dg(k));
        end if;
      end loop;
    end if;
    return res;
  end Maximum_Power;

  function Maximum_Power
             ( t : QuadDobl_Complex_Polynomials.Term ) return integer32 is

    use QuadDobl_Complex_Polynomials;

    res : integer32 := -1;

  begin
    if t.dg /= null then
      res := integer32(t.dg(t.dg'first));
      for k in t.dg'first+1..t.dg'last loop
        if integer32(t.dg(k)) > res
         then res := integer32(t.dg(k));
        end if;
      end loop;
    end if;
    return res;
  end Maximum_Power;

  function Maximum_Power
             ( p : Standard_Complex_Polynomials.Poly ) return integer32 is

    res : integer32 := -1;

    procedure Visit_Term ( t : in Standard_Complex_Polynomials.Term;
                           continue : out boolean ) is

      deg : constant integer32 := Maximum_Power(t);

    begin
      if deg > res
       then res := deg;
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Maximum_Power;

  function Maximum_Power
             ( p : DoblDobl_Complex_Polynomials.Poly ) return integer32 is

    res : integer32 := -1;

    procedure Visit_Term ( t : in DoblDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

      deg : constant integer32 := Maximum_Power(t);

    begin
      if deg > res
       then res := deg;
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is
      new DoblDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Maximum_Power;

  function Maximum_Power
             ( p : QuadDobl_Complex_Polynomials.Poly ) return integer32 is

    res : integer32 := -1;

    procedure Visit_Term ( t : in QuadDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

      deg : constant integer32 := Maximum_Power(t);

    begin
      if deg > res
       then res := deg;
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is
      new QuadDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Maximum_Power;

  function Maximum_Power
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return integer32 is

    res : integer32 := Maximum_Power(p(p'first));
    deg : integer32;

  begin
    for k in p'first+1..p'last loop
      deg := Maximum_Power(p(k));
      if deg > res
       then res := deg;
      end if;
    end loop;
    return res;
  end Maximum_Power;

  function Maximum_Power
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return integer32 is

    res : integer32 := Maximum_Power(p(p'first));
    deg : integer32;

  begin
    for k in p'first+1..p'last loop
      deg := Maximum_Power(p(k));
      if deg > res
       then res := deg;
      end if;
    end loop;
    return res;
  end Maximum_Power;

  function Maximum_Power
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return integer32 is

    res : integer32 := Maximum_Power(p(p'first));
    deg : integer32;

  begin
    for k in p'first+1..p'last loop
      deg := Maximum_Power(p(k));
      if deg > res
       then res := deg;
      end if;
    end loop;
    return res;
  end Maximum_Power;

end Maximum_Power_Degrees;
