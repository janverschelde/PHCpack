with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;

package body System_Vector_Convertors is

  function Convert ( t : Standard_Complex_Polynomials.Term )
                   return Standard_Complex_Monomials.Monomial is

    res : constant Standard_Complex_Monomials.Monomial
        := Standard_Complex_Monomials.Create(t.cf,t.dg.all);

  begin
    return res;
  end Convert;

  function Convert ( t : DoblDobl_Complex_Polynomials.Term )
                   return DoblDobl_Complex_Monomials.Monomial is

    res : constant DoblDobl_Complex_Monomials.Monomial
        := DoblDobl_Complex_Monomials.Create(t.cf,t.dg.all);

  begin
    return res;
  end Convert;

  function Convert ( t : QuadDobl_Complex_Polynomials.Term )
                   return QuadDobl_Complex_Monomials.Monomial is

    res : constant QuadDobl_Complex_Monomials.Monomial
        := QuadDobl_Complex_Monomials.Create(t.cf,t.dg.all);

  begin
    return res;
  end Convert;

  function Is_Zero ( v : Standard_Natural_Vectors.Vector )
                   return boolean is

  begin
    for i in v'range loop
      if v(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  function Count_Constant
             ( p : Standard_Complex_Polynomials.Poly ) return integer32 is

    res : integer32 := 0;

    procedure Visit_Term ( t : in Standard_Complex_Polynomials.Term;
                           continue : out boolean ) is

    begin
      if Is_Zero(t.dg.all)
       then continue := false; res := 1;
       else continue := true;
      end if;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Count_Constant;

  function Count_Constant
             ( p : DoblDobl_Complex_Polynomials.Poly ) return integer32 is

    res : integer32 := 0;

    procedure Visit_Term ( t : in DoblDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

    begin
      if Is_Zero(t.dg.all)
       then continue := false; res := 1;
       else continue := true;
      end if;
    end Visit_Term;
    procedure Visit_Terms is
      new DoblDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Count_Constant;

  function Count_Constant
             ( p : QuadDobl_Complex_Polynomials.Poly ) return integer32 is

    res : integer32 := 0;

    procedure Visit_Term ( t : in QuadDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

    begin
      if Is_Zero(t.dg.all)
       then continue := false; res := 1;
       else continue := true;
      end if;
    end Visit_Term;
    procedure Visit_Terms is
      new QuadDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Count_Constant;

  function Convert ( p : Standard_Complex_Polynomials.Poly )
                   return Standard_Monomial_Vectors.Polynomial is

    cnt : constant integer32 := Count_Constant(p);
    nbr : constant integer32
        := integer32(Standard_Complex_Polynomials.Number_of_Terms(p));
    dim : constant integer32
        := integer32(Standard_Complex_Polynomials.Number_of_Unknowns(p));
    res : Standard_Monomial_Vectors.Polynomial(dim,nbr-cnt);
    idx : integer32 := 0;

    procedure Visit_Term ( t : in Standard_Complex_Polynomials.Term;
                           continue : out boolean ) is

    begin
      if Is_Zero(t.dg.all) then
        res.cff0 := t.cf;
      else
        idx := idx + 1;
        res.mons(idx) := Standard_Complex_Monomials.Create(t.cf,t.dg.all);
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    res.cff0 := Standard_Complex_Numbers.Create(0.0);
    Visit_Terms(p);
    Standard_Monomial_Vectors.Power_Update(res);
    return res;
  end Convert;

  function Convert ( p : DoblDobl_Complex_Polynomials.Poly )
                   return DoblDobl_Monomial_Vectors.Polynomial is

    cnt : constant integer32 := Count_Constant(p);
    nbr : constant integer32
        := integer32(DoblDobl_Complex_Polynomials.Number_of_Terms(p));
    dim : constant integer32
        := integer32(DoblDobl_Complex_Polynomials.Number_of_Unknowns(p));
    res : DoblDobl_Monomial_Vectors.Polynomial(dim,nbr-cnt);
    idx : integer32 := 0;

    procedure Visit_Term ( t : in DoblDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

    begin
      if Is_Zero(t.dg.all) then
        res.cff0 := t.cf;
      else
        idx := idx + 1;
        res.mons(idx) := DoblDobl_Complex_Monomials.Create(t.cf,t.dg.all);
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is
      new DoblDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    res.cff0 := DoblDobl_Complex_Numbers.Create(integer32(0));
    Visit_Terms(p);
    DoblDobl_Monomial_Vectors.Power_Update(res);
    return res;
  end Convert;

  function Convert ( p : QuadDobl_Complex_Polynomials.Poly )
                   return QuadDobl_Monomial_Vectors.Polynomial is

    cnt : constant integer32 := Count_Constant(p);
    nbr : constant integer32
        := integer32(QuadDobl_Complex_Polynomials.Number_of_Terms(p));
    dim : constant integer32
        := integer32(QuadDobl_Complex_Polynomials.Number_of_Unknowns(p));
    res : QuadDobl_Monomial_Vectors.Polynomial(dim,nbr-cnt);
    idx : integer32 := 0;

    procedure Visit_Term ( t : in QuadDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

    begin
      if Is_Zero(t.dg.all) then
        res.cff0 := t.cf;
      else
        idx := idx + 1;
        res.mons(idx) := QuadDobl_Complex_Monomials.Create(t.cf,t.dg.all);
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is
      new QuadDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    res.cff0 := QuadDobl_Complex_Numbers.Create(integer32(0));
    Visit_Terms(p);
    QuadDobl_Monomial_Vectors.Power_Update(res);
    return res;
  end Convert;

  function Convert ( p : Standard_Complex_Polynomials.Poly )
                   return Standard_Monomial_Vectors.Link_to_Polynomial is

    mpv : constant Standard_Monomial_Vectors.Polynomial := Convert(p);
    res : Standard_Monomial_Vectors.Link_to_Polynomial;

  begin
    res := new Standard_Monomial_Vectors.Polynomial'(mpv);
    return res;
  end Convert;

  function Convert ( p : DoblDobl_Complex_Polynomials.Poly )
                   return DoblDobl_Monomial_Vectors.Link_to_Polynomial is

    mpv : constant DoblDobl_Monomial_Vectors.Polynomial := Convert(p);
    res : DoblDobl_Monomial_Vectors.Link_to_Polynomial;

  begin
    res := new DoblDobl_Monomial_Vectors.Polynomial'(mpv);
    return res;
  end Convert;

  function Convert ( p : QuadDobl_Complex_Polynomials.Poly )
                   return QuadDobl_Monomial_Vectors.Link_to_Polynomial is

    mpv : constant QuadDobl_Monomial_Vectors.Polynomial := Convert(p);
    res : QuadDobl_Monomial_Vectors.Link_to_Polynomial;

  begin
    res := new QuadDobl_Monomial_Vectors.Polynomial'(mpv);
    return res;
  end Convert;

  function Convert ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                   return Standard_Polynomial_Vectors.System is

    p1 : constant Standard_Complex_Polynomials.Poly := p(p'first);
    dim : constant integer32
        := integer32(Standard_Complex_Polynomials.Number_of_Unknowns(p1));
    nbr : constant integer32 := p'last;
    res : Standard_Polynomial_Vectors.System(dim,nbr);

  begin
    for i in p'range loop
      res.pols(i) := Convert(p(i));
    end loop;
    Standard_Polynomial_Vectors.Power_Update(res);
    return res;
  end Convert;

  function Convert ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                   return DoblDobl_Polynomial_Vectors.System is

    p1 : constant DoblDobl_Complex_Polynomials.Poly := p(p'first);
    dim : constant integer32
        := integer32(DoblDobl_Complex_Polynomials.Number_of_Unknowns(p1));
    nbr : constant integer32 := p'last;
    res : DoblDobl_Polynomial_Vectors.System(dim,nbr);

  begin
    for i in p'range loop
      res.pols(i) := Convert(p(i));
    end loop;
    DoblDobl_Polynomial_Vectors.Power_Update(res);
    return res;
  end Convert;

  function Convert ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                   return QuadDobl_Polynomial_Vectors.System is

    p1 : constant QuadDobl_Complex_Polynomials.Poly := p(p'first);
    dim : constant integer32
        := integer32(QuadDobl_Complex_Polynomials.Number_of_Unknowns(p1));
    nbr : constant integer32 := p'last;
    res : QuadDobl_Polynomial_Vectors.System(dim,nbr);

  begin
    for i in p'range loop
      res.pols(i) := Convert(p(i));
    end loop;
    QuadDobl_Polynomial_Vectors.Power_Update(res);
    return res;
  end Convert;

  function Convert ( p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys )
                   return Standard_Polynomial_Vectors.Link_to_System is

    res : Standard_Polynomial_Vectors.Link_to_System:= null;

    use Standard_Complex_Poly_Systems;

  begin
    if p /= null then
      declare
        s : constant Standard_Polynomial_Vectors.System := Convert(p.all);
      begin
        res := new Standard_Polynomial_Vectors.System'(s);
      end;
    end if;
    return res;
  end Convert;

  function Convert ( p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys )
                   return DoblDobl_Polynomial_Vectors.Link_to_System is

    res : DoblDobl_Polynomial_Vectors.Link_to_System := null;

    use DoblDobl_Complex_Poly_Systems;

  begin
    if p /= null then
      declare
        s : constant DoblDobl_Polynomial_Vectors.System := Convert(p.all);
      begin
        res := new DoblDobl_Polynomial_Vectors.System'(s);
      end;
    end if;
    return res;
  end Convert;

  function Convert ( p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys )
                   return QuadDobl_Polynomial_Vectors.Link_to_System is

    res : QuadDobl_Polynomial_Vectors.Link_to_System := null;

    use QuadDobl_Complex_Poly_Systems;

  begin
    if p /= null then
      declare
        s : constant QuadDobl_Polynomial_Vectors.System := Convert(p.all);
      begin
        res := new QuadDobl_Polynomial_Vectors.System'(s);
      end;
    end if;
    return res;
  end Convert;

end System_Vector_Convertors;
