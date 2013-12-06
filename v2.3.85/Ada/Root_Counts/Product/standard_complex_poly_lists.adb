package body Standard_Complex_Poly_Lists is

  function Create ( p : Poly ) return Prod_Poly is

    lp : Lists_of_Standard_Complex_Polynomials.List;
    res : Prod_Poly;

  begin
    Lists_of_Standard_Complex_Polynomials.Construct(p,lp);
    res := Prod_Poly(lp);
    return res;
  end Create;

  function Number_of_Factors ( p : Prod_Poly ) return natural32 is
  begin
    return Lists_of_Standard_Complex_Polynomials.Length_Of
             (Lists_of_Standard_Complex_Polynomials.List(p));
  end Number_of_Factors;

  function Number_of_Unknowns ( p : Prod_Poly ) return natural32 is

    use Lists_of_Standard_Complex_Polynomials;

  begin
    if Is_Null(Lists_of_Standard_Complex_Polynomials.List(p))
     then return 0;
     else return Number_of_Unknowns(Head_Of(p));
    end if;
  end Number_of_Unknowns;

  function Expand ( p : Prod_Poly ) return Poly is

    res : Poly := Null_Poly;
    tmp : Prod_Poly;

    use Lists_of_Standard_Complex_Polynomials;

  begin
    if not Is_Null(p) then
      Copy(Head_Of(p),res);
      tmp := Tail_Of(p);
      while not Is_Null(tmp) loop
        Mul(res,Head_Of(tmp));
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    return res;
  end Expand;

  procedure Shallow_Clear ( p : in out Prod_Poly ) is
  begin
    Lists_of_Standard_Complex_Polynomials.Clear
      (Lists_of_Standard_Complex_Polynomials.List(p));
  end Shallow_Clear;

  procedure Deep_Clear ( p : in out Prod_Poly ) is

    tmp : Prod_Poly := p;
    i : Poly;

    use Lists_of_Standard_Complex_Polynomials;

  begin
    while not Is_Null(tmp) loop
      i := Head_Of(tmp);
      Clear(i);
      tmp := Tail_Of(tmp);
    end loop;
    Shallow_Clear(p);
  end Deep_Clear;

end Standard_Complex_Poly_Lists;
