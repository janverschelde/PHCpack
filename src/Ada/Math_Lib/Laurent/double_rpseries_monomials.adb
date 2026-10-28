with unchecked_deallocation;

package body Double_rpSeries_Monomials is

-- CONSTRUCTORS :

  function Make ( cff : Double_Real_Power_Series.Series;
                  deg : Standard_Integer_Vectors.Vector )
                return Monomial is

    res : Monomial(deg'last,cff.tdx);

  begin
    res.cff := cff;
    res.deg := deg;
    return res;
  end Make;

  function Make ( cff : Double_Real_Power_Series.Series;
                  deg : Standard_Integer_Vectors.Vector )
                return Link_to_Monomial is

    res_rep : constant Monomial(deg'last,cff.tdx) := Make(cff,deg);
    res : constant Link_to_Monomial := new Monomial'(res_rep);

  begin
    return res;
  end Make;

  function Make ( x : complex_number; dim,tdx : integer32 )
                return Monomial is

    res : Monomial(dim,tdx);

  begin
    res.cff := Double_Real_Power_Series.make(x,tdx);
    res.deg := (1..dim => 0);
    return res;
  end Make;

  function Make ( x : complex_number; dim,tdx : integer32 )
                return Link_to_Monomial is

    res_rep : constant Monomial(dim,tdx) := Make(x,dim,tdx);
    res : constant Link_to_Monomial := new Monomial'(res_rep);

  begin
    return res;
  end Make;

-- SELECTORS :

  function Eval ( mon : Monomial; x : Double_rpSeries_Vectors.Vector )
                return Double_Real_Power_Series.Link_to_Series is

    res : Double_Real_Power_Series.Link_to_Series;

  begin
    res := Double_Real_Power_Series.make(mon.cff.cff,mon.cff.pwt);
    for i in mon.deg'range loop
      if mon.deg(i) /= 0 then
        declare
          use Double_Real_Power_Series;
          pwr : Link_to_Series := x(i)**mon.deg(i);
        begin
          Double_Real_Power_Series.mul(res,pwr);
          Double_Real_Power_Series.clear(pwr);
        end;
      end if;
    end loop;
    return res;
  end Eval;

  function Eval ( mon : Link_to_Monomial; x : Double_rpSeries_Vectors.Vector )
                return Double_Real_Power_Series.Link_to_Series is

    res : Double_Real_Power_Series.Link_to_Series;

  begin
    if mon = null
     then res := Double_Real_Power_Series.Make(integer32(0),x(x'first).tdx);
     else res := Eval(mon.all,x);
    end if;
    return res;
  end Eval;

-- DESTRUCTOR :

  procedure Clear ( mon : in out Link_to_Monomial ) is

    procedure free is new unchecked_deallocation(Monomial,Link_to_Monomial);

  begin
    if mon /= null
     then free(mon);
    end if;
  end Clear;

end Double_rpSeries_Monomials;
