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

-- DESTRUCTOR :

  procedure Clear ( mon : in out Link_to_Monomial ) is

    procedure free is new unchecked_deallocation(Monomial,Link_to_Monomial);

  begin
    if mon /= null
     then free(mon);
    end if;
  end Clear;

end Double_rpSeries_Monomials;
