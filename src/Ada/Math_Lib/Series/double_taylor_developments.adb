with Binomial_Coefficients;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;

package body Double_Taylor_Developments is

  function Double_Taylor_Coefficients 
             ( deg : integer32; alpha, point : double_float )
             return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(0..deg) := (0..deg => 0.0);
    cff : double_float := alpha;

  begin
    res(0) := point**alpha;
    res(1) := alpha*(point**(alpha-1.0));
    for k in 2..deg loop
      cff := cff*(alpha-1.0)/double_float(k);
      res(k) := cff*(point**(alpha-double_float(k)));
    end loop;
    return res;
  end Double_Taylor_Coefficients;

  function Double_Taylor_Value
             ( cff : Standard_Floating_Vectors.Vector;
               arg, point : double_float ) return double_float is

    res : double_float := cff(0);
    inc : double_float := 1.0;
    argminpoint : constant double_float := arg - point;

  begin
    for k in 1..cff'last loop
      inc := inc*argminpoint;
      res := res + cff(k)*inc;
    end loop;
    return res;
  end Double_Taylor_Value;

  function Double_Taylor_Value
             ( cff : Standard_Floating_Vectors.Vector;
               arg : double_float ) return double_float is

    res : double_float := cff(0);
    inc : double_float := 1.0;

  begin
    for k in 1..cff'last loop
      inc := inc*arg;
      res := res + cff(k)*inc;
    end loop;
    return res;
  end Double_Taylor_Value;

  function Double_Taylor_Expansion
             ( cff : Standard_Floating_Vectors.Vector;
               point : double_float )
             return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(cff'range) := (cff'range => 0.0);
    bincff : double_float;
    expcff : natural;

  begin
    for k in cff'range loop
      for ell in 0..k loop
        expcff := natural(k - ell);
        bincff := Binomial_Coefficients.binomial(k,ell)*(point**expcff);
        if expcff mod 2 = 0
         then res(ell) := res(ell) + cff(k)*bincff;
         else res(ell) := res(ell) - cff(k)*bincff;
        end if;
      end loop;
    end loop;
    return res;
  end Double_Taylor_Expansion;

end Double_Taylor_Developments;
