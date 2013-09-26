with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Common_Divisors;           use Multprec_Common_Divisors;

package body Multprec_Integer_Norms is

  function gcd ( v : Multprec_Integer_Vectors.Vector )
               return Integer_Number is

    tmp : Multprec_Integer_Vectors.Vector(v'range);
    res,gcd_res : Integer_Number;

  begin
    for i in v'range loop
      if v(i) < 0
       then tmp(i) := -v(i);
       else Copy(v(i),tmp(i));
      end if;
    end loop;
    Copy(tmp(tmp'first),res);
    for i in (tmp'first+1)..tmp'last loop
      gcd_res := gcd(res,tmp(i));
      Copy(gcd_res,res); Clear(gcd_res);
      exit when (Equal(res,1));
    end loop;
    for i in tmp'range loop
      Clear(tmp(i));
    end loop;
    return res;
  end gcd;

  procedure Normalize ( v : in out Multprec_Integer_Vectors.Vector ) is

    g : Integer_Number := gcd(v);

  begin
    if (not Equal(g,0)) and then (not Equal(g,1)) then
      for i in v'range loop
        Div(v(i),g);
      end loop;
    end if;
    Clear(g);
  end Normalize;

end Multprec_Integer_Norms;
