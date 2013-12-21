with Standard_Common_Divisors;           use Standard_Common_Divisors;
with Standard64_Common_Divisors;         use Standard64_Common_Divisors;

package body Standard_Integer_Norms is

  function gcd ( v : Standard_Integer_Vectors.Vector ) return integer32 is

    tmp : Standard_Integer_Vectors.Vector(v'range);
    res : integer32;

  begin
    for i in v'range loop
      if v(i) < 0
       then tmp(i) := -v(i);
       else tmp(i) :=  v(i);
      end if;
    end loop;
    res := tmp(tmp'first);
    for i in (tmp'first+1)..tmp'last loop
      res := gcd(res,tmp(i));
      exit when (res = 1);
    end loop;
    return res;
  end gcd;

  function gcd ( v : Standard_Integer64_Vectors.Vector ) return integer64 is

    tmp : Standard_Integer64_Vectors.Vector(v'range);
    res : integer64;

  begin
    for i in v'range loop
      if v(i) < 0
       then tmp(i) := -v(i);
       else tmp(i) :=  v(i);
      end if;
    end loop;
    res := tmp(tmp'first);
    for i in (tmp'first+1)..tmp'last loop
      res := gcd(res,tmp(i));
      exit when (res = 1);
    end loop;
    return res;
  end gcd;

  procedure Normalize ( v : in out Standard_Integer_Vectors.Vector ) is

    g : constant integer32 := gcd(v);

  begin
    if (g /= 0) and then (g /= 1) then
      for i in v'range loop
        v(i) := v(i)/g;
      end loop;
    end if;
  end Normalize;

  procedure Normalize ( v : in out Standard_Integer64_Vectors.Vector ) is

    g : constant integer64 := gcd(v);

  begin
    if (g /= 0) and then (g /= 1) then
      for i in v'range loop
        v(i) := v(i)/g;
      end loop;
    end if;
  end Normalize;

end Standard_Integer_Norms;
