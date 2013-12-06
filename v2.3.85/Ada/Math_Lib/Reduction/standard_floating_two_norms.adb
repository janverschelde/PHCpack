with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;

package body Standard_Floating_Two_Norms is

  function Norm2 ( v : Vector ) return double_float is

    res : double_float := 0.0;

  begin
    for i in v'range loop
      res := res + v(i)*v(i);
    end loop;
    res := SQRT(res);
    return res;
  end Norm2;

  procedure Normalize ( v : in out Vector ) is

    nrm : constant double_float := Norm2(v);

  begin
    for i in v'range loop
      v(i) := v(i)/nrm;
    end loop;
  end Normalize;

  procedure Normalize( v : in out VecVec ) is
  begin
    for i in v'range loop
      Normalize(v(i).all);
    end loop;
  end Normalize;

end Standard_Floating_Two_Norms;
