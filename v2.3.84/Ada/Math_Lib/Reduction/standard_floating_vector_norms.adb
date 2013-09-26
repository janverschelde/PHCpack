with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;

package body Standard_Floating_Vector_Norms is

  function Norm2 ( v : Vector ) return double_float is

    res : double_float := 0.0;

  begin
    for i in v'range loop
       res := res + v(i)*v(i);
    end loop;
    res := SQRT(res);
    return res;
  end Norm2;

  function Sum_Norm ( v : Vector ) return double_float is

    res : double_float := 0.0;

  begin
    for i in v'range loop
      res := res + abs(v(i));
    end loop;
    return res;
  end Sum_Norm;

  function Max_Norm ( v : Vector ) return double_float is

    res : double_float := abs(v(v'first));
    tmp : double_float;

  begin
    for i in v'first+1..v'last loop
      tmp := abs(v(i));
      if tmp > res
       then res := tmp;
      end if;
    end loop;
    return res;
  end Max_Norm;

end Standard_Floating_Vector_Norms;
