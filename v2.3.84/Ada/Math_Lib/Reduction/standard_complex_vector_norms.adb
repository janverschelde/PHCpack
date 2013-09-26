with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package body Standard_Complex_Vector_Norms is

  function Norm2 ( v : Vector ) return double_float is

    res : double_float := 0.0;

  begin
    for i in v'range loop
      res := res + REAL_PART(v(i))*REAL_PART(v(i))
                 + IMAG_PART(v(i))*IMAG_PART(v(i));
    end loop;
    res := SQRT(res);
    return res;
  end Norm2;

  function Sum_Norm ( v : Vector ) return double_float is

    res : double_float := 0.0;

  begin
    for i in v'range loop
      res := res + AbsVal(v(i));
    end loop;
    return res;
  end Sum_Norm;

  function Max_Norm ( v : Vector ) return double_float is

    res : double_float := AbsVal(v(v'first));
    tmp : double_float;

  begin
    for i in v'first+1..v'last loop
      tmp := AbsVal(v(i));
      if tmp > res
       then res := tmp;
      end if;
    end loop;
    return res;
  end Max_Norm;

end Standard_Complex_Vector_Norms;
