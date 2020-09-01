with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with PentDobl_Mathematical_Functions;    use PentDobl_Mathematical_Functions;

package body Penta_Double_Vector_Norms is

  function Norm2 ( v : Vector ) return penta_double is

    res : penta_double := Create(0.0);

  begin
    for i in v'range loop
       res := res + v(i)*v(i);
    end loop;
    res := SQRT(res);
    return res;
  end Norm2;

  function Sum_Norm ( v : Vector ) return penta_double is

    res : penta_double := Create(0.0);

  begin
    for i in v'range loop
      res := res + AbsVal(v(i));
    end loop;
    return res;
  end Sum_Norm;

  function Max_Norm ( v : Vector ) return penta_double is

    res : penta_double := AbsVal(v(v'first));
    tmp : penta_double;

  begin
    for i in v'first+1..v'last loop
      tmp := AbsVal(v(i));
      if tmp > res
       then res := tmp;
      end if;
    end loop;
    return res;
  end Max_Norm;

end Penta_Double_Vector_Norms;
