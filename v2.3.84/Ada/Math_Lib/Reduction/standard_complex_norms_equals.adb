with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;

package body Standard_Complex_Norms_Equals is

  function Conjugated_Inner_Product ( v,w : Vector ) return Complex_Number is

    res : Complex_Number := Create(0.0);

  begin
    for i in v'range loop
      res := res + Conjugate(v(i))*w(i);
    end loop;
    return res;
  end Conjugated_Inner_Product;

  function Norm2 ( v : Vector ) return double_float is

    res : double_float := 0.0;
  
  begin
    for i in v'range loop
      res := res + REAL_PART(v(i))*REAL_PART(v(i))
                 + IMAG_PART(v(i))*IMAG_PART(v(i));
    end loop;
    return SQRT(res);
  end Norm2;

  procedure Normalize ( v : in out Vector ) is

    nrm : constant double_float := Norm2(v);
    d : Complex_Number; 

  begin
    if nrm + 1.0 /= 1.0 then
      d := Create(nrm);
      for i in v'range loop
        v(i) := v(i)/d;
      end loop;
    end if;
  end Normalize;

  function Max_Norm ( v : Vector ) return double_float is

    res : double_float := AbsVal(v(v'first));

  begin
    for i in v'first+1..v'last loop
      declare
        abstmp : constant double_float := AbsVal(v(i));
      begin
        if abstmp > res
         then res := abstmp;
        end if;
      end;
    end loop;
    return res;
  end Max_Norm;

  function Sum_Norm ( v : Vector ) return double_float is

    res : double_float := AbsVal(v(v'first));

  begin
    for i in v'first+1..v'last loop
      res := res + AbsVal(v(i));
    end loop;
    return res;
  end Sum_Norm;

  function Max_Norm ( m : Matrix ) return double_float is

    res : double_float := 0.0;

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        declare
          abstmp : constant double_float := AbsVal(m(i,j));
        begin
          if abstmp > res
           then res := abstmp;
          end if;
        end;
      end loop;
    end loop;
    return res;
  end Max_Norm;

  function Equal ( x,y : Complex_Number; tol : double_float ) return boolean is

    dif : constant Complex_Number := x-y;
    absdif : constant double_float := AbsVal(dif);
    res : constant boolean := (absdif < tol);

  begin
    return res;
  end Equal;

  function Equal ( x,y : Vector; tol : double_float ) return boolean is
  begin
    for i in x'range loop
      if not Equal(x(i),y(i),tol)
       then return false;
      end if;
    end loop;
    return true;
  end Equal;

end Standard_Complex_Norms_Equals;
