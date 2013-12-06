with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Mathematical_Functions;    use Multprec_Mathematical_Functions;

package body Multprec_Complex_Norms_Equals is

  function Norm2 ( v : Vector ) return Floating_Number is

    res : Floating_Number := Create(integer(0));
    wrk : Floating_Number;
  
  begin
    for i in v'range loop
      wrk := REAL_PART(v(i));
      Mul(wrk,wrk);
      Add(res,wrk);
      Clear(wrk);
      wrk := IMAG_PART(v(i));
      Mul(wrk,wrk);
      Add(res,wrk);
      Clear(wrk);
    end loop;
    wrk := SQRT(res);
    Clear(res);
    return wrk;
  end Norm2;

  function Max_Norm ( v : Vector ) return Floating_Number is

    res : Floating_Number := AbsVal(v(v'first));

  begin
    for i in v'first+1..v'last loop
      declare
        abstmp : Floating_Number := AbsVal(v(i));
      begin
        if abstmp > res
         then Copy(abstmp,res);
        end if;
        Clear(abstmp);
      end;
    end loop;
    return res;
  end Max_Norm;

  function Sum_Norm ( v : Vector ) return Floating_Number is

    res : Floating_Number := AbsVal(v(v'first));

  begin
    for i in v'first+1..v'last loop
      declare
        abstmp : Floating_Number := AbsVal(v(i));
      begin
        Add(res,abstmp);
        Clear(abstmp);
      end;
    end loop;
    return res;
  end Sum_Norm;

  function Equal ( x,y : Complex_Number;
                   tol : Floating_Number ) return boolean is

    dif : Complex_Number := x-y;
    absdif : Floating_Number := AbsVal(dif);
    res : constant boolean := (absdif < tol);

  begin
    Clear(dif);
    Clear(absdif);
    return res;
  end Equal;

  function Equal ( x,y : Vector; tol : Floating_Number ) return boolean is
  begin
    for i in x'range loop
      if not Equal(x(i),y(i),tol)
       then return false;
      end if;
    end loop;
    return true;
  end Equal;

end Multprec_Complex_Norms_Equals;
