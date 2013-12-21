with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package body Generic_Norms_Equals is

  function Max_Norm ( v : Vector ) return number is

    res : number := AbsVal(v(v'first));

  begin
    for i in v'first+1..v'last loop
      declare
        abstmp : number := AbsVal(v(i));
      begin
        if abstmp > res
         then Copy(abstmp,res);
        end if;
        Clear(abstmp);
      end;
    end loop;
    return res;
  end Max_Norm;

  function Sum_Norm ( v : Vector ) return number is

    res : number := AbsVal(v(v'first));

  begin
    for i in v'first+1..v'last loop
      declare
        abstmp : number := AbsVal(v(i));
      begin
        Add(res,abstmp);
        Clear(abstmp);
      end;
    end loop;
    return res;
  end Sum_Norm;

  function Max_Norm ( m : Matrix ) return number is

    res : number := Create(0);
 
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        declare
          abstmp : number := AbsVal(m(i,j));
        begin
          if abstmp > res
           then Copy(abstmp,res);
          end if;
          Clear(abstmp);
        end;
      end loop;
    end loop;
    return res;
  end Max_Norm;

  function Equal ( x,y,tol : number ) return boolean is

    dif : number := x-y;
    absdif : number := AbsVal(dif);
    res : constant boolean := (absdif < tol);

  begin
    Clear(dif);
    Clear(absdif);
    return res;
  end Equal;

  function Equal ( x,y : Vector; tol : number ) return boolean is
  begin
    for i in x'range loop
      if not Equal(x(i),y(i),tol)
       then return false;
      end if;
    end loop;
    return true;
  end Equal;

end Generic_Norms_Equals;
