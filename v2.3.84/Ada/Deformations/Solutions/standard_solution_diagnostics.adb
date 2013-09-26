with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;

package body Standard_Solution_Diagnostics is

  function Is_Real ( sol : Solution; tol : double_float ) return boolean is
  begin
    for i in sol.v'range loop
      if ABS(IMAG_PART(sol.v(i))) > tol
       then return false;
      end if;
    end loop;
    return true;
  end Is_Real;

  function Equal ( s1,s2 : Solution; tol : double_float ) return boolean is
  begin
    for i in s1.v'range loop
      if not Equal(s1.v(i),s2.v(i),tol)
       then return false;
      end if;
    end loop;
    return true;
  end Equal;

  function Is_Clustered ( sol : Solution; nb : natural32; sols : Solution_List;
                          tol : double_float ) return natural32 is

    tmp : Solution_List := sols;
    cnt : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      if cnt /= nb then
        if Equal(sol,Head_Of(tmp).all,tol)
         then return cnt;
        end if;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return nb;
  end Is_Clustered;

  function Is_Clustered
             ( sol : Solution; nb : natural32; sols : Solution_Array;
               tol : double_float ) return natural32 is
  begin
    for i in sols'range loop
      if i /= integer32(nb) then
        if Equal(sol,sols(i).all,tol)
         then return natural32(i);
        end if;
      end if;
    end loop;
    return nb;
  end Is_Clustered;

  function Multiplicity ( sol : Solution; sols : Solution_List; 
                          tol : double_float ) return natural32 is

    tmp : Solution_List := sols;
    cnt : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      if Equal(sol,Head_Of(tmp).all,tol)
       then cnt := cnt + 1;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return cnt;
  end Multiplicity;

  function Multiplicity ( sol : Solution; sols : Solution_Array;
                          tol : double_float ) return natural32 is

    cnt : natural32 := 0;

  begin
    for i in sols'range loop
      if Equal(sol,sols(i).all,tol)
       then cnt := cnt + 1;
      end if;
    end loop;
    return cnt;
  end Multiplicity;

  function At_Infinity ( sol : Solution; prj : boolean;
                         tol : double_float ) return boolean is
  begin
    if prj then
      if AbsVal(sol.v(sol.v'last)) < 1.0/tol
       then return true;
       else return false;
      end if;
    else
      for i in 1..sol.n loop
        if AbsVal(sol.v(i)) > tol
         then return true;
        end if;
      end loop;
      return false;
    end if;
  end At_Infinity;

end Standard_Solution_Diagnostics;
