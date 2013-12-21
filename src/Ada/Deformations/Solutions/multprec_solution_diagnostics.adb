with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;

package body Multprec_Solution_Diagnostics is

  function Is_Real ( sol : Solution; tol : Floating_Number ) return boolean is

    res : boolean := true;

  begin
    for i in sol.v'range loop
      declare
        imapar : Floating_Number := IMAG_PART(sol.v(i));
        abstmp : Floating_Number := AbsVal(imapar);
      begin
        if abstmp > tol
         then res := false;
        end if;
        Clear(imapar); Clear(abstmp);
      end;
      exit when not res;
    end loop;
    return res;
  end Is_Real;

  function Equal ( s1,s2 : Solution; tol : Floating_Number ) return boolean is
  begin
    for i in s1.v'range loop
      if not Equal(s1.v(i),s2.v(i),tol)
       then return false;
      end if;
    end loop;
    return true;
  end Equal;

  function Is_Clustered ( sol : Solution; nb : natural32;
                          sols : Solution_List; tol : Floating_Number )
                        return natural32 is

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

  function Is_Clustered ( sol : Solution; nb : natural32;
                          sols : Solution_Array; tol : Floating_Number )
                        return natural32 is
  begin
    for i in sols'range loop
      if natural32(i) /= nb then
        if Equal(sol,sols(i).all,tol)
         then return natural32(i);
        end if;
      end if;
    end loop;
    return nb;
  end Is_Clustered;

  function Multiplicity ( sol : Solution; sols : Solution_List; 
                          tol : Floating_Number ) return natural32 is

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
                          tol : Floating_Number ) return natural32 is

    cnt : natural32 := 0;

  begin
    for i in sols'range loop
      if Equal(sol,sols(i).all,tol)
       then cnt := cnt + 1;
      end if;
    end loop;
    return cnt;
  end Multiplicity;

end Multprec_Solution_Diagnostics;
