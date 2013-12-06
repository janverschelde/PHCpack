with Standard_Integer_Numbers;            use Standard_Integer_Numbers;

package body Permutations is

  function Is_Permutation ( p : Permutation ) return boolean is
  begin
    for i in p'range loop
      if ((p(i) = 0) or else (p(i) < -p'last) or else (p(i) > p'last)) then
        return false;
      else
        for j in p'first..(i-1) loop
          if ((p(i) = p(j)) or else (p(i) = -p(j)))
           then return false;
          end if;
        end loop;
      end if;
    end loop;
    return true;
  end Is_Permutation;

  function Equal ( p1,p2 : Permutation ) return boolean is
  begin
    if (p1'first /= p2'first) or else (p1'last /= p2'last) then
      return false;
    else
      for i in p1'range loop
        if p1(i) /= p2(i)
         then return false;
        end if;
      end loop;
      return true;
    end if;
  end Equal;

  function "*" ( p1,p2 : Permutation ) return Permutation is

    r : Permutation(p1'range);

  begin
    for i in r'range loop
      if p2(i) >= 0
       then r(i) := p1(p2(i));
       else r(i) := -p1(-p2(i));
      end if;
    end loop;
    return r;
  end "*";

  function inv ( p : Permutation ) return Permutation is

    r : Permutation(p'range);

  begin
    for i in r'range loop
      if p(i) >= 0
       then r(p(i)) := i;
       else r(-p(i)) := -i;
      end if;
    end loop;
    return r;
  end inv;

end Permutations;
