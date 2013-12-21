with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Set_Structure;

package body Supporting_Set_Structure is

  function Is_Supporting
             ( p : Poly; i : natural32; verbose : boolean ) return boolean is

    res : boolean := true;
    n : constant natural32 := Number_of_Unknowns(p);
    deg,cnt : integer32;

  begin
    for k in 1..n loop
      deg := Degree(p,integer32(k));
      if verbose then
        put("Degree of p in variable "); put(k,1);
        put(" : "); put(deg,1);
      end if;
      cnt := 0;
      for j in 1..Set_Structure.Number_of_Sets(i) loop
        if Set_Structure.Is_In(i,j,k)
         then cnt := cnt + 1;
        end if;
      end loop;
      if verbose
       then put(" count : "); put(cnt,1); new_line;
      end if;
      if res then
        if cnt < deg
         then res := false;
        end if;
      end if;
    end loop;
    return res;
  end Is_Supporting;

  function Is_Supporting ( p : Poly_Sys; verbose : boolean ) return boolean is

    res : boolean := true;
    sup : boolean;

  begin
    if Set_Structure.Empty then
      return false;
    else
      if natural32(p'last) /= Set_Structure.Dimension then
        return false;
      else
        for i in p'range loop
          sup := Is_Supporting(p(i),natural32(i),verbose);
          if res
           then res := sup;
          end if;
        end loop;
        return res;
      end if;
    end if;
  end Is_Supporting;

end Supporting_Set_Structure;
