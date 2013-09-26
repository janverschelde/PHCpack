with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;

procedure ts_subsets is

-- DESCRIPTION :
--   Generates all subsets of k elements from a set of n elements.

  function Complement ( n : integer32; v : Vector ) return Vector is

  -- DESCRIPTION :
  --   Returns the complement of the vector w.r.t. the set 1..n.

    res : Vector(1..n-v'length);
    cnt : integer32 := 0;
    found : boolean;

  begin
    for i in 1..n loop
      found := false;
      for j in v'range loop
        if v(j) = i
         then found := true; exit;
        end if;
      end loop;
      if not found  then
        cnt := cnt + 1;
        res(cnt) := i;
      end if;
    end loop;
    return res;
  end Complement;

  procedure Enumerate ( start,i,n : in integer32; accu : in out Vector ) is

  -- DESCRIPTION :
  --   Enumerates all subsets of 1..n, of size accu'length, starting to
  --   fill up accu(i) with entries in start..n.

  begin
    if i > accu'last then
      put("Subset : "); put(accu);
      put("  Complement : "); put(Complement(n,accu)); new_line;
    else
      for l in start..n loop
        accu(i) := l;
        Enumerate(l+1,i+1,n,accu);
      end loop;
    end if;
  end Enumerate;

  procedure Main is

    k,n : integer32 := 0;

  begin
    put("Give the cardinality of whole set : "); get(n);
    put("Give the cardinality of subset : "); get(k);
    declare
      acc : Vector(1..k);
    begin
      Enumerate(1,1,n,acc);
    end;
  end Main;

begin
  Main;
end ts_subsets;
