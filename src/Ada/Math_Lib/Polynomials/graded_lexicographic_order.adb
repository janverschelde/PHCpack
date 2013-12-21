with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;

package body Graded_Lexicographic_Order is

-- FOR STANDARD NATURAL VECTORS :

  function "<" ( v1,v2 : Standard_Natural_Vectors.Vector ) return boolean is

    use Standard_Natural_Vectors;
    s1,s2 : natural32;

  begin
    s1 := Sum(v1);
    s2 := Sum(v2);
    if s1 < s2 then
      return true;
    elsif s1 > s2 then
      return false;
    elsif v1'first /= v2'first or else v1'last /= v2'last then
      raise CONSTRAINT_ERROR;
    else
      for i in v1'range loop
        if v1(i) < v2(i) then
          return true;
        elsif v1(i) > v2(i) then
          return false;
        end if;
      end loop;
      return false;  -- v1 = v2
    end if;
  end "<";

  function "<" ( v1,v2 : Standard_Natural_Vectors.Link_to_Vector )
               return boolean is

    use Standard_Natural_Vectors;

  begin
    if v2 = null then
      return false;
    elsif v1 = null then
      if Sum(v2) > 0
       then return true;
       else return false;
      end if;
    else
      return v1.all < v2.all;
    end if;
  end "<";

  function ">" ( v1,v2 : Standard_Natural_Vectors.Vector ) return boolean is

    use Standard_Natural_Vectors;

    s1,s2 : natural32;

  begin
    s1 := Sum(v1);
    s2 := Sum(v2);
    if s1 < s2 then
      return false;
    elsif s1 > s2 then
      return true;
    elsif v1'first /= v2'first or else v1'last /= v2'last then
      raise CONSTRAINT_ERROR;
    else
      for i in v1'range loop
        if v1(i) < v2(i) then
          return false;
        elsif v1(i) > v2(i) then
          return true;
        end if;
      end loop;
      return false;  -- v1 = v2
    end if;
  end ">";

  function ">" ( v1,v2 : Standard_Natural_Vectors.Link_to_Vector )
               return boolean is

    use Standard_Natural_Vectors;

  begin
    if v1 = null then
      return false;
    elsif v2 = null then
      if Sum(v1) > 0
       then return true;
       else return false;
      end if;
    else
      return v1.all > v2.all;
    end if;
  end ">";

-- FOR STANDARD INTEGER VECTORS :

  function "<" ( v1,v2 : Standard_Integer_Vectors.Vector ) return boolean is

    use Standard_Integer_Vectors;

    s1,s2 : integer32;

  begin
    s1 := Sum(v1);
    s2 := Sum(v2);
    if s1 < s2 then
      return true;
    elsif s1 > s2 then
      return false;
    elsif v1'first /= v2'first or else v1'last /= v2'last then
      raise CONSTRAINT_ERROR;
    else
      for i in v1'range loop
        if v1(i) < v2(i) then
          return true;
        elsif v1(i) > v2(i) then
          return false;
        end if;
      end loop;
      return false;  -- v1 = v2
    end if;
  end "<";

  function "<" ( v1,v2 : Standard_Integer_Vectors.Link_to_Vector )
               return boolean is

    use Standard_Integer_Vectors;

  begin
    if v2 = null then
      return false;
    elsif v1 = null then
      if Sum(v2) > 0
       then return true;
       else return false;
      end if;
    else
      return v1.all < v2.all;
    end if;
  end "<";

  function ">" ( v1,v2 : Standard_Integer_Vectors.Vector ) return boolean is

    use Standard_Integer_Vectors;

    s1,s2 : integer32;

  begin
    s1 := Sum(v1);
    s2 := Sum(v2);
    if s1 < s2 then
      return false;
    elsif s1 > s2 then
      return true;
    elsif v1'first /= v2'first or else v1'last /= v2'last then
      raise CONSTRAINT_ERROR;
    else
      for i in v1'range loop
        if v1(i) < v2(i) then
          return false;
        elsif v1(i) > v2(i) then
          return true;
        end if;
      end loop;
      return false;  -- v1 = v2
    end if;
  end ">";

  function ">" ( v1,v2 : Standard_Integer_Vectors.Link_to_Vector )
               return boolean is

    use Standard_Integer_Vectors;

  begin
    if v1 = null then
      return false;
    elsif v2 = null then
      if Sum(v1) > 0
       then return true;
       else return false;
      end if;
    else
      return v1.all > v2.all;
    end if;
  end ">";

end Graded_Lexicographic_Order;
