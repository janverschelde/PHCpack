with unchecked_deallocation;

package body Generic_Vectors is

-- COMPARISON AND COPYING :

  function Equal ( v1,v2 : Vector ) return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      return false;
    else
      for i in v1'range loop
        if not equal(v1(i),v2(i))
         then return false;
        end if;
      end loop;
      return true;
    end if;
  end Equal;

  procedure Copy ( v1: in Vector; v2 : in out Vector ) is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      raise CONSTRAINT_ERROR;
    else
      Clear(v2);
      for i in v1'range loop
        copy(v1(i),v2(i));
      end loop;
    end if;
  end Copy;

-- ARITHMETIC AS FUNCTIONS :

  function "+" ( v1,v2 : Vector ) return Vector is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      raise CONSTRAINT_ERROR;
    else
      declare
        res : Vector(v1'range);
      begin
        for i in v1'range loop
          res(i) := v1(i) + v2(i);
        end loop;
        return res;
      end;
    end if;
  end "+";

  function "+" ( v : Vector ) return Vector is

    res : Vector(v'range);

  begin
    for i in v'range loop
      res(i) := +v(i);
    end loop;
    return res;
  end "+";
  
  function "-" ( v : Vector ) return Vector is

    res : Vector(v'range);

  begin
    for i in v'range loop
      res(i) := -v(i);
    end loop;
    return res;
  end "-";

  function "-" ( v1,v2 : Vector ) return Vector is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      raise CONSTRAINT_ERROR;
    else
      declare
        res : Vector(v1'range);
      begin
        for i in v1'range loop
          res(i) := v1(i) - v2(i);
        end loop;
        return res;
      end;
    end if;
  end "-";

  function "*" ( v : Vector; a : number ) return Vector is

    res : Vector(v'range);

  begin
    for i in v'range loop
      res(i) := v(i) * a;
    end loop;
    return res;
  end "*";

  function "*" ( a : number; v : Vector ) return Vector is
  begin
    return v*a;
  end "*";

  function "*" ( v1,v2 : Vector ) return number is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      raise CONSTRAINT_ERROR;
    else
      declare 
        temp,sum : number;
      begin
        if v1'first <= v1'last then
          sum := v1(v1'first)*v2(v2'first);
          for i in v1'first+1..v1'last loop
            temp := v1(i)*v2(i);
            Add(sum,temp);
            Clear(temp);
          end loop;
        end if;
        return sum;
      end;
    end if;
  end "*";

  function "*" ( v1,v2 : Vector ) return Vector is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      raise CONSTRAINT_ERROR;
    else
      declare
        res : Vector(v1'range);
      begin
        for i in v1'range loop
          res(i) := v1(i)*v2(i);
        end loop;
        return res;
      end;      
    end if;
  end "*";

  function Sum ( v : Vector ) return number is

    res : number;

  begin
    if v'first > v'last then
      return Ring.zero;
    else
      Copy(v(v'first),res);
      for i in v'first+1..v'last loop
        Add(res,v(i));
      end loop;
      return res;
    end if;
  end Sum;

-- ARITHMETIC AS PROCEDURES :

  procedure Add ( v1 : in out Vector; v2 : in Vector ) is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      raise CONSTRAINT_ERROR;
    else
      for i in v1'range loop
        Add(v1(i),v2(i));
      end loop;
    end if;
  end Add;

  procedure Min ( v : in out Vector ) is
  begin
    for i in v'range loop
      Min(v(i));
    end loop;
  end Min;

  procedure Sub ( v1 : in out Vector; v2 : in Vector ) is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      raise CONSTRAINT_ERROR;
    else
      for i in v1'range loop
        Sub(v1(i),v2(i));
      end loop;
    end if;
  end Sub;

  procedure Mul ( v : in out Vector; a : in number ) is
  begin
    for i in v'range loop
      Mul(v(i),a);
    end loop;
  end Mul;

  procedure Mul ( v1 : in out Vector; v2 : in Vector ) is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last then
      raise CONSTRAINT_ERROR;
    else
      for i in v1'range loop
        Mul(v1(i),v2(i));
      end loop;
    end if;
  end Mul;

-- DESTRUCTOR :

  procedure Clear ( v : in out Vector ) is
  begin
    for i in v'range loop
      Clear(v(i));
    end loop;
  end Clear;

-- OPERATIONS ON POINTERS TO VECTORS :

-- COMPARISON AND COPYING :

  function Equal ( v1,v2 : Link_to_Vector ) return boolean is
  begin
    if (v1 = null) and (v2 = null) then
      return true;
    elsif (v1 = null) or (v2 = null) then
      return false;
    else
      return Equal(v1.all,v2.all);
    end if;
  end Equal;

  procedure Copy ( v1: in Link_to_Vector; v2 : in out Link_to_Vector ) is
  begin
    Clear(v2);
    if v1 /= null then
      v2 := new Vector(v1'range);
      for i in v1'range loop
        v2(i) := v1(i);
      end loop;
    end if;
  end Copy;

-- ARITHMETIC AS FUNCTIONS :

  function "+" ( v1,v2 : Link_to_Vector ) return Link_to_Vector is
  begin
    if v1 = null then
      return v2;
    elsif v2 = null then
      return v1;
    else
      return new Vector'(v1.all + v2.all);
    end if;
  end "+";

  function "+" ( v : Link_to_Vector ) return Link_to_Vector is
  begin
    if v = null
     then return v;
     else return new Vector'(+v.all);
    end if;
  end "+";

  function "-" ( v : Link_to_Vector ) return Link_to_Vector is
  begin
    if v = null
     then return v;
     else return new Vector'(-v.all);
    end if;
  end "-";

  function "-" ( v1,v2 : Link_to_Vector ) return Link_to_Vector is
  begin
    if v2 = null then
      return v1;
    elsif v1 = null then
      return -v2;
    else
      return new Vector'(v1.all - v2.all);
    end if;
  end "-";

  function "*" ( v : Link_to_Vector; a : number ) return Link_to_Vector is
  begin
    if v = null 
     then return null;
     else return new Vector'(v.all*a);
    end if;
  end "*";

  function "*" ( a : number; v : Link_to_Vector ) return Link_to_Vector is
  begin
    return v*a;
  end "*";

  function "*" ( v1,v2 : Link_to_Vector ) return number is
  begin
    return v1.all*v2.all;
  end "*";

  function "*" ( v1,v2 : Link_to_Vector ) return Link_to_Vector is
  begin
    if (v1 = null) or (v2 = null)
     then return null;
     else return new Vector'(v1.all*v2.all);
    end if;
  end "*";

  function Sum ( v : Link_to_Vector ) return number is
  begin
    if v = null
     then return Ring.zero;
     else return Sum(v.all);
    end if;
  end Sum;

-- ARITHMETIC AS PROCEDURES :

  procedure Add ( v1 : in out Link_to_Vector; v2 : in Link_to_Vector ) is
  begin
    if v2 = null then
      null;
    elsif v1 = null then
      Copy(v1 => v2,v2 => v1);
    else
      Add(v1.all,v2.all);
    end if;
  end Add;

  procedure Min ( v : in Link_to_Vector ) is
  begin
    if v = null
     then null;
     else Min(v.all);
    end if;
  end Min;

  procedure Sub ( v1 : in out Link_to_Vector; v2 : in Link_to_Vector ) is
  begin
    if v2 = null then
      null;
    elsif v1 = null then
      v1 := new Vector'(v2.all);
      Min(v1.all);
    else
      Sub(v1.all,v2.all);
    end if;
  end Sub;

  procedure Mul ( v : in Link_to_Vector; a : in number ) is
  begin
    if v /= null
     then Mul(v.all,a);
    end if;
  end Mul;

  procedure Mul ( v1 : in out Link_to_Vector; v2 : in Link_to_Vector ) is
  begin
    if v2 = null then
      null;
    elsif v1 = null then
      Clear(v1);
    else
      Mul(v1.all,v2.all);
    end if;
  end Mul;

-- DESTRUCTOR :

  procedure Clear ( v : in out Link_to_Vector ) is

    procedure free is new unchecked_deallocation(Vector,Link_to_Vector);

  begin
    if v /= null then
      Clear(v.all);
      free(v);
    end if;
  end Clear;

end Generic_Vectors;
