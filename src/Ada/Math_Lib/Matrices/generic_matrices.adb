with unchecked_deallocation;

package body Generic_Matrices is

-- COMPARISON AND COPYING :

  function Equal ( a,b : Matrix ) return boolean is
  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        if not Equal(a(i,j),b(i,j))
         then return false;
        end if;
      end loop;
    end loop;
    return true;
  exception
    when CONSTRAINT_ERROR => return false;
  end Equal;

  procedure Copy ( a : in Matrix; b : in out Matrix ) is
  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        Copy(a(i,j),b(i,j));
      end loop;
    end loop;
  end Copy;

-- MATRIX-MATRIX OPERATIONS :

  function Transpose ( a : Matrix ) return Matrix is

    res : Matrix(a'range(2),a'range(1));

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        res(j,i) := a(i,j);
      end loop;
    end loop;
    return res;
  end Transpose;

  function "+" ( a,b : Matrix ) return Matrix is

    res : Matrix(a'range(1),a'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := a(i,j) + b(i,j);
      end loop;
    end loop;
    return res;
  end "+";

  function "-" ( a,b : Matrix ) return Matrix is

    res : Matrix(a'range(1),a'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := a(i,j) - b(i,j);
      end loop;
    end loop;
    return res;
  end "-";

  function "+" ( a : Matrix ) return Matrix is

    res : Matrix(a'range(1),a'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := +a(i,j);
      end loop;
    end loop;
    return res;
  end "+";

  function "-" ( a : Matrix ) return Matrix is

    res : Matrix(a'range(1),a'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := -a(i,j);
      end loop;
    end loop;
    return res;
  end "-";

  function "*" ( a,b : Matrix ) return Matrix is

    res : Matrix(a'range(1),b'range(2));
    acc : number;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := a(i,a'first(2))*b(b'first(1),j);
        for k in a'first(2)+1..a'last(2) loop
          acc := a(i,k)*b(k,j);
          Add(res(i,j),acc);
          Clear(acc);
        end loop;
      end loop;
    end loop;
    return res;
  end "*";

  procedure Add ( a : in out Matrix; b : in Matrix ) is
  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        Add(a(i,j),b(i,j));
      end loop;
    end loop;
  end Add;

  procedure Sub ( a : in out Matrix; b : in Matrix ) is
  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        Sub(a(i,j),b(i,j));
      end loop;
    end loop;
  end Sub;

  procedure Mul1 ( a : in out Matrix; b : in Matrix ) is

    temp : Vector(a'range(2));
    acc : number;

  begin
    for i in a'range(1) loop
      for j in b'range(2) loop
        temp(j) := a(i,a'first(2))*b(b'first(1),j);
        for k in a'first(2)+1..a'last(2) loop
          acc :=  a(i,k)*b(k,j);
          Add(temp(j),acc);
          Clear(acc);
        end loop;
      end loop;
      for j in a'range(2) loop
        Copy(temp(j),a(i,j));
      end loop;
    end loop;
  end Mul1;

  procedure Mul2 ( a : in Matrix; b : in out Matrix ) is

    temp : Vector(a'range(1));
    acc : number;

  begin
    for i in b'range(2) loop
      for j in a'range(1) loop
        temp(j) := a(j,a'first(1))*b(a'first(1),i);
        for k in a'first(1)+1..a'last(1) loop
          acc := a(j,k)*b(k,i);
          Add(temp(j),acc);
          Clear(acc);
        end loop;
      end loop;
      for j in b'range(1) loop
        Copy(temp(j),b(j,i));
      end loop;
    end loop;
  end Mul2;

-- MATRIX-VECTOR OPERATIONS :

  function "*" ( a : Matrix; v : Vector ) return Vector is

    res : Vector(a'range(1));
    acc : number;

  begin
    for i in res'range loop
      res(i) := a(i,a'first(2))*v(v'first);
      for j in a'first(2)+1..a'last(2) loop
        acc := a(i,j)*v(j);
        Add(res(i),acc);
        Clear(acc);
      end loop;
    end loop;
    return res;
  end "*";

  function "*" ( v : Vector; a : Matrix ) return Vector is

    res : Vector(a'range(2));
    acc : number;

  begin
    for j in res'range loop
      res(j) := v(v'first)*a(a'first(1),j);
      for i in a'first(1)+1..a'last(1) loop
        acc := v(i)*a(i,j);
        Add(res(j),acc);
        Clear(acc);
      end loop;
    end loop;
    return res;
  end "*";

  procedure Mul ( a : in Matrix; v : in out Vector ) is

    tv : Vector(v'range);
    acc : number;

  begin
    for i in v'range loop
      tv(i) := a(i,a'first(2))*v(v'first);
      for j in a'first(2)+1..a'last(2) loop
        acc := a(i,j)*v(j);
        Add(tv(i),a(i,j)*v(j));
        Clear(acc);
      end loop;
    end loop;
    for i in v'range loop
      v(i) := tv(i);
    end loop;
  end Mul;

  procedure Mul ( v : in out Vector; a : in Matrix ) is

    tv : Vector(v'range);
    acc : number;

  begin
    for j in v'range loop
      tv(j) := v(v'first)*a(a'first(1),j);
      for i in a'first(1)+1..a'last(1) loop
        acc := v(j)*a(i,j);
        Add(tv(j),acc);
        Clear(acc);
      end loop;
    end loop;
    for i in v'range loop
      v(i) := tv(i);
    end loop;
  end Mul;

-- SCALING A MATRIX :

  function "*" ( x : number; a : Matrix ) return Matrix is

    res : Matrix(a'range(1),a'range(2));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := x*a(i,j);
      end loop;
    end loop;
    return res;
  end "*";

  function "*" ( a : Matrix; x : number ) return Matrix is

    res : constant Matrix(a'range(1),a'range(2)) := x*a;

  begin
    return res;
  end "*";

  procedure Mul ( a : in out Matrix; x : in number ) is
  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        Mul(a(i,j),x);
      end loop;
    end loop;
  end Mul;

-- DESTRUCTORS :

  procedure Clear ( a : in out Matrix ) is
  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        Clear(a(i,j));
      end loop;
    end loop;
  end Clear;

  procedure Clear ( a : in out Link_to_Matrix ) is

    procedure free is new unchecked_deallocation(Matrix,Link_to_Matrix);

  begin
    if a /= null
     then Clear(a.all);
    end if;
    free(a);
  end Clear;

end Generic_Matrices;
