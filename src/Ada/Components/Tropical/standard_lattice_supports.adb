package body Standard_Lattice_Supports is

  function Inner_Product
              ( u,v : Standard_Integer64_Vectors.Vector ) return integer64 is

    res : integer64 := 0;

  begin
    for i in u'range loop
      res := res + u(i)*v(i);
    end loop;
    return res;
  end Inner_Product;

  function Inner_Product
              ( A : Standard_Integer64_Matrices.Matrix; k : integer32;
                v : Standard_Integer64_Vectors.Vector ) return integer64 is

    res : integer64 := 0;

  begin
    for i in A'range(1) loop
      res := res + A(i,k)*v(i);
    end loop;
    return res;
  end Inner_Product;

  function Inner_Product
              ( A : Standard_Integer64_Matrices.Matrix; i,j : integer32 )
              return integer64 is

    res : integer64 := 0;

  begin
    for k in A'range(1) loop
      res := res + A(k,i)*A(k,j);
    end loop;
    return res;
  end Inner_Product;

  function Inner_Products
              ( A : Standard_Integer64_Matrices.Matrix;
                v : Standard_Integer64_Vectors.Vector )
              return Standard_Integer64_Vectors.Vector is

    res : Standard_Integer64_Vectors.Vector(A'range(2));

  begin
    for k in A'range(2) loop
      res(k) := Inner_Product(A,k,v);
    end loop;
    return res;
  end Inner_Products;

  function Minimum ( A : Standard_Integer64_Matrices.Matrix;
                     v : Standard_Integer64_Vectors.Vector )
                   return integer64 is

    res : integer64 := Inner_Product(A,A'first(2),v);
    wrk : integer64;

  begin
    for k in A'first(2)+1..A'last(2) loop
      wrk := Inner_Product(A,k,v);
      if wrk < res
       then res := wrk;
      end if;
    end loop;
    return res;
  end Minimum;

  function Support ( A : Standard_Integer64_Matrices.Matrix;
                     v : Standard_Integer64_Vectors.Vector )
                   return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(A'range(2));
    ind : integer32 := A'first(2);
    min : integer64 := Inner_Product(A,A'first(2),v);
    wrk : integer64;

  begin
    res(ind) := A'first(2);
    for k in A'first(2)+1..A'last(2) loop
      wrk := Inner_Product(A,k,v);
      if wrk = min then
        ind := ind + 1;
        res(ind) := k;
      elsif wrk < min then
        ind := A'first(2);
        res(ind) := k;
        min := wrk;
      end if;
    end loop;
    return res(res'first..ind);
  end Support;

  function Support ( A,B : Standard_Integer64_Matrices.Matrix;
                     v : Standard_Integer64_Vectors.Vector )
                   return Standard_Integer_VecVecs.VecVec is

    res : Standard_Integer_VecVecs.VecVec(1..2);

  begin
    res(1) := new Standard_Integer_Vectors.Vector'(Support(A,v));
    res(2) := new Standard_Integer_Vectors.Vector'(Support(B,v));
    return res;
  end Support;

  procedure Inner ( A : in Standard_Integer64_Matrices.Matrix;
                    i,j : in integer32;
                    v : in out Standard_Integer64_Vectors.Vector ) is

    p : integer32 := A'first(2)-1;
    s : constant integer64 := Inner_Product(A,i,v);
    t : integer64;
 
  begin
    while p < A'last(2) loop
      p := p + 1;
      if ((p /= i) and (p /= j)) then
        t := Inner_Product(A,p,v);
        if t < s
         then Standard_Integer64_Vectors.Min(v);
        end if;
        exit when (t /= s);
      end if;
    end loop;
  end Inner;

  procedure Inner ( A : in Standard_Integer64_Matrices.Matrix;
                    i,j : in integer32;
                    f : in Standard_Integer_Vectors.Vector;
                    v : in out Standard_Integer64_Vectors.Vector ) is

    s : constant integer64 := Inner_Product(A,i,v);
    t : integer64;

  begin
    for k in f'range loop
      if f(k) /= i and f(k) /= j then
        t := Inner_Product(A,f(k),v);
        if t < s
         then Standard_Integer64_Vectors.Min(v);
        end if;
        exit when (t /= s);
      end if;
    end loop;
  end inner;

  procedure Inner ( A : in Standard_Integer64_Matrices.Matrix;
                    i,j,k : in integer32;
                    v : in out Standard_Integer64_Vectors.Vector ) is

    p : integer32 := A'first(2)-1;
    s : constant integer64 := Inner_Product(A,i,v);
    t : integer64;

  begin
    while p < A'last(2) loop
      p := p + 1;
      if ((p /= i) and (p /= j) and (p /= k)) then
        t := Inner_Product(A,p,v);
        if t < s
         then Standard_Integer64_Vectors.Min(v);
        end if;
        exit when (t /= s);
      end if;
    end loop;
  end Inner;

  function Support_Points
              ( A : Standard_Integer64_Matrices.Matrix;
                s : Standard_Integer_Vectors.Vector )
              return Standard_Integer64_Matrices.Matrix is

    res : Standard_Integer64_Matrices.Matrix(A'range(1),s'range);

  begin
    for j in s'range loop
      for i in A'range(1) loop
        res(i,j) := A(i,s(j));
      end loop;
    end loop;
    return res;
  end Support_Points;

  function Equal ( A,B : Standard_Integer64_Matrices.Matrix;
                   i,j : integer32 ) return boolean is
  begin
    for k in A'range(1) loop
      if A(k,i) /= B(k,j)
       then return false;
      end if;
    end loop;
    return true;
  end Equal;

  function Indices ( A,V : Standard_Integer64_Matrices.Matrix )
                   return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(V'range(2));

  begin
    for k in V'range(2) loop
      for j in A'range(2) loop
        if Equal(A,V,j,k)
         then res(k) := j; exit;
        end if;
      end loop;
    end loop;
    return res;
  end Indices;

  function Member ( v : Standard_Integer_Vectors.Vector;
                    e : integer32 ) return integer32 is
  begin
    for k in v'range loop
      if v(k) = e
       then return k;
      end if;
    end loop;
    return (v'first-1);
  end Member;

end Standard_Lattice_Supports;
