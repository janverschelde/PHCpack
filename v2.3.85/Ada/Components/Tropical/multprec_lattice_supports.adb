package body Multprec_Lattice_Supports is

  function Inner_Product
              ( u,v : Multprec_Integer_Vectors.Vector )
              return Integer_Number is

    res : Integer_Number := Multprec_Integer_Numbers.Create(integer32(0));
    acc : Integer_Number;

  begin
    for i in u'range loop
      if not Equal(u(i),0) then
        if not Equal(v(i),0) then
          acc := u(i)*v(i);
          Add(res,acc);
          Clear(acc);
        end if;
      end if;
    end loop;
    return res;
  end Inner_Product;

  function Inner_Product
              ( A : Multprec_Integer_Matrices.Matrix; k : integer32;
                v : Multprec_Integer_Vectors.Vector )
              return Integer_Number is

    res : Integer_Number := Multprec_Integer_Numbers.Create(integer32(0));
    acc : Integer_Number;

  begin
    for i in A'range(1) loop
      if not Equal(A(i,k),0) then
        if not Equal(v(i),0) then
          acc := A(i,k)*v(i);
          Add(res,acc);
          Clear(acc);
        end if;
      end if;
    end loop;
    return res;
  end Inner_Product;

  function Inner_Product
              ( A : Multprec_Integer_Matrices.Matrix; i,j : integer32 )
              return Integer_Number is

    res : Integer_Number := Multprec_Integer_Numbers.Create(integer32(0));
    acc : Integer_Number;

  begin
    for k in A'range(1) loop
      if not Equal(A(k,i),0) then
        if not Equal(A(i,j),0) then
          acc := A(k,i)*A(k,j);
          Add(res,acc);
          Clear(acc);
        end if;
      end if;
    end loop;
    return res;
  end Inner_Product;

  function Inner_Products
              ( A : Multprec_Integer_Matrices.Matrix;
                v : Multprec_Integer_Vectors.Vector )
              return Multprec_Integer_Vectors.Vector is

    res : Multprec_Integer_Vectors.Vector(A'range(2));

  begin
    for k in A'range(2) loop
      res(k) := Inner_Product(A,k,v);
    end loop;
    return res;
  end Inner_Products;

  function Minimum ( A : Multprec_Integer_Matrices.Matrix;
                     v : Multprec_Integer_Vectors.Vector )
                   return Integer_Number is

    res : Integer_Number := Inner_Product(A,A'first(2),v);
    wrk : Integer_Number;

  begin
    for k in A'first(2)+1..A'last(2) loop
      wrk := Inner_Product(A,k,v);
      if wrk < res
       then Copy(wrk,res);
      end if;
      Clear(wrk);
    end loop;
    return res;
  end Minimum;

  function Support ( A : Multprec_Integer_Matrices.Matrix;
                     v : Multprec_Integer_Vectors.Vector )
                   return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(A'range(2));
    ind : integer32 := A'first(2);
    min : Integer_Number := Inner_Product(A,A'first(2),v);
    wrk : Integer_Number;

  begin
    res(ind) := A'first(2);
    for k in A'first(2)+1..A'last(2) loop
      wrk := Inner_Product(A,k,v);
      if Equal(wrk,min) then
        ind := ind + 1;
        res(ind) := k;
      elsif wrk < min then
        ind := A'first(2);
        res(ind) := k;
        Copy(wrk,min);
      end if;
      Clear(wrk);
    end loop;
    Clear(min);
    return res(res'first..ind);
  end Support;

  function Support ( A,B : Multprec_Integer_Matrices.Matrix;
                     v : Multprec_Integer_Vectors.Vector )
                   return Standard_Integer_VecVecs.VecVec is

    res : Standard_Integer_VecVecs.VecVec(1..2);

  begin
    res(1) := new Standard_Integer_Vectors.Vector'(Support(A,v));
    res(2) := new Standard_Integer_Vectors.Vector'(Support(B,v));
    return res;
  end Support;

  procedure Inner ( A : in Multprec_Integer_Matrices.Matrix;
                    i,j : in integer32;
                    v : in out Multprec_Integer_Vectors.Vector ) is

    p : integer32 := A'first(2)-1;
    s : Integer_Number := Inner_Product(A,i,v);
    t : Integer_Number;
 
  begin
    while p < A'last(2) loop
      p := p + 1;
      if ((p /= i) and (p /= j)) then
        t := Inner_Product(A,p,v);
        if t < s
         then Multprec_Integer_Vectors.Min(v);
        end if;
        exit when (not Equal(s,t));
        Clear(t);
      end if;
    end loop;
    Clear(s);
  end Inner;

  procedure Inner ( A : in Multprec_Integer_Matrices.Matrix;
                    i,j : in integer32;
                    f : in Standard_Integer_Vectors.Vector;
                    v : in out Multprec_Integer_Vectors.Vector ) is

    s : Integer_Number := Inner_Product(A,i,v);
    t : Integer_Number;

  begin
    for k in f'range loop
      if f(k) /= i and f(k) /= j then
        t := Inner_Product(A,f(k),v);
        if t < s
         then Multprec_Integer_Vectors.Min(v);
        end if;
        exit when (not Equal(s,t));
        Clear(t);
      end if;
    end loop;
    Clear(s);
  end inner;

  procedure Inner ( A : in Multprec_Integer_Matrices.Matrix;
                    i,j,k : in integer32;
                    v : in out Multprec_Integer_Vectors.Vector ) is

    p : integer32 := A'first(2)-1;
    s : Integer_Number := Inner_Product(A,i,v);
    t : Integer_Number;

  begin
    while p < A'last(2) loop
      p := p + 1;
      if ((p /= i) and (p /= j) and (p /= k)) then
        t := Inner_Product(A,p,v);
        if t < s
         then Multprec_Integer_Vectors.Min(v);
        end if;
        exit when (not Equal(s,t));
        Clear(t);
      end if;
    end loop;
    Clear(s);
  end Inner;

  function Support_Points
              ( A : Multprec_Integer_Matrices.Matrix;
                s : Standard_Integer_Vectors.Vector )
              return Multprec_Integer_Matrices.Matrix is

    res : Multprec_Integer_Matrices.Matrix(A'range(1),s'range);

  begin
    for j in s'range loop
      for i in A'range(1) loop
        Copy(A(i,s(j)),res(i,j));
      end loop;
    end loop;
    return res;
  end Support_Points;

  function Equal ( A,B : Multprec_Integer_Matrices.Matrix; i,j : integer32 )
                 return boolean is
  begin
    for k in A'range(1) loop
      if not Equal(A(k,i),B(k,j))
       then return false;
      end if;
    end loop;
    return true;
  end Equal;

  function Indices ( A,V : Multprec_Integer_Matrices.Matrix )
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

end Multprec_Lattice_Supports;
