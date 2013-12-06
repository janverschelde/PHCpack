--with Multprec64_Common_Divisors;        use Multprec64_Common_Divisors;
--with Multprec_Integer64_Linear_Solvers; use Multprec_Integer64_Linear_Solvers;
with Multprec_Common_Divisors;           use Multprec_Common_Divisors;
with Multprec_Integer_Linear_Solvers;    use Multprec_Integer_Linear_Solvers;
with Multprec_Lattice_Supports;

package body Multprec_Integer_Orthogonals is

  function gcd ( A : Matrix; k : integer32 ) return Integer_Number is

    res,acc : Integer_Number;

  begin
    Copy(A(A'first(1),k),res);
    for i in A'first(1)+1..A'last(1) loop
      acc := gcd(res,A(i,k));
      Copy(acc,res); Clear(acc);
      exit when (Equal(res,1));
    end loop;
    return res;
  end gcd;

  procedure Normalize ( A : in out Matrix; k : in integer32 ) is

    g : Integer_Number := gcd(A,k);

  begin
    if not Equal(g,1) and not Equal(g,0) then
      for i in A'range(1) loop
        Div(A(i,k),g);
      end loop;
    end if;
    Clear(g);
  end Normalize;

  function Orthogonalize ( A : Matrix ) return Matrix is

    res : Matrix(A'range(1),A'range(2));
    p,q,acc : Integer_Number;

  begin
    Copy(A,res);
    Normalize(res,res'first(2));
    for k in A'first(2)+1..A'last(2) loop
      Normalize(res,k);
      for j in A'first(2)..(k-1) loop
        p := Multprec_Lattice_Supports.Inner_Product(res,j,j);
        q := Multprec_Lattice_Supports.Inner_Product(res,j,k);
        for i in res'range(1) loop
          Mul(res(i,k),p);
          acc := q*res(i,j);
          Sub(res(i,k),acc);
          Clear(acc);
        end loop;
        Normalize(res,k);
        Clear(p); Clear(q);
      end loop;
      exit when (k > A'last(1));
    end loop;
    return res;
  end Orthogonalize;

  function Complement  ( A : Matrix ) return Vector is

    res : Vector(A'range(1));
    B : Matrix(A'range(2),A'range(1));

  begin
    for i in A'range(1) loop
      res(i) := Multprec_Integer_Numbers.Create(integer32(0));
    end loop;
    for i in B'range(1) loop
      for j in B'range(2) loop
        Copy(A(j,i),B(i,j));
      end loop;
    end loop;
    Upper_Triangulate(B);
    Scale(B);
    Solve0(B,res);
    Clear(B);
    return res;
  end Complement;

end Multprec_Integer_Orthogonals;
