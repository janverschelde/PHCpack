with Standard64_Common_Divisors;         use Standard64_Common_Divisors;
with Standard_Integer64_Linear_Solvers;  use Standard_Integer64_Linear_Solvers;
with Standard_Lattice_Supports;

package body Standard_Integer_Orthogonals is

  function gcd ( A : Matrix; k : integer32 ) return integer64 is

    res : integer64 := A(A'first(1),k);

  begin
    for i in A'first(1)+1..A'last(1) loop
      res := gcd(res,A(i,k));
      exit when (res = 1);
    end loop;
    return res;
  end gcd;

  procedure Normalize ( A : in out Matrix; k : in integer32 ) is

    g : constant integer64 := gcd(A,k);

  begin
    if g /= 1 and g /= 0 then
      for i in A'range(1) loop
        A(i,k) := A(i,k)/g;
      end loop;
    end if;
  end Normalize;

  function Orthogonalize ( A : Matrix ) return Matrix is

    res : Matrix(A'range(1),A'range(2)) := A;
    p,q : integer64;

  begin
    Normalize(res,res'first(2));
    for k in A'first(2)+1..A'last(2) loop
      Normalize(res,k);
      for j in A'first(2)..(k-1) loop
        p := Standard_Lattice_Supports.Inner_Product(res,j,j);
        q := Standard_Lattice_Supports.Inner_Product(res,j,k);
        for i in res'range(1) loop
          res(i,k) := p*res(i,k) - q*res(i,j);
        end loop;
        Normalize(res,k);
      end loop;
      exit when (k > A'last(1));
    end loop;
    return res;
  end Orthogonalize;

  function Complement  ( A : Matrix ) return Vector is

    res : Vector(A'range(1)) := (A'range(1) => 0);
    B : Matrix(A'range(2),A'range(1));

  begin
    for i in B'range(1) loop
      for j in B'range(2) loop
        B(i,j) := A(j,i);
      end loop;
    end loop;
    Upper_Triangulate(B);
    Scale(B);
    Solve0(B,res);
    return res;
  end Complement;

end Standard_Integer_Orthogonals;
