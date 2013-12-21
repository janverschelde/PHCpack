with Standard_Integer_Norms;            use Standard_Integer_Norms;
with Standard_Integer64_Linear_Solvers; use Standard_Integer64_Linear_Solvers;
with Standard_Integer_Kernel;

package body Standard_Integer64_Kernel is

  function Rank_Upper ( U : Matrix ) return integer32 is

    res : integer32 := 0;
    pivot : integer32 := U'first(2);

  begin
    for i in U'range(1) loop
      while pivot <= U'last(2) and then U(i,pivot) = 0 loop
        pivot := pivot + 1;
      end loop;
      exit when pivot > U'last(2);
      res := i;
    end loop;
    return res;
  end Rank_Upper;

  procedure Pivots_in_Upper ( U : in Matrix; rank : out integer32;
                              p : out Standard_Integer_Vectors.Vector ) is

    pivot : integer32 := U'first(2);

  begin
    for i in p'range loop
      p(i) := 0;
    end loop;
    for i in U'range(1) loop
      while pivot <= U'last(2) and then U(i,pivot) = 0 loop
        pivot := pivot + 1;
      end loop;
      exit when pivot > U'last(2);
      p(i) := pivot;
      rank := i;
    end loop;
  end Pivots_in_Upper;

  procedure Normalize_Sign ( v : in out Standard_Integer64_Vectors.Vector ) is
  begin
    for i in v'range loop
      if v(i) > 0 then
        return;
      elsif v(i) < 0 then
        exit;
      end if;
    end loop;
    Standard_Integer64_Vectors.Min(v);
  end Normalize_Sign;

  procedure Upper_Kernel ( U : in Matrix; n,r : in integer32;
                           p : in Standard_Integer_Vectors.Vector;
                           V : out Link_to_Matrix ) is

    W : Matrix(1..n,1..n-r);
    B : Matrix(1..r,1..r+1);
    x : Standard_Integer64_Vectors.Vector(B'range(2));
    q : constant Standard_Integer_Vectors.Vector(1..n-r)
      := Standard_Integer_Kernel.Complement(n,r,p);

  begin
    for j in 1..r loop                -- copy pivot column j of U to B
      for i in 1..r loop
        B(i,j) := U(i,p(j));
      end loop;
    end loop;
    for k in W'range(2) loop          -- compute k-th basis vector
      for i in B'range(1) loop        -- select k-th nonpivot
        B(i,B'last(2)) := U(i,q(k));  -- as right hand side vector
      end loop;
      Solve0(B,x);
      Normalize(x);                   -- gcd of components equals 1
      Normalize_Sign(x);
      for i in 1..n loop
        W(i,k) := 0;
      end loop;
      for i in 1..r loop
        W(p(i),k) := x(i);
      end loop;
      W(q(k),k) := x(x'last);
    end loop;
    V := new Matrix'(W);
  end Upper_Kernel;

  procedure Kernel ( A : in Matrix;
                     r : out integer32; V : out Link_to_Matrix ) is

    n : constant integer32 := A'last(2);
    U : Matrix(A'range(1),A'range(2)) := A;
    pivots : Standard_Integer_Vectors.Vector(A'range(1));

  begin
    Upper_Triangulate(U);
    r := Rank_Upper(U);
    Pivots_in_Upper(U,r,pivots);
    if r < n
     then Upper_Kernel(U,n,r,pivots,V);
    end if;
  end Kernel;

end Standard_Integer64_Kernel;
