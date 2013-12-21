with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Linear_Solvers;   use Standard_Complex_Linear_Solvers;

package body Standard_Cascading_Planes is

  function Project ( x : Matrix; h : integer32 ) return Matrix is

    res : Matrix(x'range(1),x'range(2));

  begin
    for i in x'first(1)..x'first(1)+h-1 loop
      for j in x'range(2) loop
        res(i,j) := x(i,j);
      end loop;
    end loop;
    for i in x'first(1)+h..x'last(1) loop
      for j in x'range(2) loop
        res(i,j) := Create(0.0);
      end loop;
    end loop;
    return res;
  end Project;

  function Project ( x : Vector; h : integer32 ) return Vector is

    res : Vector(x'range);

  begin
    for i in x'first(1)..x'first(1)+h-1 loop
      res(i) := x(i);
    end loop;
    for i in x'first(1)+h..x'last(1) loop
      res(i) := Create(0.0);
    end loop;
    return res;
  end Project;

  function Double_to_Diagonal ( A : Matrix ) return Matrix is

    n : constant integer32 := A'last(2);
    res : Matrix(A'range(1),1..2*n);

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := A(i,j);
        res(i,n+j) := -A(i,j);
      end loop;
    end loop;
    return res;
  end Double_to_Diagonal;

  function Compute_Offset ( C : Matrix; d : Vector ) return Vector is

    n : constant integer32 := d'last;
    n2 : constant integer32 := 2*n;
    res : Vector(C'range(2));
    mat : Matrix(C'range(2),C'range(2));
    piv : Standard_Integer_Vectors.Vector(C'range(2));
    inf : integer32;

  begin
    for i in 1..n loop
      for j in 1..n loop
        if i = j
         then mat(i,j) := Create(1.0); mat(i,j+n) := Create(-1.0);
         else mat(i,j) := Create(0.0); mat(i,j+n) := Create(0.0);
        end if;
      end loop;
      for j in C'range(2) loop
        mat(i+n,j) := C(i,j);
      end loop;
      res(i) := Create(0.0);
      res(i+n) := -d(i);
    end loop;
    lufac(mat,n2,piv,inf);
    lusolve(mat,n2,piv,res);
    return res;
  end Compute_Offset;

  procedure Shift_Offset ( p : in out Matrix; b : in Vector ) is
  begin
    for i in p'range(1) loop
      p(i,0) := b(i);
    end loop;
  end Shift_Offset;

  function Target_Space ( n,n2,apb,b : integer32; dA,BB,CC : Matrix;
                          d : Vector ) return Matrix is

    PC : constant Matrix(1..n,1..n2) := Project(CC,b);
    BPC : constant Matrix(1..apb,1..n2) := BB*PC;
    Pd : constant Vector(1..n) := Project(d,b);
    BPd : constant Vector(1..apb) := BB*Pd;
    eqs : Matrix(1..apb,0..n2);

  begin
    for i in eqs'range(1) loop
      eqs(i,0) := BPd(i);
    end loop;
    for i in eqs'range(1) loop
      for j in 1..n2 loop
        eqs(i,j) := dA(i,j) + BPC(i,j);
      end loop;
    end loop;
    return eqs;
  end Target_Space;

  function Start_Space ( g1,g2 : Matrix ) return Matrix is

    n : constant integer32 := g1'last(1);
    n2 : constant integer32 := 2*n;
    res : Matrix(1..n2,0..g1'last(2)+g2'last(2));

  begin
    for i in 1..n loop
      res(i,0) := g1(i,0);
      res(i+n,0) := g2(i,0);
      for j in 1..g1'last(2) loop
        res(i,j) := g1(i,j);
        res(i+n,j) := Create(0.0);
      end loop;
      for j in 1..g2'last(2) loop
        res(i,g1'last(2)+j) := Create(0.0);
        res(i+n,g1'last(2)+j) := g2(i,j);
      end loop;
    end loop;
    return res;
  end Start_Space;

end Standard_Cascading_Planes;
