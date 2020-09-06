with Standard_Integer_Vectors;
with Standard_Integer64_Vectors;
with Standard_Integer_Linear_Solvers;
with Standard_Integer64_Linear_Solvers;

package body Standard_Integer_Matrix_Inverse is

  function Is_Identity 
             ( A : Standard_Integer_Matrices.Matrix ) return boolean is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        if i = j then
          if A(i,j) /= 1
           then return false;
          end if;
        else
          if A(i,j) /= 0 
           then return false;
          end if;
        end if;
      end loop;
    end loop;
    return true;
  end Is_Identity;

  function Is_Identity 
             ( A : Standard_Integer64_Matrices.Matrix ) return boolean is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        if i = j then
          if A(i,j) /= 1
           then return false;
          end if;
        else
          if A(i,j) /= 0 
           then return false;
          end if;
        end if;
      end loop;
    end loop;
    return true;
  end Is_Identity;

  function Is_Inverse_Pair
             ( A,B : Standard_Integer_Matrices.Matrix ) return boolean is
 
    use Standard_Integer_Matrices;
    AB : constant Matrix(A'range(1),B'range(2)) := A*B;
    BA : constant Matrix(B'range(1),A'range(2)) := B*A;

  begin 
    if not Is_Identity(AB) or not Is_Identity(BA)
     then return false;
     else return true;
    end if;
  end Is_Inverse_Pair;

  function Is_Inverse_Pair
             ( A,B : Standard_Integer64_Matrices.Matrix ) return boolean is
 
    use Standard_Integer64_Matrices;
    AB : constant Matrix(A'range(1),B'range(2)) := A*B;
    BA : constant Matrix(B'range(1),A'range(2)) := B*A;

  begin 
    if not Is_Identity(AB) or not Is_Identity(BA)
     then return false;
     else return true;
    end if;
  end Is_Inverse_Pair;

  function Inverse_of_Unimodular_Upper
             ( A : in Standard_Integer_Matrices.Matrix )
             return Standard_Integer_Matrices.Matrix is

    use Standard_Integer_Matrices;
    use Standard_Integer_Linear_Solvers;

    res : Matrix(A'range(1),A'range(2));
    B : Matrix(A'range(1),A'first(2)..A'last(2)+1);
    x : Standard_Integer_Vectors.Vector(B'range(2));
 
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        B(i,j) := A(i,j);
      end loop;
    end loop;
    for k in A'range(2) loop
      for i in B'range(1) loop
        B(i,B'last(2)) := 0;
      end loop;
      B(k,B'last(2)) := -1;
      Solve0(B,x);
      if x(x'last) > 0 then
        for i in res'range(1) loop
          res(i,k) := x(i);
        end loop;
      else
        for i in res'range(1) loop
          res(i,k) := -x(i);
        end loop;
      end if;
    end loop;
    return res;
  end Inverse_of_Unimodular_Upper;

  function Inverse_of_Unimodular_Upper
             ( A : in Standard_Integer64_Matrices.Matrix )
             return Standard_Integer64_Matrices.Matrix is

    use Standard_Integer64_Matrices;
    use Standard_Integer64_Linear_Solvers;

    res : Matrix(A'range(1),A'range(2));
    B : Matrix(A'range(1),A'first(2)..A'last(2)+1);
    x : Standard_Integer64_Vectors.Vector(B'range(2));
 
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        B(i,j) := A(i,j);
      end loop;
    end loop;
    for k in A'range(2) loop
      for i in B'range(1) loop
        B(i,B'last(2)) := 0;
      end loop;
      B(k,B'last(2)) := -1;
      Solve0(B,x);
      if x(x'last) > 0 then
        for i in res'range(1) loop
          res(i,k) := x(i);
        end loop;
      else
        for i in res'range(1) loop
          res(i,k) := -x(i);
        end loop;
      end if;
    end loop;
    return res;
  end Inverse_of_Unimodular_Upper;
 
  procedure Inverse ( A : in Standard_Integer_Matrices.Matrix;
                      determinant : out integer32;
                      B : out Standard_Integer_Matrices.Matrix ) is

    use Standard_Integer_Matrices;
  
    d : integer32;
    U : Matrix(A'range(1),A'range(1));
    W : Matrix(A'range(1),A'range(2)) := A;
    R : Matrix(A'range(1),A'range(1));
  
  begin
    Standard_Integer_Linear_Solvers.Upper_Triangulate(U,W);
    d := W(W'first(1),W'first(2));
    for i in W'first(1)+1..W'last(1) loop
      d := d*W(i,i);
    end loop;
    determinant := d;
    if d*d = 1 then
      R := Inverse_of_Unimodular_Upper(W);
      B := R*U;
    end if;
  end Inverse;
 
  procedure Inverse ( A : in Standard_Integer64_Matrices.Matrix;
                      determinant : out integer64;
                      B : out Standard_Integer64_Matrices.Matrix ) is

    use Standard_Integer64_Matrices;
  
    d : integer64;
    U : Matrix(A'range(1),A'range(1));
    W : Matrix(A'range(1),A'range(2)) := A;
    R : Matrix(A'range(1),A'range(1));
  
  begin
    Standard_Integer64_Linear_Solvers.Upper_Triangulate(U,W);
    d := W(W'first(1),W'first(2));
    for i in W'first(1)+1..W'last(1) loop
      d := d*W(i,i);
    end loop;
    determinant := d;
    if d*d = 1 then
      R := Inverse_of_Unimodular_Upper(W);
      B := R*U;
    end if;
  end Inverse;

end Standard_Integer_Matrix_Inverse;
