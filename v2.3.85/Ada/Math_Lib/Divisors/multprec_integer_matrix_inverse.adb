with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Integer_Vectors;
with Multprec_Integer_Linear_Solvers;    use Multprec_Integer_Linear_Solvers;  

package body Multprec_Integer_Matrix_Inverse is

  function Is_Identity 
             ( A : Multprec_Integer_Matrices.Matrix ) return boolean is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        if i = j then
          if not Equal(A(i,j),1)
           then return false;
          end if;
        else
          if not Equal(A(i,j),0)
           then return false;
          end if;
        end if;
      end loop;
    end loop;
    return true;
  end Is_Identity;

  function Is_Inverse_Pair
             ( A,B : Multprec_Integer_Matrices.Matrix ) return boolean is
 
    res : boolean;
    use Multprec_Integer_Matrices;
    AB : Matrix(A'range(1),B'range(2)) := A*B;
    BA : Matrix(B'range(1),A'range(2)) := B*A;

  begin 
    res := Is_Identity(AB) and Is_Identity(BA);
    Clear(AB); Clear(BA);
    return res;
  end Is_Inverse_Pair;

  function Inverse_of_Unimodular_Upper
             ( A : in Multprec_Integer_Matrices.Matrix )
             return Multprec_Integer_Matrices.Matrix is

    use Multprec_Integer_Matrices;

    res : Matrix(A'range(1),A'range(2));
    B : Matrix(A'range(1),A'first(2)..A'last(2)+1);
    x : Multprec_Integer_Vectors.Vector(B'range(2));
 
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        Copy(A(i,j),B(i,j));
      end loop;
    end loop;
    for k in A'range(2) loop
      for i in B'range(1) loop
        Clear(B(i,B'last(2)));
      end loop;
      for i in B'range(1) loop
        B(i,B'last(2)) := Create(integer(0));
      end loop;
      B(k,B'last(2)) := Create(-integer(1));
      Solve0(B,x);
      if x(x'last) > 0 then
        for i in res'range(1) loop
          Copy(x(i),res(i,k)); Clear(x(i));
        end loop;
      else
        for i in res'range(1) loop
          Copy(x(i),res(i,k)); Clear(x(i)); Min(res(i,k));
        end loop;
      end if;
    end loop;
    return res;
  end Inverse_of_Unimodular_Upper;
 
  procedure Inverse ( A : in Multprec_Integer_Matrices.Matrix;
                      determinant : out Integer_Number;
                      B : out Multprec_Integer_Matrices.Matrix ) is

    use Multprec_Integer_Matrices;
  
    d : Integer_Number;
    U,W,R : Matrix(A'range(1),A'range(1));
  
  begin
    Copy(A,W);
    Upper_Triangulate(U,W);
    Copy(W(W'first(1),W'first(2)),d);
    for i in W'first(1)+1..W'last(1) loop
      Mul(d,W(i,i));
    end loop;
    Copy(d,determinant);
    if Equal(d,1) or Equal(d,-1) then
      R := Inverse_of_Unimodular_Upper(W);
      B := R*U; Clear(R);
    end if;
    Clear(d); Clear(W); Clear(U);
  end Inverse;

end Multprec_Integer_Matrix_Inverse;
