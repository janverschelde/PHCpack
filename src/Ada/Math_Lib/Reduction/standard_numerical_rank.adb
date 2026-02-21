with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors_IO;       use Standard_Complex_Vectors_IO;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Complex_Singular_Values;  use Standard_Complex_Singular_Values;

package body Standard_Numerical_Rank is

  function Numerical_Rank
             ( S : Vector; tol : double_float;
               vrblvl : integer32 := 0 ) return natural32 is

    jump : constant double_float := SQRT(tol);

  begin
    if vrblvl > 0 then
      put_line("-> in Standard_Numerical_Rank.numerical_rank 0 ...");   
      put_line("The singular values : "); put_line(S);
    end if;
    if REAL_PART(S(S'first)) < tol then
      return 0;
    else
      for i in S'first..S'last-1 loop
        if REAL_PART(S(i+1))/REAL_PART(S(i)) < jump
         then return natural32(i);
        end if;
      end loop;
      return natural32(S'last);
    end if;
  end Numerical_Rank;

  procedure Numerical_Rank 
              ( A : in out Matrix; tol : in double_float;
                S : out Vector; U,V : out Matrix;
                rco : out double_float; rank : out natural32;
                vrblvl : in integer32 := 0 ) is

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    e : Vector(1..p);
    info : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in Standard_Numerical_Rank.numerical_rank 1 ...");   
    end if;
    SVD(A,n,p,S,e,U,V,11,info);
    rco := REAL_PART(S(S'last))/REAL_PART(S(S'first));
    rank := Numerical_Rank(S,tol,vrblvl-1);
  end Numerical_Rank;

  function Numerical_Rank
             ( A : Matrix; tol : double_float;
               vrblvl : integer32 := 0 ) return natural32 is

    X : Matrix(A'range(1),A'range(2)) := A;
    U : Matrix(A'range(1),A'range(1));
    V : Matrix(A'range(2),A'range(2));
    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    m : constant integer32 := Min0(n+1,p);
    S : Vector(1..m);
    e : Vector(1..p);
    info : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in Standard_Numerical_Rank.numerical_rank 2 ...");   
    end if;
    SVD(X,n,p,S,e,U,V,11,info);
    return Numerical_Rank(S,tol,vrblvl-1);
  end Numerical_Rank;

end Standard_Numerical_Rank;
