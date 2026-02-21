with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Multprec_Complex_Numbers;          use Multprec_Complex_Numbers;
with Multprec_Complex_Vectors_IO;       use Multprec_Complex_Vectors_IO;
with Multprec_Complex_Singular_Values;  use Multprec_Complex_Singular_Values;

package body Multprec_Numerical_Rank is

  function Numerical_Rank
             ( S : Vector; tol : double_float;
               vrblvl : integer32 := 0 ) return natural32 is

    jump : constant double_float := SQRT(tol);
    tmp : Floating_Number;

  begin
    if vrblvl > 0 then
      put_line("-> in Multprec_Numerical_Rank.numerical_rank 0 ...");
      put_line("The singular values :"); put_line(S);
    end if;
    tmp := REAL_PART(S(S'first));
    if tmp < tol then
      Clear(tmp); return 0;
    else
      Clear(tmp);
      for i in S'first..S'last-1 loop
        tmp := REAL_PART(S(i+1))/REAL_PART(S(i));
        if tmp < jump
         then Clear(tmp); return natural32(i);
        end if;
        Clear(tmp);
      end loop;
      return natural32(S'last);
    end if;
  end Numerical_Rank;

  procedure Numerical_Rank
              ( A : in out Matrix; tol : in double_float;
                S : out Vector; U,V : out Matrix;
                rco : out Floating_Number; rank : out natural32;
                vrblvl : in integer32 := 0 ) is

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    e : Vector(1..p);
    info : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in Multprec_Numerical_Rank.numerical_rank 1 ...");
    end if;
    SVD(A,n,p,S,e,U,V,11,info);
    rco := REAL_PART(S(S'last))/REAL_PART(S(S'first));
    rank := Numerical_Rank(S,tol,vrblvl-1);
  end Numerical_Rank;

  function Numerical_Rank
             ( A : Matrix; tol : double_float;
               vrblvl : integer32 := 0 ) return natural32 is

    res : natural32;
    X : Matrix(A'range(1),A'range(2));
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
      put_line("-> in Multprec_Numerical_Rank.numerical_rank 2 ...");
    end if;
    Copy(A,X);
    SVD(X,n,p,S,e,U,V,11,info);
    res := Numerical_Rank(S,tol,vrblvl-1);
    Clear(U); Clear(V);
    Clear(e); Clear(S);
    return res;
  end Numerical_Rank;

end Multprec_Numerical_Rank;
