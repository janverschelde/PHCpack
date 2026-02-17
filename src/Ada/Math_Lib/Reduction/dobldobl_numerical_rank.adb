with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with DoblDobl_Complex_Singular_Values;  use DoblDobl_Complex_Singular_Values;

package body DoblDobl_Numerical_Rank is

  function Numerical_Rank
             ( S : Vector; tol : double_float ) return natural32 is

    jump : constant double_float := SQRT(tol);

  begin
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
                rco : out double_double; rank : out natural32 ) is

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    e : Vector(1..p);
    info : integer32;

  begin
    SVD(A,n,p,S,e,U,V,11,info);
    rco := REAL_PART(S(S'last))/REAL_PART(S(S'first));
    rank := Numerical_Rank(S,tol);
  end Numerical_Rank;

  function Numerical_Rank
             ( A : Matrix; tol : double_float ) return natural32 is

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
    SVD(X,n,p,S,e,U,V,11,info);
    return Numerical_Rank(S,tol);
  end Numerical_Rank;

end DoblDobl_Numerical_Rank;
