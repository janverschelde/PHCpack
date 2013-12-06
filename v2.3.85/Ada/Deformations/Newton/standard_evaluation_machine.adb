with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Singular_Values;  use Standard_Complex_Singular_Values;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Numerical_Derivatives;    use Standard_Numerical_Derivatives;
with Evaluated_Minors;                  use Evaluated_Minors;

package body Standard_Evaluation_Machine is

  p_eval : Eval_Poly;
  s_eval : Link_to_Eval_Poly_Sys;
  count_deflations : integer32;

  procedure Initialize ( p : in Poly ) is
  begin
    p_eval := Create(p);
  end Initialize;

  procedure Initialize ( p : in Poly_Sys ) is
  begin
    s_eval := new Eval_Poly_Sys'(Create(p));
  end Initialize;

  procedure Deflate is
  begin
    count_deflations := count_deflations + 1;
  end Deflate;

  function Evaluate ( x : Vector ) return Complex_Number is
  begin
    return Eval(p_eval,x);
  end Evaluate;

  function Poly_Eval ( x : Vector ) return Vector is
  begin
    return Eval(s_eval.all,x);
  end Poly_Eval;

  function Jacobian_Determinant
             ( n : integer32; x : Vector ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the determinant of the Jacobian matrix at x.

    h : constant double_float := 0.001;
    jm : constant Matrix(1..n,x'range) := Diff3(Poly_Eval'access,n,x,h);

  begin
    return Determinant(jm);
  end Jacobian_Determinant;

  function Inverse_Condition_Number
             ( n : integer32; x : Vector ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the inverse of the condition number at x,
  --   for a system of n equations.

    res : Complex_Number;
    h : constant double_float := 0.001;
    jm : Matrix(1..n,x'range) := Diff3(Poly_Eval'access,n,x,h);
    u : Matrix(1..n,1..n);
    v : Matrix(x'range,x'range);
    d : constant integer32 := x'length;
    m : constant integer32 := Min0(n+1,d);
    e : Vector(1..d);
    s : Vector(1..m);
    info : integer32;

  begin
    SVD(jm,n,d,s,e,u,v,11,info);
    res := s(m)/s(1);
    return res;
  end Inverse_Condition_Number;

  function Evaluate ( x : Vector ) return Vector is

    y : constant Vector := Eval(s_eval.all,x);

  begin
    if count_deflations = 0 then
      return y;
    else
      declare
        res : Vector(y'first..y'last+1);
        n : constant integer32 := s_eval'last;
      begin
        res(y'range) := y;
        -- res(n+1) := Inverse_Condition_Number(n,x);
        res(n+1) := Jacobian_Determinant(n,x);
        return res;
      end;
    end if;
  end Evaluate;

  procedure Clear is
  begin
    Clear(p_eval);
    Clear(s_eval);
  end Clear;

begin
  count_deflations := 0;
end Standard_Evaluation_Machine;
