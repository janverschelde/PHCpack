with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Multprec_Complex_Numbers;          use Multprec_Complex_Numbers;
with Multprec_Complex_Vectors_io;       use Multprec_Complex_Vectors_io;
with Multprec_Complex_Norms_Equals;     use Multprec_Complex_Norms_Equals;
with Multprec_Complex_Singular_Values;  use Multprec_Complex_Singular_Values;
with Multprec_Numerical_Rank;           use Multprec_Numerical_Rank;

package body Multprec_Complex_Newton_Steps is

  function Inverse_Condition_Number ( s : Vector ) return Floating_Number is

    res : Floating_Number := REAL_PART(s(s'last));
    den : Floating_Number := REAL_PART(s(s'first));
    one : Floating_Number := Create(natural32(1));
    den_one : Floating_Number := one + den;

  begin
    if Equal(one,den_one) then
      Clear(res);
      res := Create(natural32(0));
    else
      Div(res,den);
    end if;
    Clear(den); Clear(den_one);
    return res;
  end Inverse_Condition_Number;

  procedure Silent_Newton_Step
                ( n : in natural32; z : in out Vector; tol : in double_float;
                  err,rco,res : out Floating_Number; rank : out natural32 ) is

    y : Vector(1..integer32(n)) := f(z);
    ejm : Matrix(y'range,z'range) := jm(z);
    u : Matrix(y'range,y'range);
    v : Matrix(z'range,z'range);
    p : constant integer32 := z'length;
    m : constant integer32 := Min0(integer32(n)+1,p);
    e : Vector(1..p);
    s : Vector(1..m);
    info : integer32;
    dz : Vector(z'range);

  begin
    SVD(ejm,integer32(n),p,s,e,u,v,11,info);
    rco := Inverse_Condition_Number(s);
    rank := Numerical_Rank(s,tol);
    Min(y);
    dz := Solve(u,v,s,y);
    Clear(s); Clear(y);
    err := Max_Norm(dz);
    Add(z,dz);
    Clear(dz);
    y := f(z);
    res := Max_Norm(y);
    Clear(y);
  end Silent_Newton_Step;

  procedure Reporting_Newton_Step
                ( file : in file_type;
                  n : in natural32; z : in out Vector; tol : in double_float;
                  err,rco,res : out Floating_Number; rank : out natural32 ) is

    y : Vector(1..integer32(n)) := f(z);
    ejm : Matrix(y'range,z'range) := jm(z);
    u : Matrix(y'range,y'range);
    v : Matrix(z'range,z'range);
    p : constant integer32 := z'length;
    m : constant integer32 := Min0(integer32(n)+1,p);
    e : Vector(1..p);
    s : Vector(1..m);
    info : integer32;
    dz : Vector(z'range);
    rtmp : Floating_Number;

  begin
    SVD(ejm,integer32(n),p,s,e,u,v,11,info);
    put_line(file,"The singular values : ");
    put_line(file,s);
    rco := Inverse_Condition_Number(s);
    rank := Numerical_Rank(s,tol);
    Min(y);
    dz := Solve(u,v,s,y);
    Clear(s); Clear(y);
    err := Max_Norm(dz);
    Add(z,dz);
    Clear(dz);
    y := f(z);
    res := Max_Norm(y);
    Clear(y);
  end Reporting_Newton_Step;

  procedure Silent_Newton_Step_with_Singular_Values
                ( n : in natural32; z : in out Vector; tol : in double_float;
                  err,rco,res : out Floating_Number;
                  s : out Vector; rank : out natural32 ) is

    y : Vector(1..integer32(n)) := f(z);
    ejm : Matrix(y'range,z'range) := jm(z);
    u : Matrix(y'range,y'range);
    v : Matrix(z'range,z'range);
    p : constant integer32 := z'length;
   -- m : constant integer32 := Min0(integer32(n)+1,p);
    e : Vector(1..p);
    info : integer32;
    dz : Vector(z'range);
    rtmp : Floating_Number;

  begin
    SVD(ejm,integer32(n),p,s,e,u,v,11,info);
    rco := Inverse_Condition_Number(s);
    rank := Numerical_Rank(s,tol);
    Min(y);
    dz := Solve(u,v,s,y);
    Clear(y);
    err := Max_Norm(dz);
    Add(z,dz);
    Clear(dz);
    y := f(z);
    res := Max_Norm(y);
    Clear(y);
  end Silent_Newton_Step_with_Singular_Values;

  procedure Reporting_Newton_Step_with_Singular_Values
                ( file : in file_type;
                  n : in natural32; z : in out Vector; tol : in double_float;
                  err,rco,res : out Floating_Number;
                  s : out Vector; rank : out natural32 ) is

    y : Vector(1..integer32(n)) := f(z);
    ejm : Matrix(y'range,z'range) := jm(z);
    u : Matrix(y'range,y'range);
    v : Matrix(z'range,z'range);
    p : constant integer32 := z'length;
   -- m : constant integer32 := Min0(integer32(n)+1,p);
    e : Vector(1..p);
    info : integer32;
    dz : Vector(z'range);
    rtmp : Floating_Number;

  begin
    SVD(ejm,integer32(n),p,s,e,u,v,11,info);
    put_line(file,"The singular values : ");
    put_line(file,s);
    rco := Inverse_Condition_Number(s);
    rank := Numerical_Rank(s,tol);
    Min(y);
    dz := Solve(u,v,s,y);
    Clear(y);
    err := Max_Norm(dz);
    Add(z,dz);
    Clear(dz);
    y := f(z);
    res := Max_Norm(y);
    Clear(y);
  end Reporting_Newton_Step_with_Singular_Values;

end Multprec_Complex_Newton_Steps;
