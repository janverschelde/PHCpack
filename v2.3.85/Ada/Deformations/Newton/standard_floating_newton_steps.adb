with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_Singular_Values;  use Standard_Complex_Singular_Values;
with Standard_Numerical_Rank;           use Standard_Numerical_Rank;

package body Standard_Floating_Newton_Steps is

  procedure Silent_Newton_Step
                ( n : in natural32; z : in out Vector; tol : double_float;
                  err,rco,res : out double_float; rank : out natural32 ) is

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
    SVD(ejm,n,p,s,e,u,v,11,info);
    rco := REAL_PART(s(s'last))/REAL_PART(s(s'first));
    rank := Numerical_Rank(s,tol);
    dz := Solve(u,v,s,-y);
    err := Max_Norm(dz);
    Add(z,dz);
    y := f(z);
    res := Max_Norm(y);
  end Silent_Newton_Step;

  procedure Reporting_Newton_Step
                ( file : in file_type;
                  n : in natural32; z : in out Vector; tol : in double_float;
                  err,rco,res : out double_float; rank : out natural32 ) is

    y : Vector(1..integer32(n)) := f(z);
    ejm : Matrix(y'range,z'range) := jm(z);
    u : Matrix(y'range,y'range);
    v : Matrix(z'range,z'range);
    p : constant natural := z'length;
    m : constant natural := Min0(integer32(n)+1,p);
    e : Vector(1..p);
    s : Vector(1..m);
    info : integer32;
    dz : Vector(z'range);

  begin
    SVD(ejm,n,p,s,e,u,v,11,info);
    put_line(file,"The singular values : ");
    put_line(file,s);
    rco := REAL_PART(s(s'last))/REAL_PART(s(s'first));
    rank := Numerical_Rank(s,tol);
    dz := Solve(u,v,s,-y);
    err := Max_Norm(dz);
    Add(z,dz);
    y := f(z);
    res := Max_Norm(y);
  end Reporting_Newton_Step;

  procedure Silent_Newton_Step_with_Singular_Values
                ( n : in natural32; z : in out Vector; tol : in double_float;
                  err,rco,res : out double_float;
                  s : out Vector; rank : out natural32 ) is

    y : Vector(1..integer32(n)) := f(z);
    ejm : Matrix(y'range,z'range) := jm(z);
    u : Matrix(y'range,y'range);
    v : Matrix(z'range,z'range);
    p : constant integer32 := z'length;
    m : constant integer32 := Min0(integer32(n)+1,p);
    e : Vector(1..p);
    info : integer;
    dz : Vector(z'range);

  begin
    SVD(ejm,n,p,s,e,u,v,11,info);
    rco := REAL_PART(s(s'last))/REAL_PART(s(s'first));
    rank := Numerical_Rank(s,tol);
    dz := Solve(u,v,s,-y);
    err := Max_Norm(dz);
    Add(z,dz);
    y := f(z);
    res := Max_Norm(y);
  end Silent_Newton_Step_with_Singular_Values;

  procedure Reporting_Newton_Step_with_Singular_Values
                ( file : in file_type;
                  n : in natural32; z : in out Vector; tol : in double_float;
                  err,rco,res : out double_float;
                  s : out Vector; rank : out natural32 ) is

    y : Vector(1..integer32(n)) := f(z);
    ejm : Matrix(y'range,z'range) := jm(z);
    u : Matrix(y'range,y'range);
    v : Matrix(z'range,z'range);
    p : constant integer32 := z'length;
    m : constant integer32 := Min0(integer32(n)+1,p);
    e : Vector(1..p);
    info : integer;
    dz : Vector(z'range);

  begin
    SVD(ejm,n,p,s,e,u,v,11,info);
    put_line(file,"The singular values : ");
    put_line(file,s);
    rco := REAL_PART(s(s'last))/REAL_PART(s(s'first));
    rank := Numerical_Rank(s,tol);
    dz := Solve(u,v,s,-y);
    err := Max_Norm(dz);
    Add(z,dz);
    y := f(z);
    res := Max_Norm(y);
--  exception
--    when others
--      => put_line("exception caught by "
--           & "Reporting_Newton_Step_with_Singular_Values");
--         raise;
  end Reporting_Newton_Step_with_Singular_Values;

end Standard_Floating_Newton_Steps;
