with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_Singular_Values;  use Standard_Complex_Singular_Values;
with Symbol_Table;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Random_Polynomials;       use Standard_Random_Polynomials;
with Standard_Evaluation_Machine;       use Standard_Evaluation_Machine;
with Standard_Numerical_Derivatives;    use Standard_Numerical_Derivatives;

procedure ts_numdif is

-- DESCRIPTION :
--   Interactive development of numerical differentiation.

  procedure Compute_Derivatives
              ( n : in integer32; p : in Poly; x : in Vector ) is

  -- DESCRIPTION :
  --   Computes all n partial derivatives of the polynomial p,
  --   evaluated at the point x.

    h,error : double_float := 0.0;
    dp : Poly;
    num_dpx,sym_dpx : Complex_Number;

  begin
    Standard_Evaluation_Machine.Initialize(p);
    new_line; put("Give step : "); get(h);
    put("Derivatives with step h = "); put(h,3);
    put_line(" :");
    for i in 1..n loop
      dp := Diff(p,i);
      sym_dpx := Eval(dp,x);
      put("  formal : "); put(sym_dpx); new_line;
      num_dpx := Diff1(Standard_Evaluation_Machine.Evaluate'access,x,i,h);
      put("numeric1 : "); put(num_dpx); 
      error := AbsVal(sym_dpx-num_dpx);
      put("  error : "); put(error,3); new_line;
      num_dpx := Diff2(Standard_Evaluation_Machine.Evaluate'access,x,i,h);
      put("numeric2 : "); put(num_dpx); 
      error := AbsVal(sym_dpx-num_dpx);
      put("  error : "); put(error,3); new_line;
      num_dpx := Diff3(Standard_Evaluation_Machine.Evaluate'access,x,i,h);
      put("numeric3 : "); put(num_dpx); 
      error := AbsVal(sym_dpx-num_dpx);
      put("  error :"); put(error,3); new_line;
      Clear(dp);
    end loop;
    Standard_Evaluation_Machine.Clear;
  end Compute_Derivatives;

  procedure Compute_Derivatives
              ( n : in integer32; p : in Poly_Sys;
                ejm : in Eval_Jaco_Mat; x : in Vector ) is

    sym_dpx : constant Matrix(1..n,1..n) := Eval(ejm,x);
    num_dpx : Matrix(1..n,1..n);
    h,error : double_float := 0.0;

  begin
    Standard_Evaluation_Machine.Initialize(p);
    new_line; put("Give step : "); get(h);
    put("Derivatives with step h = "); put(h,3); new_line;
    num_dpx := Diff3(Standard_Evaluation_Machine.Evaluate'access,n,x,h);
    for i in 1..n loop
      for j in 1..n loop
        put(" formal("); put(i,1); put(",");
                         put(j,1); put("):"); put(sym_dpx(i,j)); new_line;
        put("numeric("); put(i,1); put(",");
                         put(j,1); put("):"); put(num_dpx(i,j)); 
        error := AbsVal(sym_dpx(i,j)-num_dpx(i,j));
        put("  error :"); put(error,3); new_line;
      end loop;
    end loop;
    Standard_Evaluation_Machine.Clear;
  end Compute_Derivatives;

  procedure SVD_Differences_Newton_Step
              ( file : in file_type; n : in integer32;
                x,fx : in out Vector; err,res : out double_float ) is

  -- DESCRIPTION :
  --   Does one step with Newton's method on a system with n equations,
  --   starting at the vector x.
  --   Divided differences (with extrapolation) are used to approximate
  --   the derivatives.  SVD is used to solve the linear systems.

  -- ON ENTRY :
  --   file    for intermediate output and diagnostics;
  --   n       the number of equations in the system;
  --   x       initial approximation for a root;
  --   fx      vector of range 1..n, function value at x.

  -- ON RETURN :
  --   x       improved approximation for a root;
  --   fx      vector of range 1..n, contains f(x);
  --   err     max norm of the correction term dx;
  --   res     max norm of the residual f(x).

  -- REQUIRED :
  --   The Standard_Evaluation_Machine is properly initialized
  --   or updated to return a vector of range 1..n each time the
  --   Evaluate function is called.

    h : constant double_float := 0.001;
    jp : Matrix(1..n,x'range);
    u : Matrix(1..n,1..n);
    v : Matrix(x'range,x'range);
    d : constant integer32 := x'length;
    m : constant integer32 := Min0(n+1,d);
    e : Vector(1..d);
    s : Vector(1..m);
    info : integer32;
    dx : Vector(x'range);

  begin
    jp := Diff3(Standard_Evaluation_Machine.Evaluate'access,n,x,h);
    put_line(file,"The Jacobi matrix : "); put(file,jp,3);
    SVD(jp,n,d,s,e,u,v,11,info);
    put_line(file,"The singular values :"); put_line(file,s);
    dx := Solve(u,v,s,-fx);
    err := Max_Norm(dx);
    Add(x,dx);
    put_line(file,"The new approximation : "); put_line(file,x);
    fx := Standard_Evaluation_Machine.Evaluate(x);
    put_line(file,"The function values : "); put_line(file,fx);
    res := Max_Norm(fx);
    put(file,"Error : "); put(file,err,3);
    put(file,"  and residual : "); put(file,res,3); new_line(file);
  end SVD_Differences_Newton_Step;

  procedure Call_Newton
              ( n : in integer32; p : in Poly_Sys; x : in out Vector ) is

    max : constant integer32 := 1;  -- max #deflations
    dim : integer32 := n;
    fp : Vector(1..n+max);
    err,res : double_float;
    ans : character;

  begin
    Standard_Evaluation_Machine.Initialize(p);
    fp(1..n) := Standard_Evaluation_Machine.Evaluate(x);
    loop
      SVD_Differences_Newton_Step(Standard_Output,dim,x,fp(1..dim),err,res);
      put("Do you want one more step ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when ans /= 'y';
      if dim - n < max then
        put("Do you want deflation ? ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          Standard_Evaluation_Machine.Deflate;
          dim := dim + 1;
          fp(dim) := Create(0.0);
        end if;
      end if;
    end loop;
    Standard_Evaluation_Machine.Clear;
  end Call_Newton;

  procedure Differentiate_at_Random_Point ( n : in integer32; p : in Poly ) is

  -- DESCRIPTION :
  --   Generates a random point and computes all partial derivatives.

    x : constant Vector(1..n) := Random_Vector(1,n);

  begin
    put_line("-> a random point : "); put_line(x);
    Compute_Derivatives(n,p,x);
  end Differentiate_at_Random_Point;

  procedure Differentiate_at_Random_Point
              ( n : in integer32; p : in Poly_Sys;
                ejm : in Eval_Jaco_Mat ) is

  -- DESCRIPTION :
  --   Generates a random point and computes the Jacobian matrix.

    x : constant Vector(1..n) := Random_Vector(1,n);

  begin
    put_line("-> a random point : "); put_line(x);
    Compute_Derivatives(n,p,ejm,x);
  end Differentiate_at_Random_Point;

  procedure Differentiate_at_Given_Point ( n : in integer32; p : in Poly ) is

  -- DESCRIPTION :
  --   Prompts the user for a point and computes all partial derivatives.

    x : Vector(1..n);

  begin
    put("Give "); put(n,1); put_line(" complex numbers : "); get(x);
    put_line("-> your point : "); put_line(x);
    Compute_Derivatives(n,p,x);
  end Differentiate_at_Given_Point;

  procedure Differentiate_at_Given_Point
              ( n : in integer32; p : in Poly_Sys;
                ejm : in Eval_Jaco_Mat ) is

  -- DESCRIPTION :
  --   Prompts the user for a point and computes the Jacobian matrix.

    x : Vector(1..n);

  begin
    put("Give "); put(n,1); put_line(" complex numbers : "); get(x);
    put_line("-> your point : "); put_line(x);
    Compute_Derivatives(n,p,ejm,x);
  end Differentiate_at_Given_Point;

  function Generate_Polynomial ( n : natural32 ) return Poly is

    d,m : natural32 := 0;
    p : Poly;
    ans : character;

  begin
    new_line;
    put_line("MENU to generate polynomial : ");
    put_line("  0. give your own polynomial;");
    put_line("  1. generate random sparse polynomial;");
    put_line("  2. generate random dense polynomial;");
    put("Type 0, 1, or 2 to choose : ");
    Ask_Alternative(ans,"012");
    if ans = '0' then
      put("Give a polynomial : "); get(p);
      put("-> your polynomial : "); put(p); new_line;
    else
      put("Give degree : "); get(d);
      if ans = '2'
       then p := Random_Dense_Poly(n,d,0);
       else put("Give number of terms : "); get(m);
            p := Random_Sparse_Poly(n,d,m,0);
      end if;
      put("-> p = "); put_line(p);
    end if;
    return p;
  end Generate_Polynomial;

  function Generate_System return Poly_Sys is

    ans : character;
    n : natural32 := 0;

  begin
    new_line;
    put_line("MENU to generate polynomial system :");
    put_line("  0. give your own polynomial system;");
    put_line("  1. generate random sparse polynomial system;");
    put_line("  2. generate random dense polynomial system;");
    put("Type 0, 1, or 2 to choose : "); 
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' =>
        declare
          lp : Link_to_Poly_Sys;
        begin
          get(lp); 
          return lp.all;
        end;
      when others =>
        put("Give dimension of the system : "); get(n);
        declare 
          p : Poly_Sys(1..integer32(n));
          d,m : natural32 := 0;
        begin
          put("Give degree of each equation : "); get(d);
          if ans = '2' then
            for i in p'range loop
              p(i) := Random_Dense_Poly(n,d,0);
            end loop;
          else
            put("Give number of terms : "); get(m);
            for i in p'range loop
              p(i) := Random_Sparse_Poly(n,d,m,0);
            end loop;
          end if;
          return p;
        end;
    end case;
  end Generate_System;

  procedure Differentiate_One_Polynomial is

    n : natural32 := 0;
    p : Poly;
    ans : character;

  begin
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    p := Generate_Polynomial(n);
    loop
      new_line;
      put_line("MENU to test numerical differentiation :");
      put_line("  0. exit this program;");
      put_line("  1. differentiate at random point;");
      put_line("  2. differentiate at a given point;");
      put("Type 0, 1, or 2 to choose : ");
      Ask_Alternative(ans,"012");
      exit when (ans = '0');
      if ans = '1'
       then Differentiate_at_Random_Point(integer32(n),p);
       else Differentiate_at_Given_Point(integer32(n),p);
      end if;
    end loop;
  end Differentiate_One_Polynomial;

  procedure Differentiate_Polynomial_System is

    p : constant Poly_Sys := Generate_System;
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    ejm : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
    ans : character;

  begin
    put_line("The polynomial system : "); put(p);
    loop
      new_line;
      put_line("MENU to test numerical differentiation :");
      put_line("  0. exit this program;");
      put_line("  1. differentiate at random point;");
      put_line("  2. differentiate at a given point;");
      put("Type 0, 1, or 2 to choose : ");
      Ask_Alternative(ans,"012");
      exit when (ans = '0');
      if ans = '1'
       then Differentiate_at_Random_Point(p'last,p,ejm);
       else Differentiate_at_Given_Point(p'last,p,ejm);
      end if;
    end loop;
  end Differentiate_Polynomial_System;

  procedure Test_Newton is

    p : constant Poly_Sys := Generate_System;
    n : constant integer32 := p'last;
    x : Vector(1..n);

  begin
    put_line("The polynomial system : "); put(p);
    put("Give "); put(n,1);
    put_line(" complex numbers as start point : ");
    get(x);
    Call_Newton(n,p,x);
  end Test_Newton;

  procedure Main is

    ans : character;

  begin
    put_line("Choose one of the following : ");
    put_line("  1. numerical differentiation of one polynomial;");
    put_line("  2. numerical differentiation of polynomial system;");
    put_line("  3. Newton's method on polynomial system.");
    put("Type 1, 2, or 3 to choose : ");
    Ask_Alternative(ans,"123");
    case ans is
      when '1' => Differentiate_One_Polynomial;
      when '2' => Differentiate_Polynomial_System;
      when others => Test_Newton;
    end case;
  end Main;

begin
  new_line;
  put_line("Numerical differentiation of polynomials...");
  new_line;
  Main;
end ts_numdif;
