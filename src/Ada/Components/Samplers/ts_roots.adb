with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Numbers_Polar;     use Standard_Complex_Numbers_Polar;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Polynomial_Roots;                   use Polynomial_Roots;

procedure ts_roots is

  tol : constant double_float := 1.0E-8;

  procedure Random_Affine_Solve
             ( n : in integer32; p : in Vector; fail : out boolean ) is

    s : Vector(1..n);
    r : double_float;

  begin
    put_line("The coefficients of a random polynomial : ");
    put_line(p);
    Affine_Solve(Standard_Output,p,s,r);
    fail := (r > tol);
  end Random_Affine_Solve;

  procedure Random_Projective_Solve
              ( n : in integer32; p : in Vector; fail : out boolean ) is

    s0,s1 : Vector(1..n);
    affres,prores : double_float;

  begin
    put_line("The coefficients of a random polynomial : ");
    put_line(p);
    Projective_Solve(Standard_Output,p,s0,s1,prores);
    fail := (prores > tol);
    for i in s1'range loop
       s1(i) := s1(i)/s0(i);
    end loop;
    put_line("Evaluation of p at affine roots :");
    Test_Affine_Roots(Standard_Output,p,s1,affres);
  end Random_Projective_Solve;

  procedure Main is

    n,m : integer32 := 0;
    affine_fail,projective_fail : boolean;

  begin
    new_line;
    put_line("Testing root finding of univariate polynomials.");
    new_line;
    put("Give the number of random tests : "); get(m);
    put("Give the degree : "); get(n);
    for i in 1..m loop
      declare
        p : Vector(0..n) := Random_Vector(0,n);
      begin
        Random_Affine_Solve(n,p,affine_fail);
        Random_Projective_Solve(n,p,projective_fail);
        if affine_fail then
          put("Failure in affine ");
          if projective_fail
           then put_line("and projective ");
          end if;
          put_line("solver.");
        end if;
      end;
      exit when projective_fail;
    end loop;
    if not projective_fail
     then put("Tested "); put(m,1); put_line(" cases successfully.");
    end if;
  end Main;

begin
  Main;
end ts_roots;
