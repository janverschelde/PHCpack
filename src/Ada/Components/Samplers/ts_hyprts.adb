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
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Random_Polynomials;       use Standard_Random_Polynomials;
with Polynomial_Roots;                  use Polynomial_Roots;
with Hypersurface_Roots;                use Hypersurface_Roots;

procedure ts_hyprts is

  tol : constant double_float := 1.0E-8;

  procedure Random_Poly ( n : out natural32; p : out Poly ) is

  -- DESCRIPTION :
  --   Interactive generation of a random polynomial.

    d,m : natural32 := 0;

  begin
    put_line("Generating a random multivariate complex polynomial...");
    n := 0;
    put("  Give the number of variables : "); get(n);
    put("  Give the degree : "); get(d);
    put("  Give the number of terms : "); get(m);
    p := Random_Sparse_Poly(n,d,m,0);
    put_line("The random polynomial : ");
    put_line(p);
  end Random_Poly;

  procedure Test_Substitution ( n : in integer32; p : in Poly ) is

  -- DESCRIPTION :
  --   Compares the evaluation of p at t*v with the substitution.

    v : constant Vector(1..n) := Random_Vector(1,n);
    t : constant Complex_Number := Random1;
    c : constant Vector(0..Degree(p)) := Substitute(p,v);
    x : constant Vector(1..n) := t*v;
    y1 : constant Complex_Number := Eval(p,x);
    y2 : constant Complex_Number := Eval(c,t);

  begin
   -- put_line("Substituted coefficient vector : "); put_line(c);
    put("Evaluation at random line  : "); put(y1); new_line;
    put("Evaluation by substitution : "); put(y2); new_line;
  end Test_Substitution;

  procedure Test_Affine_Solving
               ( n : in integer32; p : in Poly; fail : out boolean ) is

  -- DESCRIPTION :
  --   Solves p(t*v) = 0 after explicit substitution.

    v : constant Vector(1..n) := Random_Vector(1,n);
    c : constant Vector(0..Degree(p)) := Substitute(p,v);
    s : Vector(1..c'last);
    res : double_float;

  begin
    Affine_Solve(Standard_Output,c,s,res);
    fail := (res > tol);
  end Test_Affine_Solving;

  procedure Test_Projective_Solving
               ( n : in integer32; p : in Poly; fail : out boolean ) is

  -- DESCRIPTION :
  --   Solves p(t*v) = 0 after explicit substitution.

    v : constant Vector(1..n) := Random_Vector(1,n);
    c : constant Vector(0..Degree(p)) := Substitute(p,v);
    s0,s1 : Vector(1..c'last);
    res,abseva : double_float;
    x : Vector(1..n);
    t,y : Complex_Number;

  begin
    Projective_Solve(Standard_Output,c,s0,s1,res);
    fail := (res > tol);
    res := 0.0;
    for i in s1'range loop
      t := s1(i)/s0(i);
      put(i,3); put(" : ");
      put(t); put(" : ");
      x := t*v;
      y := Eval(p,x);
      abseva := AbsVal(y);
      put(abseva,3); new_line;
      if abseva > res
       then res := abseva;
      end if;
    end loop;
    put("Maximal residual at affine roots : ");
    put(res,3); new_line;
  end Test_Projective_Solving;

  procedure Test_Affine_Path_Tracker
                ( n : in integer32; p : in Poly; fail : out boolean ) is

    v0 : constant Vector(1..n) := Random_Vector(1,n);
    v1 : constant Vector(1..n) := Random_Vector(1,n);
    c : Vector(0..Degree(p)) := Substitute(p,v0);
    s : Vector(1..c'last);
    res : double_float;

  begin
    Affine_Solve(Standard_Output,c,s,res);
    fail := (res > tol);
    if not fail then
      Affine_Track_Moving_Line(p,v0,v1,s);
      c := Substitute(p,v1);
      Test_Affine_Roots(Standard_Output,c,s,res);
      fail := (res > tol);
    end if;
  end Test_Affine_Path_Tracker;

  procedure Test_Projective_Path_Tracker
                ( n : in integer32; p : in Poly; fail : out boolean ) is

    v0 : constant Vector(1..n) := Random_Vector(1,n);
    v1 : constant Vector(1..n) := Random_Vector(1,n);
    c : Vector(0..Degree(p)) := Substitute(p,v0);
    s0,s1 : Vector(1..c'last);
    res : double_float;

  begin
    Projective_Solve(Standard_Output,c,s0,s1,res);
    fail := (res > tol);
    if not fail then
      Projective_Track_Moving_Line(p,v0,v1,s0,s1);
      c := Substitute(p,v1);
      Test_Projective_Roots(Standard_Output,c,s0,s1,res);
      fail := (res > tol);
    end if;
  end Test_Projective_Path_Tracker;

  procedure Find_Permutation ( x1,x2 : in Vector ) is

    absdif : double_float;

  begin
    for i in x1'range loop
      put("Path "); put(i,1); put(" -> ");
      for j in x2'range loop
        absdif := AbsVal(x1(i) - x2(j));
        if absdif < tol
         then put(j,1); put_line("."); exit;
        end if;
      end loop;
    end loop;
  end Find_Permutation;

  procedure Test_Monodromy_Loop
                ( n : in integer32; p : in Poly; fail : out boolean ) is

  -- DESCRIPTION :
  --   Tests one loop with the monodromy: v0 -> v1 -> v2 -> v0.

    v0 : constant Vector(1..n) := Random_Vector(1,n);
    v1 : constant Vector(1..n) := Random_Vector(1,n);
    v2 : constant Vector(1..n) := Random_Vector(1,n);
    c : constant Vector(0..Degree(p)) := Substitute(p,v0);
    sA,sB : Vector(1..c'last);
    res : double_float;

  begin
    Affine_Solve(Standard_Output,c,sA,res);
    sB := sA;
    fail := (res > tol);
    if not fail then
      Affine_Track_Moving_Line(p,v0,v1,sB);
      Affine_Track_Moving_Line(p,v1,v2,sB);
      Affine_Track_Moving_Line(p,v2,v0,sB);
      Test_Affine_Roots(Standard_Output,c,sB,res);
      fail := (res > tol);
      Find_Permutation(sA,sB);
    end if;
  end Test_Monodromy_Loop;

  procedure Random_Tester ( nb,k : in integer32 ) is

  -- DESCRIPTION :
  --   Performs as many random tests as the number nb.
  --   Tests the solver, affine path tracker, projective path tracker,
  --   or monodromy loop, depending whether k = 1, 2, 3, or 4.

    n,d,m : integer32 := 0;
    p : Poly;
    fail : boolean;

  begin
    put_line("Generation of random polynomials...");
    put("  Give the number of variables : "); get(n);
    put("  Give the degree : "); get(d);
    put("  Give the number of terms : "); get(m);
    for i in 1..nb loop
      p := Random_Sparse_Poly(natural32(n),natural32(d),natural32(m),0);
      put_line(p);
      case k is
        when 1 => Test_Affine_Solving(n,p,fail);
        when 2 => Test_Projective_Solving(n,p,fail);
        when 3 => Test_Affine_Path_Tracker(n,p,fail);
        when 4 => Test_Projective_Path_Tracker(n,p,fail);
        when 5 => Test_Monodromy_Loop(n,p,fail);
        when others => put_line("Invalid value for k."); return;
      end case;
      exit when fail;
    end loop;
    if fail
     then put_line("Ended in failure.");
     else put("Tested "); put(nb,1); put_line(" cases successfully.");
    end if;
  end Random_Tester;

  procedure Main is

    n : natural32 := 0;
    p : Poly;
    ans : character;

  begin
    new_line;
    put_line("Tracing points on a hypersurface.");
    new_line;
    put_line("Choose one of the following :");
    put_line("  1. test substitution of random vector;");
    put_line("  2. random test on solving in affine space;");
    put_line("  3. random test on solving in projective space;");
    put_line("  4. move random line and solutions in affine space;");
    put_line("  5. move random line and solutions in projective space;");
    put_line("  6. execute monodromy loop on random polynomial.");
    put("Type 1, 2, 3, 4, 5 or 6 to choose : ");
    Ask_Alternative(ans,"123456");
    new_line;
    case ans is
      when '1' => Random_Poly(n,p);
                  Test_Substitution(integer32(n),p);
      when others =>
             put("Give the number of random tests : "); get(n);
             case ans is
               when '2' => Random_Tester(integer32(n),1);
               when '3' => Random_Tester(integer32(n),2);
               when '4' => Random_Tester(integer32(n),3);
               when '5' => Random_Tester(integer32(n),4);
               when '6' => Random_Tester(integer32(n),5);
               when others => put_line("invalid option.");
             end case;
    end case;
  end Main;

begin
  Main;
end ts_hyprts;
