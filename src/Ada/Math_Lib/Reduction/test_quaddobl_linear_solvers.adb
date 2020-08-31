with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Quad_Double_Vectors_io;
with Quad_Double_Vector_Norms;
with QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vector_Norms;
with QuadDobl_Random_Vectors;
with QuadDobl_Random_Matrices;
with Quad_Double_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;
with Test_LU_Decompositions;

package body Test_QuadDobl_Linear_Solvers is

  procedure Run_Quad_Double_Linear_Solver
              ( A : in Quad_Double_Matrices.Matrix;
                b : in Quad_Double_Vectors.Vector ) is

    use Quad_Double_Vectors;
    use Quad_Double_Vectors_io;
    use Quad_Double_Matrices;
    use Quad_Double_Linear_Solvers;
    use Test_LU_Decompositions;

    AA : Matrix(A'range(1),A'range(2)) := A;
    n : constant integer32 := A'last(1);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    x : Vector(b'range) := b;
    res : Vector(b'range);
    tol : constant quad_double := create(1.0E-32);
    maxerr : quad_double;
    fail : boolean;

  begin
    new_line;
    lufac(AA,n,ipvt,info);
    put("info = "); put(info,1); new_line;
    lusolve(AA,n,ipvt,x);
    put_line("The solution :"); put_line(x);
    res := b - A*x;
    put_line("The residual :"); put_line(res);
    new_line;
    put_line("Testing the LU factorization ...");
    new_line;
    Test_Decomposition(n,A,AA,ipvt,tol,true,maxerr,fail);
    put("largest error : "); put(maxerr,3);
    put(" < "); put(tol,3); put(" ? ");
    if fail
     then put_line("bug!?");
     else put_line("okay.");
    end if;
  end Run_Quad_Double_Linear_Solver;

  procedure Run_QuadDobl_Complex_Linear_Solver
              ( A : in QuadDobl_Complex_Matrices.Matrix;
                b : in QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Vectors_io;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Linear_Solvers;
    use Test_LU_Decompositions;

    AA : Matrix(A'range(1),A'range(2)) := A;
    n : constant integer32 := A'last(1);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    x : Vector(b'range) := b;
    res : Vector(b'range);
    tol : constant quad_double := create(1.0E-32);
    maxerr : quad_double;
    fail : boolean;

  begin
    new_line;
    lufac(AA,n,ipvt,info);
    put("info = "); put(info,1); new_line;
    lusolve(AA,n,ipvt,x);
    put_line("The solution :"); put_line(x);
    res := b - A*x;
    put_line("The residual :"); put_line(res);
    new_line;
    put_line("Testing the LU factorization ...");
    new_line;
    Test_Decomposition(n,A,AA,ipvt,tol,true,maxerr,fail);
    put("largest error : "); put(maxerr,3);
    put(" < "); put(tol,3); put(" ? ");
    if fail
     then put_line("bug!?");
     else put_line("okay.");
    end if;
  end Run_QuadDobl_Complex_Linear_Solver;

  procedure Test_Quad_Double_Linear_Solver is

    use Quad_Double_Vectors;
    use Quad_Double_Matrices;

    n,m : integer32 := 0;

  begin
    new_line;
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      A : constant Matrix(1..n,1..m)
        := QuadDobl_Random_Matrices.Random_Matrix(natural32(n),natural32(m));
      b : constant Vector(1..n)
        := QuadDobl_Random_Vectors.Random_Vector(1,n);
    begin
      Run_Quad_Double_Linear_Solver(A,b);
    end;
  end Test_Quad_Double_Linear_Solver;

  procedure Test_QuadDobl_Complex_Linear_Solver is

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Matrices;

    n,m : integer32 := 0;

  begin
    new_line;
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      A : constant Matrix(1..n,1..m)
        := QuadDobl_Random_Matrices.Random_Matrix(natural32(n),natural32(m));
      b : constant Vector(1..n)
        := QuadDobl_Random_Vectors.Random_Vector(1,n);
    begin
      Run_QuadDobl_Complex_Linear_Solver(A,b);
    end;
  end Test_QuadDobl_Complex_Linear_Solver;

  procedure Test_Performance_Quad_Double_Linear_Solver
              ( d,k : in integer32; output,solve : in boolean ) is

    use Quad_Double_Vectors;
    use Quad_Double_Matrices;
    use Quad_Double_Linear_Solvers;

    A : constant Matrix(1..d,1..d)
      := QuadDobl_Random_Matrices.Random_Matrix(natural32(d),natural32(d));
    b : constant Vector(1..d) := QuadDobl_Random_Vectors.Random_Vector(1,d);
    ipvt : Standard_Integer_Vectors.Vector(1..d);
    info : integer32;
    AA : Matrix(1..d,1..d);
    x : Vector(1..d);

  begin
    for i in 1..k loop
      AA := A; x := b;
      lufac(AA,d,ipvt,info);
      if solve then
        lusolve(AA,d,ipvt,x);
        if output then
          declare
            res : Vector(1..d);
            nrm : quad_double;
          begin  
            res := b - A*x;
            nrm := Quad_Double_Vector_Norms.Max_Norm(res);
            put("Max norm of residual : "); put(nrm,3); new_line;
          end;
        end if;
      end if;
    end loop;
  end Test_Performance_Quad_Double_Linear_Solver;

  procedure Test_Performance_QuadDobl_Complex_Linear_Solver
              ( d,k : in integer32; output,solve : in boolean ) is

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Linear_Solvers;

    A : constant Matrix(1..d,1..d)
      := QuadDobl_Random_Matrices.Random_Matrix(natural32(d),natural32(d));
    b : constant Vector(1..d) := QuadDobl_Random_Vectors.Random_Vector(1,d);
    ipvt : Standard_Integer_Vectors.Vector(1..d);
    info : integer32;
    AA : Matrix(1..d,1..d);
    x : Vector(1..d);

  begin
    for i in 1..k loop
      AA := A; x := b;
      lufac(AA,d,ipvt,info);
      if solve then
        lusolve(AA,d,ipvt,x);
        if output then
          declare
            res : Vector(1..d);
            nrm : quad_double;
          begin  
            res := b - A*x;
            nrm := QuadDobl_Complex_Vector_Norms.Max_Norm(res);
            put("Max norm of residual : "); put(nrm,3); new_line;
          end;
        end if;
      end if;
    end loop;
  end Test_Performance_QuadDobl_Complex_Linear_Solver;

  procedure Quad_Double_Performance_Test is

    n,d : integer32 := 0;
    otp,outer,solve : boolean := false;
    ans : character;
    timer : Timing_Widget;

  begin
    new_line;
    put("Give number of times : "); get(n);
    put("Give the dimension : "); get(d);
    put("Do you want to solve (if not, LU only) ? (y/n) ");
    Ask_Yes_or_No(ans);
    solve := (ans = 'y');
    if solve then
      put("Do you want to the see the residual ? (y/n) ");
      Ask_Yes_or_No(ans);
      otp := (ans = 'y');
    end if;
    put("Generate new random matrix each time ? (y/n) ");
    Ask_Yes_or_No(ans);
    outer := (ans = 'y');
    tstart(timer);
    if outer then
      for i in 1..n loop
        Test_Performance_Quad_Double_Linear_Solver(d,1,otp,solve);
      end loop;
    else
      Test_Performance_Quad_Double_Linear_Solver(d,n,otp,solve);
    end if;
    tstop(timer);
    new_line;
    put("Solved "); put(n,1); put(" linear systems");
    put(" of dimension "); put(d,1); put_line("."); new_line;
    print_times(standard_output,timer,"performance test");
  end Quad_Double_Performance_Test;

  procedure QuadDobl_Complex_Performance_Test is

    n,d : integer32 := 0;
    otp,outer,solve : boolean := false;
    ans : character;
    timer : Timing_Widget;

  begin
    new_line;
    put("Give number of times : "); get(n);
    put("Give the dimension : "); get(d);
    put("Do you want to solve (if not, LU only) ? (y/n) ");
    Ask_Yes_or_No(ans);
    solve := (ans = 'y');
    if solve then
      put("Do you want to the see the residual ? (y/n) ");
      Ask_Yes_or_No(ans);
      otp := (ans = 'y');
    end if;
    put("Generate new random matrix each time ? (y/n) ");
    Ask_Yes_or_No(ans);
    outer := (ans = 'y');
    tstart(timer);
    if outer then
      for i in 1..n loop
        Test_Performance_QuadDobl_Complex_Linear_Solver(d,1,otp,solve);
      end loop;
    else
      Test_Performance_QuadDobl_Complex_Linear_Solver(d,n,otp,solve);
    end if;
    tstop(timer);
    new_line;
    put("Solved "); put(n,1); put(" linear systems");
    put(" of dimension "); put(d,1); put_line("."); new_line;
    print_times(standard_output,timer,"performance test");
  end QuadDobl_Complex_Performance_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of matrices of double doubles.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program;");
      put_line("  1. test linear solver on quad double numbers;");
      put_line("  2. test linear solver on complex quad double numbers;");
      put_line("  3. performance test on quad double arithmetic;");
      put_line("  4. performance test on quad double complex arithmetic.");
      put("Enter 0, 1, 2, 3, or 4 to make your choice : ");
      Ask_Alternative(ans,"0123456");
      exit when (ans = '0');
      case ans is
        when '1' => Test_Quad_Double_Linear_Solver;
        when '2' => Test_QuadDobl_Complex_Linear_Solver;
        when '3' => Quad_Double_Performance_Test;
        when '4' => QuadDobl_Complex_Performance_Test;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_QuadDobl_Linear_Solvers;
