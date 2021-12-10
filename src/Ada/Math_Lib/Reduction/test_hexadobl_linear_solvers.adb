with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with Hexa_Double_Numbers_io;             use Hexa_Double_Numbers_io;
with Standard_Integer_Vectors;
with Hexa_Double_Vectors_io;
with HexaDobl_Complex_Vectors_io;
with HexaDobl_Random_Vectors;
with HexaDobl_Random_Matrices;
with Hexa_Double_Linear_Solvers;
with HexaDobl_Complex_Linear_Solvers;
with Test_LU_Decompositions;

package body Test_HexaDobl_Linear_Solvers is

  procedure Run_Hexa_Double_Linear_Solver
              ( A : in Hexa_Double_Matrices.Matrix;
                b : in Hexa_Double_Vectors.Vector ) is

    use Hexa_Double_Vectors;
    use Hexa_Double_Vectors_io;
    use Hexa_Double_Matrices;
    use Hexa_Double_Linear_Solvers;
    use Test_LU_Decompositions;

    AA : Matrix(A'range(1),A'range(2)) := A;
    n : constant integer32 := A'last(1);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    x : Vector(b'range) := b;
    res : Vector(b'range);
    tol : constant hexa_double := create(1.0E-16);
    maxerr : hexa_double;
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
  end Run_Hexa_Double_Linear_Solver;

  procedure Run_HexaDobl_Complex_Linear_Solver
              ( A : in HexaDobl_Complex_Matrices.Matrix;
                b : in HexaDobl_Complex_Vectors.Vector ) is

    use HexaDobl_Complex_Vectors;
    use HexaDobl_Complex_Vectors_io;
    use HexaDobl_Complex_Matrices;
    use HexaDobl_Complex_Linear_Solvers;
    use Test_LU_Decompositions;

    AA : Matrix(A'range(1),A'range(2)) := A;
    n : constant integer32 := A'last(1);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    x : Vector(b'range) := b;
    res : Vector(b'range);
    tol : constant hexa_double := create(1.0E-16);
    maxerr : hexa_double;
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
  end Run_HexaDobl_Complex_Linear_Solver;

  procedure Test_Hexa_Double_Linear_Solver is

    use Hexa_Double_Vectors;
    use Hexa_Double_Matrices;

    n,m : integer32 := 0;

  begin
    new_line;
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      A : constant Matrix(1..n,1..m)
        := HexaDobl_Random_Matrices.Random_Matrix(natural32(n),natural32(m));
      b : constant Vector(1..n)
        := HexaDobl_Random_Vectors.Random_Vector(1,n);
    begin
      Run_Hexa_Double_Linear_Solver(A,b);
    end;
  end Test_Hexa_Double_Linear_Solver;

  procedure Test_HexaDobl_Complex_Linear_Solver is

    use HexaDobl_Complex_Vectors;
    use HexaDobl_Complex_Matrices;

    n,m : integer32 := 0;

  begin
    new_line;
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      A : constant Matrix(1..n,1..m)
        := HexaDobl_Random_Matrices.Random_Matrix(natural32(n),natural32(m));
      b : constant Vector(1..n)
        := HexaDobl_Random_Vectors.Random_Vector(1,n);
    begin
      Run_HexaDobl_Complex_Linear_Solver(A,b);
    end;
  end Test_HexaDobl_Complex_Linear_Solver;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Interactive testing of hexa double LU factorization.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program;");
      put_line("  1. test hexa double linear solver;");
      put_line("  2. test hexa double complex linear solver.");
      put("Type 1 or 2 to select a test, or 0 to exit : ");
      Ask_Alternative(ans,"012");
      exit when (ans = '0');
      case ans is
        when '1' => Test_Hexa_Double_Linear_Solver;
        when '2' => Test_HexaDobl_Complex_Linear_Solver;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_HexaDobl_Linear_Solvers;
