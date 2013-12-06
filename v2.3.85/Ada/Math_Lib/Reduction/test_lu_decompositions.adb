with text_io;                            use text_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Floating_Linear_Solvers;
with Double_Double_Matrices_io;          use Double_Double_Matrices_io;
with Double_Double_Linear_Solvers;
with Quad_Double_Matrices_io;            use Quad_Double_Matrices_io;
with Quad_Double_Linear_Solvers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_Polar;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_Polar;
with DoblDobl_Complex_Matrices_io;       use DoblDobl_Complex_Matrices_io;
with DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_Polar;
with QuadDobl_Complex_Matrices_io;       use QuadDobl_Complex_Matrices_io;
with QuadDobl_Complex_Linear_Solvers;

package body Test_LU_Decompositions is

  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in Standard_Floating_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in double_float; output : in boolean;
                maxerr : out double_float; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Standard_Floating_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant Standard_Floating_Matrices.Matrix(1..n,1..n)
       := Standard_Floating_Linear_Solvers.Permute(P,A);
    L : Standard_Floating_Matrices.Matrix(1..n,1..n)
      := Standard_Floating_Linear_Solvers.Lower_Diagonal(LU);
    U : constant Standard_Floating_Matrices.Matrix(1..n,1..n)
      := Standard_Floating_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : Standard_Floating_Matrices.Matrix(1..n,1..n);
    err : double_float;
    use Standard_Floating_Matrices;

  begin
    Standard_Floating_Linear_Solvers.Permute_Lower(L,ipvt);
    LtimesU := L*U;
    if output then
      put("ipvt = "); put(ipvt); new_line;
      put_line("The permutation matrix P :"); put(P);
      put_line("The LU decomposition :"); put(LU);
      put_line("The multipliers L :"); put(L);
      put_line("The upper triangular U :"); put(U);
      put_line("L*U : "); put(LtimesU);
      put_line("P*A : "); put(PA);
    end if;
    maxerr := -1.0;
    for i in A'range(1) loop
      for j in A'range(2) loop
        err := abs(PA(i,j) - LtimesU(i,j));
        if maxerr < 0.0 then 
          maxerr := err;
        elsif err > maxerr then
          maxerr := err;
        end if;
      end loop;
    end loop;
    fail := (maxerr > tol);
  end Test_Decomposition;

  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in Double_Double_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in double_double; output : in boolean;
                maxerr : out double_double; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Double_Double_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant Double_Double_Matrices.Matrix(1..n,1..n)
       := Double_Double_Linear_Solvers.Permute(P,A);
    L : Double_Double_Matrices.Matrix(1..n,1..n)
      := Double_Double_Linear_Solvers.Lower_Diagonal(LU);
    U : constant Double_Double_Matrices.Matrix(1..n,1..n)
      := Double_Double_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : Double_Double_Matrices.Matrix(1..n,1..n);
    err : double_double;
    use Double_Double_Matrices;

  begin
    Double_Double_Linear_Solvers.Permute_Lower(L,ipvt);
    LtimesU := L*U;
    if output then
      put("ipvt = "); put(ipvt); new_line;
      put_line("The permutation matrix P :"); put(P);
      put_line("The LU decomposition :"); put(LU);
      put_line("The multipliers L :"); put(L);
      put_line("The upper triangular U :"); put(U);
      put_line("L*U : "); put(LtimesU);
      put_line("P*A : "); put(PA);
    end if;
    maxerr := create(-1.0);
    for i in A'range(1) loop
      for j in A'range(2) loop
        err := abs(PA(i,j) - LtimesU(i,j));
        if maxerr < 0.0 then 
          maxerr := err;
        elsif err > maxerr then
          maxerr := err;
        end if;
      end loop;
    end loop;
    fail := (maxerr > tol);
  end Test_Decomposition;

  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in Quad_Double_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in quad_double; output : in boolean;
                maxerr : out quad_double; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Quad_Double_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant Quad_Double_Matrices.Matrix(1..n,1..n)
       := Quad_Double_Linear_Solvers.Permute(P,A);
    L : Quad_Double_Matrices.Matrix(1..n,1..n)
      := Quad_Double_Linear_Solvers.Lower_Diagonal(LU);
    U : constant Quad_Double_Matrices.Matrix(1..n,1..n)
      := Quad_Double_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : Quad_Double_Matrices.Matrix(1..n,1..n);
    err : quad_double;
    use Quad_Double_Matrices;

  begin
    Quad_Double_Linear_Solvers.Permute_Lower(L,ipvt);
    LtimesU := L*U;
    if output then
      put("ipvt = "); put(ipvt); new_line;
      put_line("The permutation matrix P :"); put(P);
      put_line("The LU decomposition :"); put(LU);
      put_line("The multipliers L :"); put(L);
      put_line("The upper triangular U :"); put(U);
      put_line("L*U : "); put(LtimesU);
      put_line("P*A : "); put(PA);
    end if;
    maxerr := create(-1.0);
    for i in A'range(1) loop
      for j in A'range(2) loop
        err := abs(PA(i,j) - LtimesU(i,j));
        if maxerr < 0.0 then 
          maxerr := err;
        elsif err > maxerr then
          maxerr := err;
        end if;
      end loop;
    end loop;
    fail := (maxerr > tol);
  end Test_Decomposition;

  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in DoblDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in double_double; output : in boolean;
                maxerr : out double_double; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := DoblDobl_Complex_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
       := DoblDobl_Complex_Linear_Solvers.Permute(P,A);
    L : DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
      := DoblDobl_Complex_Linear_Solvers.Lower_Diagonal(LU);
    U : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
      := DoblDobl_Complex_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    dif : DoblDobl_Complex_Numbers.Complex_Number;
    err : double_double;
    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Numbers_Polar;
    use DoblDobl_Complex_Matrices;

  begin
    DoblDobl_Complex_Linear_Solvers.Permute_Lower(L,ipvt);
    LtimesU := L*U;
    if output then
      put("ipvt = "); put(ipvt); new_line;
      put_line("The permutation matrix P :"); put(P);
      put_line("The LU decomposition :"); put(LU);
      put_line("The multipliers L :"); put(L);
      put_line("The upper triangular U :"); put(U);
      put_line("L*U : "); put(LtimesU);
      put_line("P*A : "); put(PA);
    end if;
    maxerr := create(-1.0);
    for i in A'range(1) loop
      for j in A'range(2) loop
        dif := PA(i,j) - LtimesU(i,j);
        err := Radius(dif);
        if maxerr < 0.0 then 
          maxerr := err;
        elsif err > maxerr then
          maxerr := err;
        end if;
      end loop;
    end loop;
    fail := (maxerr > tol);
  end Test_Decomposition;

  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in Standard_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in double_float; output : in boolean;
                maxerr : out double_float; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Standard_Complex_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
       := Standard_Complex_Linear_Solvers.Permute(P,A);
    L : Standard_Complex_Matrices.Matrix(1..n,1..n)
      := Standard_Complex_Linear_Solvers.Lower_Diagonal(LU);
    U : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
      := Standard_Complex_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : Standard_Complex_Matrices.Matrix(1..n,1..n);
    dif : Standard_Complex_Numbers.Complex_Number;
    err : double_float;
    use Standard_Complex_Numbers;
    use Standard_Complex_Numbers_Polar;
    use Standard_Complex_Matrices;

  begin
    Standard_Complex_Linear_Solvers.Permute_Lower(L,ipvt);
    LtimesU := L*U;
    if output then
      put("ipvt = "); put(ipvt); new_line;
      put_line("The permutation matrix P :"); put(P);
      put_line("The LU decomposition :"); put(LU);
      put_line("The multipliers L :"); put(L);
      put_line("The upper triangular U :"); put(U);
      put_line("L*U : "); put(LtimesU);
      put_line("P*A : "); put(PA);
    end if;
    maxerr := -1.0;
    for i in A'range(1) loop
      for j in A'range(2) loop
        dif := PA(i,j) - LtimesU(i,j);
        err := Radius(dif);
        if maxerr < 0.0 then 
          maxerr := err;
        elsif err > maxerr then
          maxerr := err;
        end if;
      end loop;
    end loop;
    fail := (maxerr > tol);
  end Test_Decomposition;

  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in QuadDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in quad_double; output : in boolean;
                maxerr : out quad_double; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := QuadDobl_Complex_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
       := QuadDobl_Complex_Linear_Solvers.Permute(P,A);
    L : QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
      := QuadDobl_Complex_Linear_Solvers.Lower_Diagonal(LU);
    U : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
      := QuadDobl_Complex_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    dif : QuadDobl_Complex_Numbers.Complex_Number;
    err : quad_double;
    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Numbers_Polar;
    use QuadDobl_Complex_Matrices;

  begin
    QuadDobl_Complex_Linear_Solvers.Permute_Lower(L,ipvt);
    LtimesU := L*U;
    if output then
      put("ipvt = "); put(ipvt); new_line;
      put_line("The permutation matrix P :"); put(P);
      put_line("The LU decomposition :"); put(LU);
      put_line("The multipliers L :"); put(L);
      put_line("The upper triangular U :"); put(U);
      put_line("L*U : "); put(LtimesU);
      put_line("P*A : "); put(PA);
    end if;
    maxerr := create(-1.0);
    for i in A'range(1) loop
      for j in A'range(2) loop
        dif := PA(i,j) - LtimesU(i,j);
        err := Radius(dif);
        if maxerr < 0.0 then 
          maxerr := err;
        elsif err > maxerr then
          maxerr := err;
        end if;
      end loop;
    end loop;
    fail := (maxerr > tol);
  end Test_Decomposition;

end Test_LU_Decompositions;
