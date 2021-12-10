with text_io;                            use text_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Floating_Linear_Solvers;
with Double_Double_Matrices_io;          use Double_Double_Matrices_io;
with Double_Double_Linear_Solvers;
with Triple_Double_Matrices_io;          use Triple_Double_Matrices_io;
with Triple_Double_Linear_Solvers;
with Quad_Double_Matrices_io;            use Quad_Double_Matrices_io;
with Quad_Double_Linear_Solvers;
with Penta_Double_Matrices_io;           use Penta_Double_Matrices_io;
with Penta_Double_Linear_Solvers;
with Octo_Double_Matrices_io;            use Octo_Double_Matrices_io;
with Octo_Double_Linear_Solvers;
with Deca_Double_Matrices_io;            use Deca_Double_Matrices_io;
with Deca_Double_Linear_Solvers;
with Hexa_Double_Matrices_io;            use Hexa_Double_Matrices_io;
with Hexa_Double_Linear_Solvers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_Polar;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_Polar;
with DoblDobl_Complex_Matrices_io;       use DoblDobl_Complex_Matrices_io;
with DoblDobl_Complex_Linear_Solvers;
with TripDobl_Complex_Numbers;
with TripDobl_Complex_Matrices_io;       use TripDobl_Complex_Matrices_io;
with TripDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_Polar;
with QuadDobl_Complex_Matrices_io;       use QuadDobl_Complex_Matrices_io;
with QuadDobl_Complex_Linear_Solvers;
with PentDobl_Complex_Numbers;
with PentDobl_Complex_Matrices_io;       use PentDobl_Complex_Matrices_io;
with PentDobl_Complex_Linear_Solvers;
with OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Matrices_io;       use OctoDobl_Complex_Matrices_io;
with OctoDobl_Complex_Linear_Solvers;
with DecaDobl_Complex_Numbers;
with DecaDobl_Complex_Matrices_io;       use DecaDobl_Complex_Matrices_io;
with DecaDobl_Complex_Linear_Solvers;
with HexaDobl_Complex_Numbers;
with HexaDobl_Complex_Matrices_io;       use HexaDobl_Complex_Matrices_io;
with HexaDobl_Complex_Linear_Solvers;

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
                A,LU : in Triple_Double_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in triple_double; output : in boolean;
                maxerr : out triple_double; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Triple_Double_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant Triple_Double_Matrices.Matrix(1..n,1..n)
       := Triple_Double_Linear_Solvers.Permute(P,A);
    L : Triple_Double_Matrices.Matrix(1..n,1..n)
      := Triple_Double_Linear_Solvers.Lower_Diagonal(LU);
    U : constant Triple_Double_Matrices.Matrix(1..n,1..n)
      := Triple_Double_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : Triple_Double_Matrices.Matrix(1..n,1..n);
    err : triple_double;
    use Triple_Double_Matrices;

  begin
    Triple_Double_Linear_Solvers.Permute_Lower(L,ipvt);
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
                A,LU : in Penta_Double_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in penta_double; output : in boolean;
                maxerr : out penta_double; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Penta_Double_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant Penta_Double_Matrices.Matrix(1..n,1..n)
       := Penta_Double_Linear_Solvers.Permute(P,A);
    L : Penta_Double_Matrices.Matrix(1..n,1..n)
      := Penta_Double_Linear_Solvers.Lower_Diagonal(LU);
    U : constant Penta_Double_Matrices.Matrix(1..n,1..n)
      := Penta_Double_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : Penta_Double_Matrices.Matrix(1..n,1..n);
    err : penta_double;
    use Penta_Double_Matrices;

  begin
    Penta_Double_Linear_Solvers.Permute_Lower(L,ipvt);
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
                A,LU : in Octo_Double_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in octo_double; output : in boolean;
                maxerr : out octo_double; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Octo_Double_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant Octo_Double_Matrices.Matrix(1..n,1..n)
       := Octo_Double_Linear_Solvers.Permute(P,A);
    L : Octo_Double_Matrices.Matrix(1..n,1..n)
      := Octo_Double_Linear_Solvers.Lower_Diagonal(LU);
    U : constant Octo_Double_Matrices.Matrix(1..n,1..n)
      := Octo_Double_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : Octo_Double_Matrices.Matrix(1..n,1..n);
    err : octo_double;
    use Octo_Double_Matrices;

  begin
    Octo_Double_Linear_Solvers.Permute_Lower(L,ipvt);
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
                A,LU : in Deca_Double_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in deca_double; output : in boolean;
                maxerr : out deca_double; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Deca_Double_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant Deca_Double_Matrices.Matrix(1..n,1..n)
       := Deca_Double_Linear_Solvers.Permute(P,A);
    L : Deca_Double_Matrices.Matrix(1..n,1..n)
      := Deca_Double_Linear_Solvers.Lower_Diagonal(LU);
    U : constant Deca_Double_Matrices.Matrix(1..n,1..n)
      := Deca_Double_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : Deca_Double_Matrices.Matrix(1..n,1..n);
    err : deca_double;
    use Deca_Double_Matrices;

  begin
    Deca_Double_Linear_Solvers.Permute_Lower(L,ipvt);
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
                A,LU : in Hexa_Double_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in hexa_double; output : in boolean;
                maxerr : out hexa_double; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := Hexa_Double_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant Hexa_Double_Matrices.Matrix(1..n,1..n)
       := Hexa_Double_Linear_Solvers.Permute(P,A);
    L : Hexa_Double_Matrices.Matrix(1..n,1..n)
      := Hexa_Double_Linear_Solvers.Lower_Diagonal(LU);
    U : constant Hexa_Double_Matrices.Matrix(1..n,1..n)
      := Hexa_Double_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : Hexa_Double_Matrices.Matrix(1..n,1..n);
    err : hexa_double;
    use Hexa_Double_Matrices;

  begin
    Hexa_Double_Linear_Solvers.Permute_Lower(L,ipvt);
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
                A,LU : in TripDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in triple_double; output : in boolean;
                maxerr : out triple_double; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := TripDobl_Complex_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant TripDobl_Complex_Matrices.Matrix(1..n,1..n)
       := TripDobl_Complex_Linear_Solvers.Permute(P,A);
    L : TripDobl_Complex_Matrices.Matrix(1..n,1..n)
      := TripDobl_Complex_Linear_Solvers.Lower_Diagonal(LU);
    U : constant TripDobl_Complex_Matrices.Matrix(1..n,1..n)
      := TripDobl_Complex_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : TripDobl_Complex_Matrices.Matrix(1..n,1..n);
    dif : TripDobl_Complex_Numbers.Complex_Number;
    err : triple_double;
    use TripDobl_Complex_Numbers;
   -- use TripDobl_Complex_Numbers_Polar;
    use TripDobl_Complex_Matrices;

  begin
    TripDobl_Complex_Linear_Solvers.Permute_Lower(L,ipvt);
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
        err := AbsVal(dif); -- Radius(dif);
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

  procedure Test_Decomposition
              ( n : in integer32;
                A,LU : in PentDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in penta_double; output : in boolean;
                maxerr : out penta_double; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := PentDobl_Complex_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant PentDobl_Complex_Matrices.Matrix(1..n,1..n)
       := PentDobl_Complex_Linear_Solvers.Permute(P,A);
    L : PentDobl_Complex_Matrices.Matrix(1..n,1..n)
      := PentDobl_Complex_Linear_Solvers.Lower_Diagonal(LU);
    U : constant PentDobl_Complex_Matrices.Matrix(1..n,1..n)
      := PentDobl_Complex_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : PentDobl_Complex_Matrices.Matrix(1..n,1..n);
    dif : PentDobl_Complex_Numbers.Complex_Number;
    err : penta_double;
    use PentDobl_Complex_Numbers;
   -- use PentDobl_Complex_Numbers_Polar;
    use PentDobl_Complex_Matrices;

  begin
    PentDobl_Complex_Linear_Solvers.Permute_Lower(L,ipvt);
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
        err := AbsVal(dif); -- Radius(dif);
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
                A,LU : in OctoDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in octo_double; output : in boolean;
                maxerr : out octo_double; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := OctoDobl_Complex_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant OctoDobl_Complex_Matrices.Matrix(1..n,1..n)
       := OctoDobl_Complex_Linear_Solvers.Permute(P,A);
    L : OctoDobl_Complex_Matrices.Matrix(1..n,1..n)
      := OctoDobl_Complex_Linear_Solvers.Lower_Diagonal(LU);
    U : constant OctoDobl_Complex_Matrices.Matrix(1..n,1..n)
      := OctoDobl_Complex_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : OctoDobl_Complex_Matrices.Matrix(1..n,1..n);
    dif : OctoDobl_Complex_Numbers.Complex_Number;
    err : octo_double;
    use OctoDobl_Complex_Numbers;
   -- use OctoDobl_Complex_Numbers_Polar;
    use OctoDobl_Complex_Matrices;

  begin
    OctoDobl_Complex_Linear_Solvers.Permute_Lower(L,ipvt);
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
        err := AbsVal(dif); -- Radius(dif);
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
                A,LU : in DecaDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in deca_double; output : in boolean;
                maxerr : out deca_double; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := DecaDobl_Complex_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant DecaDobl_Complex_Matrices.Matrix(1..n,1..n)
       := DecaDobl_Complex_Linear_Solvers.Permute(P,A);
    L : DecaDobl_Complex_Matrices.Matrix(1..n,1..n)
      := DecaDobl_Complex_Linear_Solvers.Lower_Diagonal(LU);
    U : constant DecaDobl_Complex_Matrices.Matrix(1..n,1..n)
      := DecaDobl_Complex_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : DecaDobl_Complex_Matrices.Matrix(1..n,1..n);
    dif : DecaDobl_Complex_Numbers.Complex_Number;
    err : deca_double;
    use DecaDobl_Complex_Numbers;
   -- use DecaDobl_Complex_Numbers_Polar;
    use DecaDobl_Complex_Matrices;

  begin
    DecaDobl_Complex_Linear_Solvers.Permute_Lower(L,ipvt);
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
        err := AbsVal(dif); -- Radius(dif);
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
                A,LU : in HexaDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                tol : in hexa_double; output : in boolean;
                maxerr : out hexa_double; fail : out boolean ) is

    P : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
      := HexaDobl_Complex_Linear_Solvers.Permutation_Matrix(ipvt);
    PA : constant HexaDobl_Complex_Matrices.Matrix(1..n,1..n)
       := HexaDobl_Complex_Linear_Solvers.Permute(P,A);
    L : HexaDobl_Complex_Matrices.Matrix(1..n,1..n)
      := HexaDobl_Complex_Linear_Solvers.Lower_Diagonal(LU);
    U : constant HexaDobl_Complex_Matrices.Matrix(1..n,1..n)
      := HexaDobl_Complex_Linear_Solvers.Upper_Diagonal(LU);
    LtimesU : HexaDobl_Complex_Matrices.Matrix(1..n,1..n);
    dif : HexaDobl_Complex_Numbers.Complex_Number;
    err : hexa_double;
    use HexaDobl_Complex_Numbers;
   -- use HexaDobl_Complex_Numbers_Polar;
    use HexaDobl_Complex_Matrices;

  begin
    HexaDobl_Complex_Linear_Solvers.Permute_Lower(L,ipvt);
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
        err := AbsVal(dif); -- Radius(dif);
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
