with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Floating_Two_Norms;        use Standard_Floating_Two_Norms;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Complex_Matrices;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Standard_Floating_Eigenvalues;      use Standard_Floating_Eigenvalues;

procedure ts_eigval is

-- DESCRIPTION :
--   Test on the computation of all eigenvalues and eigenvectors.

  procedure Read_Matrix ( A : in out Standard_Floating_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Allows user to interactively enter all data to a matrix.

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("Give value for ("); put(i,1); put(",");
        put(j,1); put(") : "); get(A(i,j));
      end loop;
    end loop;
  end Read_Matrix;

  function Create ( A : Standard_Floating_Matrices.Matrix )
                  return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Creates a complex matrix with same content as A.

    res : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := Create(A(i,j));
      end loop;
    end loop;
    return res;
  end Create;

  procedure Test_Real_Eigenvalue
               ( A : in Standard_Floating_Matrices.Matrix;
                 lambda : in double_float;
                 x : in Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Checks whether A*x = lambda*x holds for a real eigenvalue
  --   lambda and corresponding real eigenvector x.

    use Standard_Floating_Vectors;
    use Standard_Floating_Matrices;

    y1 : constant Vector(A'range(1)) := A*x;
    y2 : constant Vector(A'range(1)) := lambda*x;
    nrm2 : constant double_float := Norm2(y1-y2);

  begin
   -- put_line("A*x : "); put_line(y1);
   -- put_line("lambda*x : "); put_line(y2);
    put("2-norm of A*x - lambda*x : "); put(nrm2,3);
    if nrm2 < 1.0E-12
     then put_line("  Okay");
     else put_line("  BUG!");
    end if;
  end Test_Real_Eigenvalue;

  procedure Test_Complex_Eigenvalue
               ( A : in Standard_Complex_Matrices.Matrix;
                 lambda : in Complex_Number;
                 x : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Checks whether A*x = lambda*x holds for a complex eigenvalue
  --   lambda and corresponding complex eigenvector x.

    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;

    y1 : constant Vector(A'range(1)) := A*x;
    y2 : constant Vector(A'range(1)) := lambda*x;
    nrm2 : constant double_float := Norm2(y1-y2);

  begin
   -- put_line("A*x : "); put_line(y1);
   -- put_line("lambda*x : "); put_line(y2);
    put("2-norm of A*x - lambda*x : "); put(nrm2,3);
    if nrm2 < 1.0E-12
     then put_line("  Okay");
     else put_line("  BUG!");
    end if;
  end Test_Complex_Eigenvalue;

  procedure Eigen_Equation
             ( A : in Standard_Complex_Matrices.Matrix;
               L : in Standard_Complex_Vectors.Vector;
               v : in Standard_Complex_VecVecs.Vecvec ) is

  -- DESCRIPTION :
  --   Tests the eigenvalue-eigenvector equation A*x - lambda*x = 0,
  --   for all eigenvectors and eigenvalues.

  begin
    for i in L'range loop
      put("lambda("); put(i,1); put(") : "); put(L(i)); new_line;
      Test_Complex_Eigenvalue(A,L(i),v(i).all);
    end loop;
  end Eigen_Equation;

  procedure Eigen_Equation
             ( A,z : in Standard_Floating_Matrices.Matrix; 
               wr,wi : in Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Tests the eigenvalue-eigenvector equation A*x - lambda*x = 0,
  --   for all eigenvectors and eigenvalues.

    x : Standard_Floating_Vectors.Vector(A'range(1));
    v : constant Standard_Complex_VecVecs.VecVec(z'range(1)) := Create(z,wi);
    B : constant Standard_Complex_Matrices.Matrix(A'range(1),A'range(2))
      := Create(A);

  begin
    for i in wr'range loop
      put("lambda("); put(i,1); put(") : ");
      put(wr(i)); put("  "); put(wi(i));
      if wi(i) + 1.0 = 1.0 then
        put_line(" is real");
        for j in x'range loop
          x(j) := z(j,i);
        end loop;
        Test_Real_Eigenvalue(A,wr(i),x);
      else
        put_line(" is not real");
        Test_Complex_Eigenvalue(B,Create(wr(i),wi(i)),v(i).all);
      end if;
    end loop;
  end Eigen_Equation;

  procedure Compute_Eigenvalues
               ( n : in integer32; 
                 A : in Standard_Floating_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Eigenvalues will be computed twice, the second time with the
  --   option of RG that also returns the eigenvectors.

    z : Standard_Floating_Matrices.Matrix(1..n,1..n);
    wr1,wi1,wr2,wi2 : Standard_Floating_Vectors.Vector(1..n);
    realdiff,imagdiff : double_float := 0.0;
    ierr : integer32;

  begin
    RG(n,n,A,0,wr1,wi1,z,ierr);
    put("ierr : "); put(ierr,1); new_line;
    put_line("The eigenvalues : ");
    for i in 1..n loop
      put(wr1(i)); put("  "); put(wi1(i)); new_line;
    end loop;
    put_line("for testing purposes we recompute for eigenvectors...");
    RG(n,n,A,1,wr2,wi2,z,ierr);
    put("ierr : "); put(ierr,1); new_line;
    put_line("The eigenvalues : ");
    for i in 1..n loop
      put(wr2(i)); put("  "); put(wi2(i)); new_line;
      realdiff := realdiff + AbsVal(wr2(i) - wr1(i));
      imagdiff := imagdiff + AbsVal(wi2(i) - wi1(i));
    end loop;
    put("Differences real :"); put(realdiff,3);
    put(" and imaginary :"); put(imagdiff,3);
    if realdiff + imagdiff < 1.0E-12
     then put_line("  Okay");
     else put_line("  BUG!");
    end if;
    put_line("The eigenvectors : "); put(z);
    Eigen_Equation(A,z,wr2,wi2);
  end Compute_Eigenvalues;

  procedure Test_Wrappers ( A : in Standard_Floating_Matrices.Matrix ) is

    L : Standard_Complex_Vectors.Vector(A'range(1));
    v : Standard_Complex_VecVecs.VecVec(A'range(1));
    B : constant Standard_Complex_Matrices.Matrix(A'range(1),A'range(2))
      := Create(A);
    err : integer32;

  begin
    Eigenvalues(A,err,L);
    put_line("The eigenvalues : "); put_line(L);
    Eigenvectors(A,err,L,v);
    put_line("Recomputed eigenvalues : "); put_line(L);
    Eigen_Equation(B,L,v);
  end Test_Wrappers;

  procedure Main is

    ans : character;
    n : integer32 := 0;

  begin
    new_line;
    put_line("MENU to test eigenvalue/eigenvector computations :");
    put_line("  1. user types in the value for a matrix; or");
    put_line("  2. compute eigenvalues for a random matrix; or");
    put_line("  3. test wrappers for random matrices.");
    put("Type 1, 2, or 3 to choose : "); Ask_Alternative(ans,"123");
    new_line;
    put("Give the dimension : "); get(n);
    declare
      A : Standard_Floating_Matrices.Matrix(1..n,1..n);
    begin
      if ans = '1' then
        Read_Matrix(A);
      else
        A := Random_Matrix(natural32(n),natural32(n));
        if ans = '2' then
          put_line("The matrix : "); put(A);
          Compute_Eigenvalues(n,A);
        else
          Test_Wrappers(A);
        end if;
      end if;
    end;
  end Main;

begin
  Main;
end ts_eigval;
