with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Floating_Vectors;
with Double_Double_Vectors;
with Quad_Double_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Floating_Two_Norms;        use Standard_Floating_Two_Norms;
with Double_Double_Two_Norms;            use Double_Double_Two_Norms;
with Quad_Double_Two_Norms;              use Quad_Double_Two_Norms;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with DoblDobl_Complex_Vector_Norms;      use DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Double_Double_Matrices;
with Double_Double_Matrices_io;          use Double_Double_Matrices_io;
with Quad_Double_Matrices;
with Quad_Double_Matrices_io;            use Quad_Double_Matrices_io;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with DoblDobl_Random_Matrices;           use DoblDobl_Random_Matrices;
with QuadDobl_Random_Matrices;           use QuadDobl_Random_Matrices;
with Standard_Floating_Eigenvalues;
with Double_Double_Eigenvalues;
with Quad_Double_Eigenvalues;

procedure ts_eigval is

-- DESCRIPTION :
--   Test on the computation of all eigenvalues and eigenvectors.

  procedure Read_Matrix ( A : in out Standard_Floating_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Allows user to interactively enter all data to a matrix,
  --   of real numbers in standard double precision.

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("Give value for ("); put(i,1); put(",");
        put(j,1); put(") : "); get(A(i,j));
      end loop;
    end loop;
  end Read_Matrix;

  procedure Read_Matrix ( A : in out Double_Double_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Allows user to interactively enter all data to a matrix
  --   of real numbers in double double precision.

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("Give value for ("); put(i,1); put(",");
        put(j,1); put(") : "); get(A(i,j));
      end loop;
    end loop;
  end Read_Matrix;

  procedure Read_Matrix ( A : in out Quad_Double_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Allows user to interactively enter all data to a matrix
  --   of real numbers in quad double precision.

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
  --   Creates a complex matrix with the same content as A,
  --   in standard double precision.

    res : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := Standard_Complex_Numbers.Create(A(i,j));
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( A : Double_Double_Matrices.Matrix )
                  return DoblDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Creates a complex matrix with the same content as A,
  --   in double double precision.

    res : DoblDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := DoblDobl_Complex_Numbers.Create(A(i,j));
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( A : Quad_Double_Matrices.Matrix )
                  return QuadDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Creates a complex matrix with the same content as A,
  --   in quad double precision.

    res : QuadDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := QuadDobl_Complex_Numbers.Create(A(i,j));
      end loop;
    end loop;
    return res;
  end Create;

  procedure Test_Real_Eigenvalue
               ( A : in Standard_Floating_Matrices.Matrix;
                 lambda : in Standard_Floating_Numbers.double_float;
                 x : in Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Checks whether A*x = lambda*x holds for a real eigenvalue
  --   lambda and corresponding real eigenvector x,
  --   in standard double precision.

    use Standard_Floating_Numbers;
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

  procedure Test_Real_Eigenvalue
               ( A : in Double_Double_Matrices.Matrix;
                 lambda : in Double_Double_Numbers.double_double;
                 x : in Double_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Checks whether A*x = lambda*x holds for a real eigenvalue
  --   lambda and corresponding real eigenvector x,
  --   in double double precision.

    use Double_Double_Numbers;
    use Double_Double_Vectors;
    use Double_Double_Matrices;

    y1 : constant Vector(A'range(1)) := A*x;
    y2 : constant Vector(A'range(1)) := lambda*x;
    nrm2 : constant double_double := Norm2(y1-y2);

  begin
   -- put_line("A*x : "); put_line(y1);
   -- put_line("lambda*x : "); put_line(y2);
    put("2-norm of A*x - lambda*x : "); put(nrm2,3);
    if nrm2 < 1.0E-12
     then put_line("  Okay");
     else put_line("  BUG!");
    end if;
  end Test_Real_Eigenvalue;

  procedure Test_Real_Eigenvalue
               ( A : in Quad_Double_Matrices.Matrix;
                 lambda : in Quad_Double_Numbers.quad_double;
                 x : in Quad_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Checks whether A*x = lambda*x holds for a real eigenvalue
  --   lambda and corresponding real eigenvector x,
  --   in quad double precision.

    use Quad_Double_Numbers;
    use Quad_Double_Vectors;
    use Quad_Double_Matrices;

    y1 : constant Vector(A'range(1)) := A*x;
    y2 : constant Vector(A'range(1)) := lambda*x;
    nrm2 : constant quad_double := Norm2(y1-y2);

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
                 lambda : in Standard_Complex_Numbers.Complex_Number;
                 x : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Checks whether A*x = lambda*x holds for a complex eigenvalue
  --   lambda and corresponding complex eigenvector x,
  --   in standard double precision.

    use Standard_Floating_Numbers;
    use Standard_Complex_Numbers;
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

  procedure Test_Complex_Eigenvalue
               ( A : in DoblDobl_Complex_Matrices.Matrix;
                 lambda : in DoblDobl_Complex_Numbers.Complex_Number;
                 x : in DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Checks whether A*x = lambda*x holds for a complex eigenvalue
  --   lambda and corresponding complex eigenvector x,
  --   in double double precision.

    use Double_Double_Numbers;
    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;

    y1 : constant Vector(A'range(1)) := A*x;
    y2 : constant Vector(A'range(1)) := lambda*x;
    nrm2 : constant double_double := Norm2(y1-y2);

  begin
   -- put_line("A*x : "); put_line(y1);
   -- put_line("lambda*x : "); put_line(y2);
    put("2-norm of A*x - lambda*x : "); put(nrm2,3);
    if nrm2 < 1.0E-12
     then put_line("  Okay");
     else put_line("  BUG!");
    end if;
  end Test_Complex_Eigenvalue;

  procedure Test_Complex_Eigenvalue
               ( A : in QuadDobl_Complex_Matrices.Matrix;
                 lambda : in QuadDobl_Complex_Numbers.Complex_Number;
                 x : in QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Checks whether A*x = lambda*x holds for a complex eigenvalue
  --   lambda and corresponding complex eigenvector x,
  --   in double double precision.

    use Quad_Double_Numbers;
    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Matrices;

    y1 : constant Vector(A'range(1)) := A*x;
    y2 : constant Vector(A'range(1)) := lambda*x;
    nrm2 : constant quad_double := Norm2(y1-y2);

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
  --   for all eigenvectors and eigenvalues, in standard double precision.

  begin
    for i in L'range loop
      put("lambda("); put(i,1); put(") : "); put(L(i)); new_line;
      Test_Complex_Eigenvalue(A,L(i),v(i).all);
    end loop;
  end Eigen_Equation;

  procedure Eigen_Equation
             ( A : in DoblDobl_Complex_Matrices.Matrix;
               L : in DoblDobl_Complex_Vectors.Vector;
               v : in DoblDobl_Complex_VecVecs.Vecvec ) is

  -- DESCRIPTION :
  --   Tests the eigenvalue-eigenvector equation A*x - lambda*x = 0,
  --   for all eigenvectors and eigenvalues, in double double precision.

  begin
    for i in L'range loop
      put("lambda("); put(i,1); put(") : "); put(L(i)); new_line;
      Test_Complex_Eigenvalue(A,L(i),v(i).all);
    end loop;
  end Eigen_Equation;

  procedure Eigen_Equation
             ( A : in QuadDobl_Complex_Matrices.Matrix;
               L : in QuadDobl_Complex_Vectors.Vector;
               v : in QuadDobl_Complex_VecVecs.Vecvec ) is

  -- DESCRIPTION :
  --   Tests the eigenvalue-eigenvector equation A*x - lambda*x = 0,
  --   for all eigenvectors and eigenvalues, in quad double precision.

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
  --   for all eigenvectors and eigenvalues, in standard double precision.

    use Standard_Floating_Numbers;
    use Standard_Complex_Numbers;

    x : Standard_Floating_Vectors.Vector(A'range(1));
    v : constant Standard_Complex_VecVecs.VecVec(z'range(1))
      := Standard_Floating_Eigenvalues.Create(z,wi);
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

  procedure Eigen_Equation
             ( A,z : in Double_Double_Matrices.Matrix; 
               wr,wi : in Double_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Tests the eigenvalue-eigenvector equation A*x - lambda*x = 0,
  --   for all eigenvectors and eigenvalues, in double double precision.

    use Double_Double_Numbers;
    use DoblDobl_Complex_Numbers;

    x : Double_Double_Vectors.Vector(A'range(1));
    v : constant DoblDobl_Complex_VecVecs.VecVec(z'range(1))
      := Double_Double_Eigenvalues.Create(z,wi);
    B : constant DoblDobl_Complex_Matrices.Matrix(A'range(1),A'range(2))
      := Create(A);
    one : constant double_double := create(1.0);

  begin
    for i in wr'range loop
      put("lambda("); put(i,1); put(") : ");
      put(wr(i)); put("  "); put(wi(i));
      if wi(i) + one = one then
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


  procedure Eigen_Equation
             ( A,z : in Quad_Double_Matrices.Matrix; 
               wr,wi : in Quad_Double_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Tests the eigenvalue-eigenvector equation A*x - lambda*x = 0,
  --   for all eigenvectors and eigenvalues, in double double precision.

    use Quad_Double_Numbers;
    use QuadDobl_Complex_Numbers;

    x : Quad_Double_Vectors.Vector(A'range(1));
    v : constant QuadDobl_Complex_VecVecs.VecVec(z'range(1))
      := Quad_Double_Eigenvalues.Create(z,wi);
    B : constant QuadDobl_Complex_Matrices.Matrix(A'range(1),A'range(2))
      := Create(A);
    one : constant quad_double := create(1.0);

  begin
    for i in wr'range loop
      put("lambda("); put(i,1); put(") : ");
      put(wr(i)); put("  "); put(wi(i));
      if wi(i) + one = one then
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
  --   Eigenvalues will be computed twice, in standard double precision,
  --   the second time with the option that returns the eigenvectors.

    use Standard_Floating_Numbers;
    use Standard_Floating_Eigenvalues;

    z : Standard_Floating_Matrices.Matrix(1..n,1..n);
    wr1,wi1,wr2,wi2 : Standard_Floating_Vectors.Vector(1..n);
    realdiff,imagdiff : double_float := 0.0;
    tol : constant double_float := 1.0E-12;
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
    if realdiff + imagdiff < tol
     then put_line("  Okay");
     else put_line("  BUG!");
    end if;
    put_line("The eigenvectors : "); put(z);
    Eigen_Equation(A,z,wr2,wi2);
  end Compute_Eigenvalues;

  procedure Compute_Eigenvalues
               ( n : in integer32; 
                 A : in Double_Double_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Eigenvalues will be computed twice, in double double precision,
  --   the second time with the option that returns the eigenvectors.

    use Double_Double_Numbers;
    use Double_Double_Eigenvalues;

    z : Double_Double_Matrices.Matrix(1..n,1..n);
    wr1,wi1,wr2,wi2 : Double_Double_Vectors.Vector(1..n);
    realdiff,imagdiff : double_double := create(0.0);
    tol : constant double_double := create(1.0E-12);
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
    put("Differences real : "); put(realdiff,3);
    put(" and imaginary : "); put(imagdiff,3);
    if realdiff + imagdiff < tol
     then put_line("  Okay");
     else put_line("  BUG!");
    end if;
    put_line("The eigenvectors : "); put(z);
    Eigen_Equation(A,z,wr2,wi2);
  end Compute_Eigenvalues;

  procedure Compute_Eigenvalues
               ( n : in integer32; 
                 A : in Quad_Double_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Eigenvalues will be computed twice, in quad double precision,
  --   the second time with the option that returns the eigenvectors.

    use Quad_Double_Numbers;
    use Quad_Double_Eigenvalues;

    z : Quad_Double_Matrices.Matrix(1..n,1..n);
    wr1,wi1,wr2,wi2 : Quad_Double_Vectors.Vector(1..n);
    realdiff,imagdiff : quad_double := create(0.0);
    tol : constant quad_double := create(1.0E-12);
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
    put("Differences real : "); put(realdiff,3);
    put(" and imaginary : "); put(imagdiff,3);
    if realdiff + imagdiff < tol
     then put_line("  Okay");
     else put_line("  BUG!");
    end if;
    put_line("The eigenvectors : "); put(z);
    Eigen_Equation(A,z,wr2,wi2);
  end Compute_Eigenvalues;

  procedure Test_Wrappers ( A : in Standard_Floating_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   The wrappers to the EISPACK routines provide longer names and
  --   write the results in the proper data types.
  --   This test happens in standard double precision.

    L : Standard_Complex_Vectors.Vector(A'range(1));
    v : Standard_Complex_VecVecs.VecVec(A'range(1));
    B : constant Standard_Complex_Matrices.Matrix(A'range(1),A'range(2))
      := Create(A);
    err : integer32;

  begin
    Standard_Floating_Eigenvalues.Eigenvalues(A,err,L);
    put_line("The eigenvalues : "); put_line(L);
    Standard_Floating_Eigenvalues.Eigenvectors(A,err,L,v);
    put_line("Recomputed eigenvalues : "); put_line(L);
    Eigen_Equation(B,L,v);
  end Test_Wrappers;

  procedure Test_Wrappers ( A : in Double_Double_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   The wrappers to the EISPACK routines provide longer names and
  --   write the results in the proper data types.
  --   This test happens in double double precision.

    L : DoblDobl_Complex_Vectors.Vector(A'range(1));
    v : DoblDobl_Complex_VecVecs.VecVec(A'range(1));
    B : constant DoblDobl_Complex_Matrices.Matrix(A'range(1),A'range(2))
      := Create(A);
    err : integer32;

  begin
    Double_Double_Eigenvalues.Eigenvalues(A,err,L);
    put_line("The eigenvalues : "); put_line(L);
    Double_Double_Eigenvalues.Eigenvectors(A,err,L,v);
    put_line("Recomputed eigenvalues : "); put_line(L);
    Eigen_Equation(B,L,v);
  end Test_Wrappers;

  procedure Test_Wrappers ( A : in Quad_Double_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   The wrappers to the EISPACK routines provide longer names and
  --   write the results in the proper data types.
  --   This test happens in quad double precision.

    L : QuadDobl_Complex_Vectors.Vector(A'range(1));
    v : QuadDobl_Complex_VecVecs.VecVec(A'range(1));
    B : constant QuadDobl_Complex_Matrices.Matrix(A'range(1),A'range(2))
      := Create(A);
    err : integer32;

  begin
    Quad_Double_Eigenvalues.Eigenvalues(A,err,L);
    put_line("The eigenvalues : "); put_line(L);
    Quad_Double_Eigenvalues.Eigenvectors(A,err,L,v);
    put_line("Recomputed eigenvalues : "); put_line(L);
    Eigen_Equation(B,L,v);
  end Test_Wrappers;

  function Prompt_for_Test_Type return character is

  -- DESCRIPTION :
  --   Prompts the user for the menu of available tests and returns
  --   a character, which is either '1', '2', or '3', corresponding
  --   to the selected option in the test menu.

    res : character;

  begin
    new_line;
    put_line("MENU to test eigenvalue/eigenvector computations :");
    put_line("  1. user types in the value for a matrix; or");
    put_line("  2. compute eigenvalues for a random matrix; or");
    put_line("  3. test wrappers for random matrices.");
    put("Type 1, 2, or 3 to choose : "); Ask_Alternative(res,"123");
    return res;
  end Prompt_for_Test_Type;

  procedure Standard_Test ( n : in integer32; test_type : in character ) is

  -- DESCRIPTION :
  --   Tests the eigenvalue computation in standard double precision
  --   on an n-dimensional square matrix.
  --   The type of the test is defined by the character test_type.

    A : Standard_Floating_Matrices.Matrix(1..n,1..n);

  begin
    if test_type = '1' then
      Read_Matrix(A);
    else
      A := Random_Matrix(natural32(n),natural32(n));
      if test_type = '2' then
        put_line("The matrix : "); put(A);
        Compute_Eigenvalues(n,A);
      else
        Test_Wrappers(A);
      end if;
    end if;
  end Standard_Test;

  procedure DoblDobl_Test ( n : in integer32; test_type : in character ) is

  -- DESCRIPTION :
  --   Tests the eigenvalue computation in double double precision
  --   on an n-dimensional square matrix.
  --   The type of the test is defined by the character test_type.

    A : Double_Double_Matrices.Matrix(1..n,1..n);

  begin
    if test_type = '1' then
      Read_Matrix(A);
    else
      A := Random_Matrix(natural32(n),natural32(n));
      if test_type = '2' then
        put_line("The matrix : "); put(A);
        Compute_Eigenvalues(n,A);
      else
        Test_Wrappers(A);
      end if;
    end if;
  end DoblDobl_Test;

  procedure QuadDobl_Test ( n : in integer32; test_type : in character ) is

  -- DESCRIPTION :
  --   Tests the eigenvalue computation in quad double precision
  --   on an n-dimensional square matrix.
  --   The type of the test is defined by the character test_type.

    A : Quad_Double_Matrices.Matrix(1..n,1..n);

  begin
    if test_type = '1' then
      Read_Matrix(A);
    else
      A := Random_Matrix(natural32(n),natural32(n));
      if test_type = '2' then
        put_line("The matrix : "); put(A);
        Compute_Eigenvalues(n,A);
      else
        Test_Wrappers(A);
      end if;
    end if;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the type of test, the dimension, and the precision.
  --   Then launches the proper test procedure.

    tstp : constant character := Prompt_for_Test_Type;
    prec : character;
    n : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(n);
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(prec,"012");
    case prec is
      when '0' => Standard_Test(n,tstp);
      when '1' => DoblDobl_Test(n,tstp);
      when '2' => QuadDobl_Test(n,tstp);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_eigval;
