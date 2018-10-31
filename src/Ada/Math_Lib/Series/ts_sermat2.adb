with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;         use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;         use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_QR_Least_Squares;
with Standard_Dense_Series2;
with Standard_Dense_Series2_io;
with Standard_Dense_Series2_Vectors;
with Standard_Dense_Series2_Matrices;
with Standard_Dense_Matrix_Series2;
with Standard_Dense_Matrix_Series2_io;    use Standard_Dense_Matrix_Series2_io;
with Random_Series_Vectors;
with Random_Series_Matrices;
with Standard_Series_Vector_Norms2;
with Standard_Linear_Series2_Solvers;
with Standard_Least_Squares_Series2;

procedure ts_sermat2 is

-- DESCRIPTION :
--   Test on matrices of truncated dense power series.

  procedure Write ( v : in Standard_Dense_Series2_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Basic output of a vector of series in double precision.

  begin
    for i in v'range loop
      put("Component "); put(i,1); put_line(" :");
      Standard_Dense_Series2_io.put(v(i));
    end loop;
  end Write;

  procedure Write ( A : in Standard_Dense_Series2_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("Component ");
        put(i,1); put(", "); put(j,1); put_line(" :");
        Standard_Dense_Series2_io.put(A(i,j));
      end loop;
    end loop;
  end Write;

  procedure Write ( A : in Standard_Complex_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("Component ");
        put(i,1); put(", "); put(j,1); put_line(" :");
        put(A(i,j)); new_line;
      end loop;
    end loop;
  end Write;

  procedure Standard_Random_Linear_Solve ( n,degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random matrix A and right hand side vector b
  --   of dimension n and solves the linear system A*x = b,
  --   in standard double precision.
  --   The series are all of the given degree.

    use Standard_Dense_Series2_Vectors;
    use Standard_Dense_Series2_Matrices;
    use Standard_Linear_Series2_Solvers;

    A : constant Standard_Dense_Series2_Matrices.Matrix(1..n,1..n)
      := Random_Series_Matrices.Random_Series_Matrix(1,n,1,n,degree);
    wrk : Standard_Dense_Series2_Matrices.Matrix(1..n,1..n);
    x : Standard_Dense_Series2_Vectors.Vector(1..n)
      := Random_Series_Vectors.Random_Series_Vector(1,n,degree);
    b : Standard_Dense_Series2_Vectors.Vector(1..n) := A*x;
    y : Standard_Dense_Series2_Vectors.Vector(1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    rsd,err : Standard_Dense_Series2_Vectors.Vector(1..n);
    nrm_rsd,nrm_err : double_float;

  begin
    put_line("The coefficient matrix A :"); Write(A);
    put_line("The right hand side vector b :"); Write(b);
    Standard_Dense_Series2_Matrices.Copy(A,wrk); -- deep copy
    LUfac(wrk,n,ipvt,info);
    put("The pivots : "); put(ipvt); new_line;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Standard_Dense_Series2_Vectors.Copy(b,y); -- deep copy
      LUsolve(wrk,n,ipvt,y);
      put_line("The generated solution x to A*x = b :"); Write(x);
      put_line("The computed solution x to A*x = b :"); Write(y);
      rsd := b - A*y; -- backward error
      put_line("The residual b - A*x :"); Write(rsd);
      err := x - y; -- forward error
      put_line("The difference between generated and computed :"); Write(err);
      nrm_rsd := Standard_Series_Vector_Norms2.Max_Norm(rsd);
      nrm_err := Standard_Series_Vector_Norms2.Max_Norm(err);
      put("Max norm of the backward error :"); put(nrm_rsd,3); new_line;
      put("Max norm of the forward error  :"); put(nrm_err,3); new_line;
    end if;
  end Standard_Random_Linear_Solve;

  function Trunc ( A : Standard_Dense_Series2_Matrices.Matrix )
                 return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Takes the zero-th degree coefficient of every element in A.

    res : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := A(i,j).cff(0);
      end loop;
    end loop;
    return res;
  end Trunc;

  procedure Zero_QRD ( A : Standard_Dense_Series2_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Shows the QR decomposition of the zero degree terms of A.

    B : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2)) := Trunc(A);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    qraux : Standard_Complex_Vectors.Vector(A'range(2));

  begin
    Standard_Complex_QR_Least_Squares.QRD(B,qraux,ipvt,false);
    put_line("QRD on zero degree terms:"); Write(B);
  end Zero_QRD;

  procedure Solve_Normal_Equations
              ( A : in Standard_Dense_Series2_Matrices.Matrix;
                x,b : in Standard_Dense_Series2_Vectors.Vector ) is

  -- DESCRIPTION :
  --   For the overdetermined system A*x = b, formulates the normal
  --   equations A^T A*x = A^T*b and then applies LU factorization
  --   to solve the system of normal equations,
  --   in standard double precision.
  --   The computed solution is compared to the constructed solution x
  --   and the residual of the normal equations is computed as well.

    use Standard_Dense_Series2_Vectors;
    use Standard_Dense_Series2_Matrices;
    use Standard_Linear_Series2_Solvers;

    T : constant Matrix(A'range(2),A'range(1)) := Transpose(A);
    TA : constant Matrix(A'range(2),A'range(2)) := T*A;
    wrk : Matrix(A'range(2),A'range(2));
    Tb : constant Standard_Dense_Series2_Vectors.Vector(A'range(2)) := T*b;
    n : constant integer32 := A'last(2);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    sol,rsd : Standard_Dense_Series2_Vectors.Vector(1..n);

  begin
    Standard_Dense_Series2_Matrices.Copy(TA,wrk);
    LUfac(wrk,n,ipvt,info);
    put("The pivots : "); put(ipvt); new_line;
    if info /= 0 then
      put("info : "); put(info,1); new_line;
    else
      Standard_Dense_Series2_Vectors.Copy(Tb,sol);
      LUsolve(wrk,n,ipvt,sol);
      put_line("The constructed solution x to A*x = b :"); Write(x);
      put_line("The computed solution x to A*x = b :"); Write(sol);
      rsd := Tb - TA*sol;
      put_line("The residual b - A*x :"); Write(rsd);
    end if;
  end Solve_Normal_Equations;

  procedure Test_Normality
              ( qr : in Standard_Dense_Series2_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Given in qr is the orthogonal part of the QR decomposition,
  --   this procedure computes the norm of every column in qr,
  --   in standard double precision.

    use Standard_Dense_Series2;
    use Standard_Dense_Series2_io;
    use Standard_Dense_Series2_Vectors;

    col : Vector(qr'range(1));
    deg : constant integer32 := qr(qr'first(1),qr'first(2)).deg;
    nrm : Series(deg);

  begin
    new_line;
    put_line("Normality test on orthogonal part of QR decomposition ...");
    for k in qr'range(2) loop
      put("Testing column "); put(k,1); put_line(" :");
      for i in col'range loop
        col(i) := qr(i,k);
      end loop;
      nrm := Standard_Series_Vector_Norms2.Norm(col);
      put_line("The norm :"); put(nrm);
    end loop;
  end Test_Normality;

  procedure Test_Orthogonality
              ( qr : in Standard_Dense_Series2_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Given in qr is the orthogonal part of the QR decomposition,
  --   this procedures computes all products of the vectors with
  --   all other vectors to test whether they are orthogonal.

    u,v : Standard_Dense_Series2_Vectors.Vector(qr'range(1));
    deg : constant integer32 := qr(qr'first(1),qr'first(2)).deg;
    ipr : Standard_Dense_Series2.Series(deg);

  begin
    new_line;
    put_line("Orthogonality test on the QR decomposition ...");
    for i in qr'range(2) loop
      for k in u'range loop
        u(k) := qr(k,i);
      end loop;
      for j in i+1..qr'last(2) loop
        for k in v'range loop
          v(k) := qr(k,j);
        end loop;
        ipr := Standard_Series_Vector_Norms2.Inner_Product(u,v);
        put("Inner product of "); put(i,1); put(" with "); put(j,1);
        put_line(" :"); Standard_Dense_Series2_io.put(ipr);
      end loop;
    end loop;
  end Test_Orthogonality;

  procedure Test_Basis
              ( wrk,A : in Standard_Dense_Series2_Matrices.Matrix ) is

  -- DESCRIPION :
  --   Given in wrk the output of QRD on A, the orthogonal part of
  --   the QR decomposition is computed is tested for orthogonality
  --   and normality, in standard double precision.

    use Standard_Dense_Series2_Matrices;

    qr : Matrix(wrk'range(1),wrk'range(2));
    ans : character;

  begin
    Standard_Dense_Series2_Matrices.Copy(wrk,qr);
    new_line;
    put_line("Computing the orthogonal part of the QR decomposition ...");
    Standard_Least_Squares_Series2.Basis(qr,A);
    Test_Normality(qr);   
    new_line;
    put("Continue to the orthogonality test ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Test_Orthogonality(qr);
    end if;
  end Test_Basis;

  procedure QR_Solve_Least_Squares
              ( A : in Standard_Dense_Series2_Matrices.Matrix;
                x,b : in Standard_Dense_Series2_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Applies QR decomposition to the matrix A and then solves the
  --   system A*x = b with this QR decomposition.
  --   The computed solution is compared to the constructed solution
  --   and the residual is computed.

    use Standard_Dense_Series2;
    use Standard_Dense_Series2_Vectors;
    use Standard_Dense_Series2_Matrices;
    use Standard_Least_Squares_Series2;

    deg : constant integer32 := A(A'first(1),A'first(2)).deg;
    wrk : Standard_Dense_Series2_Matrices.Matrix(A'range(1),A'range(2));
    ipvt : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    qraux,sol : Standard_Dense_Series2_Vectors.Vector(A'range(2));
    info : integer32;
    n : constant integer32 := A'last(1);
    m : constant integer32 := A'last(2);
    rsd,dum,dum2,dum3 : Standard_Dense_Series2_Vectors.Vector(1..n);
    err : Standard_Dense_Series2_Vectors.Vector(1..m);
    ans : character;
    nrm_rsd,nrm_err : double_float;

  begin
    Standard_Dense_Series2_Matrices.Copy(A,wrk);
    for i in qraux'range loop
      qraux(i) := new Series'(Standard_Dense_Series2.Create(0,deg));
      sol(i) := new Series'(Standard_Dense_Series2.Create(0,deg));
    end loop;
    QRD(wrk,qraux,ipvt,false);
    put_line("The output of QRD :"); Write(wrk);
    new_line;
    put("View the QRD of the zero degree terms ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Zero_QRD(A);
    end if;
    put_line("The matrix A :"); Write(A);
    new_line;
    put("Test the orthonormality of the basis ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Test_Basis(wrk,A);
    end if;
    new_line;
    put("Continue with the least squares solving ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("the vector b before QRLS :"); Write(b);
      QRLS(wrk,n,m,qraux,b,dum2,dum3,sol,rsd,dum,110,info);
      put_line("the vector b after QRLS :"); Write(b);
      if info /= 0 then
        put("info : "); put(info,1); new_line;
      else
        put_line("The constructed solution x to A*x = b :"); Write(x);
        put_line("The computed solution x to A*x = b :"); Write(sol);
        rsd := b - A*sol; -- backward error
        put_line("The residual b - A*x :"); Write(rsd);
        err := x - sol; -- forward error
        put_line("Difference between constructed and computed :"); Write(err);
        nrm_rsd := Standard_Series_Vector_Norms2.Max_Norm(rsd);
        nrm_err := Standard_Series_Vector_Norms2.Max_Norm(err);
        put("Max norm of backward error : "); put(nrm_rsd,3); new_line;
        put("Max norm of forward error  : "); put(nrm_err,3); new_line;
      end if;
    end if;
  end QR_Solve_Least_Squares;

  procedure Standard_Random_Least_Squares ( n,m,degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-m matrix A with series of the given degree.
  --   A random solution x is generated and then the right hand side
  --   vector b is A*x.  for n > m this is an overdetermined linear
  --   system A*x = b, where the solution x is the generated vector.
  --   Computations are done in standard double precision.

    use Standard_Dense_Series2_Matrices;

    A : constant Standard_Dense_Series2_Matrices.Matrix(1..n,1..m)
      := Random_Series_Matrices.Random_Series_Matrix(1,n,1,m,degree);
    x : constant Standard_Dense_Series2_Vectors.Vector(1..m)
      := Random_Series_Vectors.Random_Series_Vector(1,m,degree);
    b : constant Standard_Dense_Series2_Vectors.Vector(1..n) := A*x;
    ans : character;

  begin
    put_line("The coefficient matrix A :"); Write(A);
    put_line("The solution x :"); Write(x);
    put_line("The right hand side vector b :"); Write(b);
    new_line;
    put("Solve the normal equations ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Solve_Normal_Equations(A,x,b);
    end if;
    new_line;
    put("Solve with a QR decomposition ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then QR_Solve_Least_Squares(A,x,b);
    end if;
  end Standard_Random_Least_Squares;

  procedure Test_Solving is

  -- DESCRIPTION :
  --   Prompts the user for the dimension of the linear system
  --   and the degree of the series, solves a random system.

    dim,deg,nrows,ncols : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the degree of the series : "); get(deg);
    new_line;
    put("Square system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give the dimension : "); get(dim);
      Standard_Random_Linear_Solve(dim,deg);
    else
      put("Give number of rows : "); get(nrows);
      put("Give number of colums : "); get(ncols);
      Standard_Random_Least_Squares(nrows,ncols,deg);
    end if;
  end Test_Solving;

  procedure Standard_Test_Conversion ( n,m,d : in integer32 ) is

  -- DESCRIPTION :
  --   Converts an n-by-m matrix of series of degree d with standard
  --   double precision complex coefficients into a matrix series.

    sA : constant Standard_Dense_Series2_Matrices.Matrix(1..n,1..m)
       := Random_Series_Matrices.Random_Series_Matrix(1,n,1,m,d);
    As : constant Standard_Dense_Matrix_Series2.Matrix 
       := Standard_Dense_Matrix_Series2.Create(sA); 

  begin
    Write(sA);
    put_line("The coefficients of the matrix series :");
    put(As);
  end Standard_Test_Conversion;

  procedure Test_Conversion is

  -- DESCRIPTION :
  --   Prompts for the number of rows, the number of columns,
  --   and the degree of the series.

    nbrows,nbcols,deg : integer32 := 0;

  begin
    new_line;
    put_line("Converting random series matrix into matrix series ...");
    put("  Give the number of rows : "); get(nbrows);
    put("  Give the number of columns : "); get(nbcols);
    put("  Give the degree of the series : "); get(deg);
    Standard_Test_Conversion(nbrows,nbcols,deg);
  end Test_Conversion;

  procedure Main is

  -- DESCRIPTION :
  --   Displays the menu of test procedures and runs the selected test.

    ans : character;

  begin
    new_line;
    put_line("MENU for testing procedures :");
    put_line("  0. test solving of a random linear system;");
    put_line("  1. convert series matrix into matrix series.");
    put("Type 0 or 1 to select a test : "); Ask_Alternative(ans,"01");
    case ans is
      when '0' => Test_Solving;
      when '1' => Test_Conversion;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_sermat2;
