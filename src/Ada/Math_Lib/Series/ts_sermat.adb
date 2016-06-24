with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers_io;         use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;         use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_QR_Least_Squares;
with Standard_Dense_Series;               use Standard_Dense_Series;
with Standard_Dense_Series_io;            use Standard_Dense_Series_io;
with Standard_Dense_Series_Vectors;
with Standard_Dense_Series_Matrices;
with Standard_Random_Series;
with Standard_Series_Vector_Norms;
with Standard_Linear_Series_Solvers;      use Standard_Linear_Series_Solvers;
with Standard_Least_Squares_Series;       use Standard_Least_Squares_Series;

procedure ts_sermat is

-- DESCRIPTION :
--   Test on matrices of truncated dense power series.

  procedure Write ( v : in Standard_Dense_Series_Vectors.Vector ) is
  begin
    for i in v'range loop
      put("Component "); put(i,1); put_line(" :");
      put(v(i));
    end loop;
  end Write;

  procedure Write ( A : in Standard_Dense_Series_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("Component ");
        put(i,1); put(", "); put(j,1); put_line(" :");
        put(A(i,j));
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

  procedure Random_Linear_Solve ( n,order : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random matrix A and right hand side vector b
  --   of dimension n and solves the linear system A*x = b.
  --   The series are all of the given order.

    A : Standard_Dense_Series_Matrices.Matrix(1..n,1..n)
      := Standard_Random_Series.Random_Series_Matrix(1,n,1,n,order);
    wrk : Standard_Dense_Series_Matrices.Matrix(1..n,1..n) := A;
    b : Standard_Dense_Series_Vectors.Vector(1..n)
      := Standard_Random_Series.Random_Series_Vector(1,n,order);
    x : Standard_Dense_Series_Vectors.Vector(1..n) := b;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    rsd : Standard_Dense_Series_Vectors.Vector(1..n);

    use Standard_Dense_Series_Vectors;
    use Standard_Dense_Series_Matrices;

  begin
    put_line("The coefficient matrix A :"); Write(A);
    put_line("The right hand side vector b :"); Write(b);
    LUfac(wrk,n,ipvt,info);
    put("The pivots : "); put(ipvt); new_line;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      LUsolve(wrk,n,ipvt,x);
      put_line("The solution x to A*x = b :"); Write(x);
      rsd := b - A*x;
      put_line("The residual b - A*x :"); Write(rsd);
    end if;
  end Random_Linear_Solve;

  function Trunc ( A : Standard_Dense_Series_Matrices.Matrix )
                 return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Takes the zero-th order coefficient of every element in A.

    res : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := A(i,j).cff(0);
      end loop;
    end loop;
    return res;
  end Trunc;

  procedure Zero_QRD ( A : Standard_Dense_Series_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Shows the QR decomposition of the zero order terms of A.

    B : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2)) := Trunc(A);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    qraux : Standard_Complex_Vectors.Vector(A'range(2));

  begin
    Standard_Complex_QR_Least_Squares.QRD(B,qraux,ipvt,false);
    put_line("QRD on zero order terms:"); Write(B);
  end Zero_QRD;

  procedure Solve_Normal_Equations
              ( A : in Standard_Dense_Series_Matrices.Matrix;
                x,b : in Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   For the overdetermined system A*x = b, formulates the normal
  --   equations A^T A*x = A^T*b and then applies LU factorization
  --   to solve the system of normal equations.
  --   The computed solution is compared to the constructed solution x
  --   and the residual of the normal equations is computed as well.

    use Standard_Dense_Series_Vectors;
    use Standard_Dense_Series_Matrices;

    T : constant Matrix(A'range(2),A'range(1)) := Transpose(A);
    TA : constant Matrix(A'range(2),A'range(2)) := T*A;
    wrk : Matrix(A'range(2),A'range(2)) := TA;
    Tb : constant Standard_Dense_Series_Vectors.Vector(A'range(2)) := T*b;
    n : constant integer32 := A'last(2);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    sol,rsd : Standard_Dense_Series_Vectors.Vector(1..n);

  begin
    LUfac(wrk,n,ipvt,info);
    put("The pivots : "); put(ipvt); new_line;
    if info /= 0 then
      put("info : "); put(info,1); new_line;
    else
      sol := Tb;
      LUsolve(wrk,n,ipvt,sol);
      put_line("The constructed solution x to A*x = b :"); Write(x);
      put_line("The computed solution x to A*x = b :"); Write(sol);
      rsd := Tb - TA*sol;
      put_line("The residual b - A*x :"); Write(rsd);
    end if;
  end Solve_Normal_Equations;

  procedure Test_Normality
              ( qr : in Standard_Dense_Series_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Given in qr is the orthogonal part of the QR decomposition,
  --   this procedures computes the norm of every column in qr.

    use Standard_Dense_Series_Vectors;

    col : Vector(qr'range(1));
    nrm : Series;

  begin
    new_line;
    put_line("Normality test on orthogonal part of QR decomposition ...");
    for k in qr'range(2) loop
      put("Testing column "); put(k,1); put_line(" :");
      for i in col'range loop
        col(i) := qr(i,k);
      end loop;
      nrm := Standard_Series_Vector_Norms.Norm(col);
      put_line("The norm :"); put(nrm);
    end loop;
  end Test_Normality;

  procedure Test_Basis
              ( wrk,A : in Standard_Dense_Series_Matrices.Matrix ) is

  -- DESCRIPION :
  --   Given in wrk the output of QRD on A, the orthogonal part of
  --   the QR decomposition is computed is tested for orthogonality
  --   and normality.

    use Standard_Dense_Series_Matrices;

    qr : Matrix(wrk'range(1),wrk'range(2)) := wrk;

  begin
    new_line;
    put_line("Computing the orthogonal part of the QR decomposition ...");
    Basis(qr,A);
    Test_Normality(qr);   
  end Test_Basis;

  procedure QR_Solve_Least_Squares
              ( A : in Standard_Dense_Series_Matrices.Matrix;
                x,b : in Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Applies QR decomposition to the matrix A and then solves the
  --   system A*x = b with this QR decomposition.
  --   The computed solution is compared to the constructed solution
  --   and the residual is computed.

    use Standard_Dense_Series_Vectors;
    use Standard_Dense_Series_Matrices;

    wrk : Standard_Dense_Series_Matrices.Matrix(A'range(1),A'range(2)) := A;
    ipvt : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    qraux,sol : Standard_Dense_Series_Vectors.Vector(A'range(2));
    info : integer32;
    n : constant integer32 := A'last(1);
    m : constant integer32 := A'last(2);
    rsd,dum,dum2,dum3 : Standard_Dense_Series_Vectors.Vector(1..n);
    ans : character;

  begin
    QRD(wrk,qraux,ipvt,false);
    put_line("The output of QRD :"); Write(wrk);
    new_line;
    put("View the QRD of the zero order terms ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Zero_QRD(A);
    end if;
   -- put_line("The matrix A :"); Write(A);
    new_line;
    put("Test the orthonormality of the basis ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Test_Basis(wrk,A);
    end if;
    QRLS(wrk,n,m,qraux,b,dum2,dum3,sol,rsd,dum,110,info);
    if info /= 0 then
      put("info : "); put(info,1); new_line;
    else
      put_line("The constructed solution x to A*x = b :"); Write(x);
      put_line("The computed solution x to A*x = b :"); Write(sol);
      rsd := b - A*sol;
      put_line("The residual b - A*x :"); Write(rsd);
    end if;
  end QR_Solve_Least_Squares;

  procedure Random_Least_Squares ( n,m,order : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-m matrix A with series of the given order.
  --   A random solution x is generated and then the right hand side
  --   vector b is A*x.  for n > m this is an overdetermined linear
  --   system A*x = b, where the solution x is the generated vector.

    use Standard_Dense_Series_Matrices;

    A : constant Standard_Dense_Series_Matrices.Matrix(1..n,1..m)
      := Standard_Random_Series.Random_Series_Matrix(1,n,1,m,order);
    x : constant Standard_Dense_Series_Vectors.Vector(1..m)
      := Standard_Random_Series.Random_Series_Vector(1,m,order);
    b : constant Standard_Dense_Series_Vectors.Vector(1..n) := A*x;
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
  end Random_Least_Squares;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension of the linear system
  --   and the order of the series.

    dim,ord,nrows,ncols : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the order of the series : "); get(ord);
    new_line;
    put("Square system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give the dimension : "); get(dim);
      Random_Linear_Solve(dim,ord);
    else
      put("Give number of rows : "); get(nrows);
      put("Give number of colums : "); get(ncols);
      Random_Least_Squares(nrows,ncols,ord);
    end if;
  end Main;

begin
  Main;
end ts_sermat;
