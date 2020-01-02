with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;         use Standard_Complex_Numbers_io;
with Double_Double_Numbers;               use Double_Double_Numbers;
with Double_Double_Numbers_io;            use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers_io;         use DoblDobl_Complex_Numbers_io;
with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with Quad_Double_Numbers_io;              use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers_io;         use QuadDobl_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;         use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_QR_Least_Squares;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_QR_Least_Squares;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_QR_Least_Squares;
with Standard_Dense_Series;
with Standard_Dense_Series_io;
with Standard_Dense_Series_Vectors;
with Standard_Dense_Series_Matrices;
with DoblDobl_Dense_Series;
with DoblDobl_Dense_Series_io;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Dense_Series_Matrices;
with QuadDobl_Dense_Series;
with QuadDobl_Dense_Series_io;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series_Matrices;
with Standard_Random_Series;
with DoblDobl_Random_Series;
with QuadDobl_Random_Series;
with Standard_Series_Vector_Norms;
with Standard_Linear_Series_Solvers;
with DoblDobl_Series_Vector_Norms;
with DoblDobl_Linear_Series_Solvers;
with QuadDobl_Series_Vector_Norms;
with QuadDobl_Linear_Series_Solvers;
with Standard_Least_Squares_Series;
with DoblDobl_Least_Squares_Series;
with QuadDobl_Least_Squares_Series;
with Standard_Dense_Matrix_Series;
with Standard_Dense_Matrix_Series_io;     use Standard_Dense_Matrix_Series_io;
with DoblDobl_Dense_Matrix_Series;
with DoblDobl_Dense_Matrix_Series_io;     use DoblDobl_Dense_Matrix_Series_io;
with QuadDobl_Dense_Matrix_Series;
with QuadDobl_Dense_Matrix_Series_io;     use QuadDobl_Dense_Matrix_Series_io;

procedure ts_sermat is

-- DESCRIPTION :
--   Test on matrices of truncated dense power series.

  procedure Write ( v : in Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Basic output of a vector of series in double precision.

  begin
    for i in v'range loop
      put("Component "); put(i,1); put_line(" :");
      Standard_Dense_Series_io.put(v(i));
    end loop;
  end Write;

  procedure Write ( v : in DoblDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Basic output of a vector of series in double double precision.

  begin
    for i in v'range loop
      put("Component "); put(i,1); put_line(" :");
      DoblDobl_Dense_Series_io.put(v(i));
    end loop;
  end Write;

  procedure Write ( v : in QuadDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Basic output of a vector of series in quad double precision.

  begin
    for i in v'range loop
      put("Component "); put(i,1); put_line(" :");
      QuadDobl_Dense_Series_io.put(v(i));
    end loop;
  end Write;

  procedure Write ( A : in Standard_Dense_Series_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("Component ");
        put(i,1); put(", "); put(j,1); put_line(" :");
        Standard_Dense_Series_io.put(A(i,j));
      end loop;
    end loop;
  end Write;

  procedure Write ( A : in DoblDobl_Dense_Series_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("Component ");
        put(i,1); put(", "); put(j,1); put_line(" :");
        DoblDobl_Dense_Series_io.put(A(i,j));
      end loop;
    end loop;
  end Write;

  procedure Write ( A : in QuadDobl_Dense_Series_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("Component ");
        put(i,1); put(", "); put(j,1); put_line(" :");
        QuadDobl_Dense_Series_io.put(A(i,j));
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

  procedure Write ( A : in DoblDobl_Complex_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("Component ");
        put(i,1); put(", "); put(j,1); put_line(" :");
        put(A(i,j)); new_line;
      end loop;
    end loop;
  end Write;

  procedure Write ( A : in QuadDobl_Complex_Matrices.Matrix ) is
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

    use Standard_Dense_Series_Vectors;
    use Standard_Dense_Series_Matrices;
    use Standard_Linear_Series_Solvers;

    A : constant Standard_Dense_Series_Matrices.Matrix(1..n,1..n)
      := Standard_Random_Series.Random_Series_Matrix(1,n,1,n,degree);
    wrk : Standard_Dense_Series_Matrices.Matrix(1..n,1..n) := A;
    x : constant Standard_Dense_Series_Vectors.Vector(1..n)
      := Standard_Random_Series.Random_Series_Vector(1,n,degree);
    b : constant Standard_Dense_Series_Vectors.Vector(1..n) := A*x;
    y : Standard_Dense_Series_Vectors.Vector(1..n) := b;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    rsd,err : Standard_Dense_Series_Vectors.Vector(1..n);
    nrm_rsd,nrm_err : double_float;

  begin
    put_line("The coefficient matrix A :"); Write(A);
    put_line("The right hand side vector b :"); Write(b);
    LUfac(wrk,n,ipvt,info);
    put("The pivots : "); put(ipvt); new_line;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      LUsolve(wrk,n,ipvt,y);
      put_line("The generated solution x to A*x = b :"); Write(x);
      put_line("The computed solution x to A*x = b :"); Write(y);
      rsd := b - A*y; -- backward error
      put_line("The residual b - A*x :"); Write(rsd);
      err := x - y; -- forward error
      put_line("The difference between generated and computed :"); Write(err);
      nrm_rsd := Standard_Series_Vector_Norms.Max_Norm(rsd);
      nrm_err := Standard_Series_Vector_Norms.Max_Norm(err);
      put("Max norm of the backward error :"); put(nrm_rsd,3); new_line;
      put("Max norm of the forward error  :"); put(nrm_err,3); new_line;
    end if;
  end Standard_Random_Linear_Solve;

  procedure DoblDobl_Random_Linear_Solve ( n,degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random matrix A and right hand side vector b
  --   of dimension n and solves the linear system A*x = b,
  --   in double double precision.
  --   The series are all of the given degree.

    use DoblDobl_Dense_Series_Vectors;
    use DoblDobl_Dense_Series_Matrices;
    use DoblDobl_Linear_Series_Solvers;

    A : constant DoblDobl_Dense_Series_Matrices.Matrix(1..n,1..n)
      := DoblDobl_Random_Series.Random_Series_Matrix(1,n,1,n,degree);
    wrk : DoblDobl_Dense_Series_Matrices.Matrix(1..n,1..n) := A;
    x : constant DoblDobl_Dense_Series_Vectors.Vector(1..n)
      := DoblDobl_Random_Series.Random_Series_Vector(1,n,degree);
    b : constant DoblDobl_Dense_Series_Vectors.Vector(1..n) := A*x;
    y : DoblDobl_Dense_Series_Vectors.Vector(1..n) := b;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    rsd,err : DoblDobl_Dense_Series_Vectors.Vector(1..n);
    nrm_rsd,nrm_err : double_double;

  begin
    put_line("The coefficient matrix A :"); Write(A);
    put_line("The right hand side vector b :"); Write(b);
    LUfac(wrk,n,ipvt,info);
    put("The pivots : "); put(ipvt); new_line;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      LUsolve(wrk,n,ipvt,y);
      put_line("The generated solution x to A*x = b :"); Write(x);
      put_line("The computed solution x to A*x = b :"); Write(y);
      rsd := b - A*y; -- backward error
      put_line("The residual b - A*x :"); Write(rsd);
      err := x - y; -- forward error
      put_line("Difference between generated and computed :"); Write(err);
      nrm_rsd := DoblDobl_Series_Vector_Norms.Max_Norm(rsd);
      nrm_err := DoblDobl_Series_Vector_Norms.Max_Norm(err);
      put("Max norm of the backward error : "); put(nrm_rsd,3); new_line;
      put("Max norm of the forward error  : "); put(nrm_err,3); new_line;
    end if;
  end DoblDobl_Random_Linear_Solve;

  procedure QuadDobl_Random_Linear_Solve ( n,degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random matrix A and right hand side vector b
  --   of dimension n and solves the linear system A*x = b,
  --   in double double precision.
  --   The series are all of the given degree.

    use QuadDobl_Dense_Series_Vectors;
    use QuadDobl_Dense_Series_Matrices;
    use QuadDobl_Linear_Series_Solvers;

    A : constant QuadDobl_Dense_Series_Matrices.Matrix(1..n,1..n)
      := QuadDobl_Random_Series.Random_Series_Matrix(1,n,1,n,degree);
    wrk : QuadDobl_Dense_Series_Matrices.Matrix(1..n,1..n) := A;
    x : constant QuadDobl_Dense_Series_Vectors.Vector(1..n)
      := QuadDobl_Random_Series.Random_Series_Vector(1,n,degree);
    b : constant QuadDobl_Dense_Series_Vectors.Vector(1..n) := A*x;
    y : QuadDobl_Dense_Series_Vectors.Vector(1..n) := b;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    rsd,err : QuadDobl_Dense_Series_Vectors.Vector(1..n);
    nrm_rsd,nrm_err : quad_double;

  begin
    put_line("The coefficient matrix A :"); Write(A);
    put_line("The right hand side vector b :"); Write(b);
    LUfac(wrk,n,ipvt,info);
    put("The pivots : "); put(ipvt); new_line;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      LUsolve(wrk,n,ipvt,y);
      put_line("The generated solution x to A*x = b :"); Write(x);
      put_line("The computed solution x to A*x = b :"); Write(y);
      rsd := b - A*y; -- backward error
      put_line("The residual b - A*x :"); Write(rsd);
      err := x - y; -- forward error
      put_line("Difference between generated and computed :"); Write(err);
      nrm_rsd := QuadDobl_Series_Vector_Norms.Max_Norm(rsd);
      nrm_err := QuadDobl_Series_Vector_Norms.Max_Norm(err);
      put("Max norm of the backward error : "); put(nrm_rsd,3); new_line;
      put("Max norm of the forward error  : "); put(nrm_err,3); new_line;
    end if;
  end QuadDobl_Random_Linear_Solve;

  function Trunc ( A : Standard_Dense_Series_Matrices.Matrix )
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

  function Trunc ( A : DoblDobl_Dense_Series_Matrices.Matrix )
                 return DoblDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Takes the zero-th degree coefficient of every element in A.

    res : DoblDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := A(i,j).cff(0);
      end loop;
    end loop;
    return res;
  end Trunc;

  function Trunc ( A : QuadDobl_Dense_Series_Matrices.Matrix )
                 return QuadDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Takes the zero-th degree coefficient of every element in A.

    res : QuadDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));

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
  --   Shows the QR decomposition of the zero degree terms of A.

    B : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2)) := Trunc(A);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    qraux : Standard_Complex_Vectors.Vector(A'range(2));

  begin
    Standard_Complex_QR_Least_Squares.QRD(B,qraux,ipvt,false);
    put_line("QRD on zero degree terms:"); Write(B);
  end Zero_QRD;

  procedure Zero_QRD ( A : DoblDobl_Dense_Series_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Shows the QR decomposition of the zero degree terms of A.

    B : DoblDobl_Complex_Matrices.Matrix(A'range(1),A'range(2)) := Trunc(A);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    qraux : DoblDobl_Complex_Vectors.Vector(A'range(2));

  begin
    DoblDobl_Complex_QR_Least_Squares.QRD(B,qraux,ipvt,false);
    put_line("QRD on zero degree terms:"); Write(B);
  end Zero_QRD;

  procedure Zero_QRD ( A : QuadDobl_Dense_Series_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Shows the QR decomposition of the zero degree terms of A.

    B : QuadDobl_Complex_Matrices.Matrix(A'range(1),A'range(2)) := Trunc(A);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    qraux : QuadDobl_Complex_Vectors.Vector(A'range(2));

  begin
    QuadDobl_Complex_QR_Least_Squares.QRD(B,qraux,ipvt,false);
    put_line("QRD on zero degree terms:"); Write(B);
  end Zero_QRD;

  procedure Solve_Normal_Equations
              ( A : in Standard_Dense_Series_Matrices.Matrix;
                x,b : in Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   For the overdetermined system A*x = b, formulates the normal
  --   equations A^T A*x = A^T*b and then applies LU factorization
  --   to solve the system of normal equations,
  --   in standard double precision.
  --   The computed solution is compared to the constructed solution x
  --   and the residual of the normal equations is computed as well.

    use Standard_Dense_Series_Vectors;
    use Standard_Dense_Series_Matrices;
    use Standard_Linear_Series_Solvers;

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

  procedure Solve_Normal_Equations
              ( A : in DoblDobl_Dense_Series_Matrices.Matrix;
                x,b : in DoblDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   For the overdetermined system A*x = b, formulates the normal
  --   equations A^T A*x = A^T*b and then applies LU factorization
  --   to solve the system of normal equations,
  --   in double double precision.
  --   The computed solution is compared to the constructed solution x
  --   and the residual of the normal equations is computed as well.

    use DoblDobl_Dense_Series_Vectors;
    use DoblDobl_Dense_Series_Matrices;
    use DoblDobl_Linear_Series_Solvers;

    T : constant Matrix(A'range(2),A'range(1)) := Transpose(A);
    TA : constant Matrix(A'range(2),A'range(2)) := T*A;
    wrk : Matrix(A'range(2),A'range(2)) := TA;
    Tb : constant DoblDobl_Dense_Series_Vectors.Vector(A'range(2)) := T*b;
    n : constant integer32 := A'last(2);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    sol,rsd : DoblDobl_Dense_Series_Vectors.Vector(1..n);

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

  procedure Solve_Normal_Equations
              ( A : in QuadDobl_Dense_Series_Matrices.Matrix;
                x,b : in QuadDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   For the overdetermined system A*x = b, formulates the normal
  --   equations A^T A*x = A^T*b and then applies LU factorization
  --   to solve the system of normal equations,
  --   in quad double precision.
  --   The computed solution is compared to the constructed solution x
  --   and the residual of the normal equations is computed as well.

    use QuadDobl_Dense_Series_Vectors;
    use QuadDobl_Dense_Series_Matrices;
    use QuadDobl_Linear_Series_Solvers;

    T : constant Matrix(A'range(2),A'range(1)) := Transpose(A);
    TA : constant Matrix(A'range(2),A'range(2)) := T*A;
    wrk : Matrix(A'range(2),A'range(2)) := TA;
    Tb : constant QuadDobl_Dense_Series_Vectors.Vector(A'range(2)) := T*b;
    n : constant integer32 := A'last(2);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    sol,rsd : QuadDobl_Dense_Series_Vectors.Vector(1..n);

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
  --   this procedure computes the norm of every column in qr,
  --   in standard double precision.

    use Standard_Dense_Series;
    use Standard_Dense_Series_io;
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

  procedure Test_Normality
              ( qr : in DoblDobl_Dense_Series_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Given in qr is the orthogonal part of the QR decomposition,
  --   this procedure computes the norm of every column in qr,
  --   in double double precision.

    use DoblDobl_Dense_Series;
    use DoblDobl_Dense_Series_io;
    use DoblDobl_Dense_Series_Vectors;

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
      nrm := DoblDobl_Series_Vector_Norms.Norm(col);
      put_line("The norm :"); put(nrm);
    end loop;
  end Test_Normality;

  procedure Test_Normality
              ( qr : in QuadDobl_Dense_Series_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Given in qr is the orthogonal part of the QR decomposition,
  --   this procedure computes the norm of every column in qr,
  --   in quad double precision.

    use QuadDobl_Dense_Series;
    use QuadDobl_Dense_Series_io;
    use QuadDobl_Dense_Series_Vectors;

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
      nrm := QuadDobl_Series_Vector_Norms.Norm(col);
      put_line("The norm :"); put(nrm);
    end loop;
  end Test_Normality;

  procedure Test_Orthogonality
              ( qr : in Standard_Dense_Series_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Given in qr is the orthogonal part of the QR decomposition,
  --   this procedures computes all products of the vectors with
  --   all other vectors to test whether they are orthogonal.

    u,v : Standard_Dense_Series_Vectors.Vector(qr'range(1));
    ipr : Standard_Dense_Series.Series;

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
        ipr := Standard_Series_Vector_Norms.Inner_Product(u,v);
        put("Inner product of "); put(i,1); put(" with "); put(j,1);
        put_line(" :"); Standard_Dense_Series_io.put(ipr);
      end loop;
    end loop;
  end Test_Orthogonality;

  procedure Test_Orthogonality
              ( qr : in DoblDobl_Dense_Series_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Given in qr is the orthogonal part of the QR decomposition,
  --   this procedures computes all products of the vectors with
  --   all other vectors to test whether they are orthogonal.

    u,v : DoblDobl_Dense_Series_Vectors.Vector(qr'range(1));
    ipr : DoblDobl_Dense_Series.Series;

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
        ipr := DoblDobl_Series_Vector_Norms.Inner_Product(u,v);
        put("Inner product of "); put(i,1); put(" with "); put(j,1);
        put_line(" :"); DoblDobl_Dense_Series_io.put(ipr);
      end loop;
    end loop;
  end Test_Orthogonality;

  procedure Test_Orthogonality
              ( qr : in QuadDobl_Dense_Series_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Given in qr is the orthogonal part of the QR decomposition,
  --   this procedures computes all products of the vectors with
  --   all other vectors to test whether they are orthogonal.

    u,v : QuadDobl_Dense_Series_Vectors.Vector(qr'range(1));
    ipr : QuadDobl_Dense_Series.Series;

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
        ipr := QuadDobl_Series_Vector_Norms.Inner_Product(u,v);
        put("Inner product of "); put(i,1); put(" with "); put(j,1);
        put_line(" :"); QuadDobl_Dense_Series_io.put(ipr);
      end loop;
    end loop;
  end Test_Orthogonality;

  procedure Test_Basis
              ( wrk,A : in Standard_Dense_Series_Matrices.Matrix ) is

  -- DESCRIPION :
  --   Given in wrk the output of QRD on A, the orthogonal part of
  --   the QR decomposition is computed is tested for orthogonality
  --   and normality, in standard double precision.

    use Standard_Dense_Series_Matrices;

    qr : Matrix(wrk'range(1),wrk'range(2)) := wrk;
    ans : character;

  begin
    new_line;
    put_line("Computing the orthogonal part of the QR decomposition ...");
    Standard_Least_Squares_Series.Basis(qr,A);
    Test_Normality(qr);   
    new_line;
    put("Continue to the orthogonality test ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Test_Orthogonality(qr);
    end if;
  end Test_Basis;

  procedure Test_Basis
              ( wrk,A : in DoblDobl_Dense_Series_Matrices.Matrix ) is

  -- DESCRIPION :
  --   Given in wrk the output of QRD on A, the orthogonal part of
  --   the QR decomposition is computed is tested for orthogonality
  --   and normality, in double double precision.

    use DoblDobl_Dense_Series_Matrices;

    qr : Matrix(wrk'range(1),wrk'range(2)) := wrk;
    ans : character;

  begin
    new_line;
    put_line("Computing the orthogonal part of the QR decomposition ...");
    DoblDobl_Least_Squares_Series.Basis(qr,A);
    Test_Normality(qr);   
    new_line;
    put("Continue to the orthogonality test ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Test_Orthogonality(qr);
    end if;
  end Test_Basis;

  procedure Test_Basis
              ( wrk,A : in QuadDobl_Dense_Series_Matrices.Matrix ) is

  -- DESCRIPION :
  --   Given in wrk the output of QRD on A, the orthogonal part of
  --   the QR decomposition is computed is tested for orthogonality
  --   and normality, in double double precision.

    use QuadDobl_Dense_Series_Matrices;

    qr : Matrix(wrk'range(1),wrk'range(2)) := wrk;
    ans : character;

  begin
    new_line;
    put_line("Computing the orthogonal part of the QR decomposition ...");
    QuadDobl_Least_Squares_Series.Basis(qr,A);
    Test_Normality(qr);   
    new_line;
    put("Continue to the orthogonality test ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Test_Orthogonality(qr);
    end if;
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
    use Standard_Least_Squares_Series;

    wrk : Standard_Dense_Series_Matrices.Matrix(A'range(1),A'range(2)) := A;
    ipvt : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    qraux,sol : Standard_Dense_Series_Vectors.Vector(A'range(2));
    info : integer32;
    n : constant integer32 := A'last(1);
    m : constant integer32 := A'last(2);
    rsd,dum,dum2,dum3,err : Standard_Dense_Series_Vectors.Vector(1..n);
    ans : character;
    nrm_rsd,nrm_err : double_float;

  begin
    QRD(wrk,qraux,ipvt,false);
    put_line("The output of QRD :"); Write(wrk);
    new_line;
    put("View the QRD of the zero degree terms ? (y/n) ");
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
    new_line;
    put("Continue with the least squares solving ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      QRLS(wrk,n,m,qraux,b,dum2,dum3,sol,rsd,dum,110,info);
      if info /= 0 then
        put("info : "); put(info,1); new_line;
      else
        put_line("The constructed solution x to A*x = b :"); Write(x);
        put_line("The computed solution x to A*x = b :"); Write(sol);
        rsd := b - A*sol; -- backward error
        put_line("The residual b - A*x :"); Write(rsd);
        err := x - sol; -- forward error
        put("Difference between constructed and computed :"); Write(err);
        nrm_rsd := Standard_Series_Vector_Norms.Max_Norm(rsd);
        nrm_err := Standard_Series_Vector_Norms.Max_Norm(err);
        put("Max norm of backward error : "); put(nrm_rsd,3); new_line;
        put("Max norm of forward error  : "); put(nrm_err,3); new_line;
      end if;
    end if;
  end QR_Solve_Least_Squares;

  procedure QR_Solve_Least_Squares
              ( A : in DoblDobl_Dense_Series_Matrices.Matrix;
                x,b : in DoblDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Applies QR decomposition to the matrix A and then solves the
  --   system A*x = b with this QR decomposition.
  --   The computed solution is compared to the constructed solution
  --   and the residual is computed.

    use DoblDobl_Dense_Series_Vectors;
    use DoblDobl_Dense_Series_Matrices;
    use DoblDobl_Least_Squares_Series;

    wrk : DoblDobl_Dense_Series_Matrices.Matrix(A'range(1),A'range(2)) := A;
    ipvt : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    qraux,sol : DoblDobl_Dense_Series_Vectors.Vector(A'range(2));
    info : integer32;
    n : constant integer32 := A'last(1);
    m : constant integer32 := A'last(2);
    rsd,dum,dum2,dum3,err : DoblDobl_Dense_Series_Vectors.Vector(1..n);
    ans : character;
    nrm_rsd,nrm_err : double_double;

  begin
    QRD(wrk,qraux,ipvt,false);
    put_line("The output of QRD :"); Write(wrk);
    new_line;
    put("View the QRD of the zero degree terms ? (y/n) ");
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
    new_line;
    put("Continue with the least squares solving ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      QRLS(wrk,n,m,qraux,b,dum2,dum3,sol,rsd,dum,110,info);
      if info /= 0 then
        put("info : "); put(info,1); new_line;
      else
        put_line("The constructed solution x to A*x = b :"); Write(x);
        put_line("The computed solution x to A*x = b :"); Write(sol);
        rsd := b - A*sol; -- backward error
        put_line("The residual b - A*x :"); Write(rsd);
        err := x - sol; -- forward error
        put_line("Difference between constructed and computed :"); Write(err);
        nrm_rsd := DoblDobl_Series_Vector_Norms.Max_Norm(rsd);
        nrm_err := DoblDobl_Series_Vector_Norms.Max_Norm(err);
        put("Max norm of backward error : "); put(nrm_rsd,3); new_line;
        put("Max norm of forward error  : "); put(nrm_err,3); new_line;
      end if;
    end if;
  end QR_Solve_Least_Squares;

  procedure QR_Solve_Least_Squares
              ( A : in QuadDobl_Dense_Series_Matrices.Matrix;
                x,b : in QuadDobl_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Applies QR decomposition to the matrix A and then solves the
  --   system A*x = b with this QR decomposition.
  --   The computed solution is compared to the constructed solution
  --   and the residual is computed.

    use QuadDobl_Dense_Series_Vectors;
    use QuadDobl_Dense_Series_Matrices;
    use QuadDobl_Least_Squares_Series;

    wrk : QuadDobl_Dense_Series_Matrices.Matrix(A'range(1),A'range(2)) := A;
    ipvt : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    qraux,sol : QuadDobl_Dense_Series_Vectors.Vector(A'range(2));
    info : integer32;
    n : constant integer32 := A'last(1);
    m : constant integer32 := A'last(2);
    rsd,dum,dum2,dum3,err : QuadDobl_Dense_Series_Vectors.Vector(1..n);
    ans : character;
    nrm_rsd,nrm_err : quad_double;

  begin
    QRD(wrk,qraux,ipvt,false);
    put_line("The output of QRD :"); Write(wrk);
    new_line;
    put("View the QRD of the zero degree terms ? (y/n) ");
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
    new_line;
    put("Continue with the least squares solving ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      QRLS(wrk,n,m,qraux,b,dum2,dum3,sol,rsd,dum,110,info);
      if info /= 0 then
        put("info : "); put(info,1); new_line;
      else
        put_line("The constructed solution x to A*x = b :"); Write(x);
        put_line("The computed solution x to A*x = b :"); Write(sol);
        rsd := b - A*sol; -- backward error
        put_line("The residual b - A*x :"); Write(rsd);
        err := x - sol; -- forward error
        put_line("Difference between constructed and computed :"); Write(err);
        nrm_rsd := QuadDobl_Series_Vector_Norms.Max_Norm(rsd);
        nrm_err := QuadDobl_Series_Vector_Norms.Max_Norm(err);
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

    use Standard_Dense_Series_Matrices;

    A : constant Standard_Dense_Series_Matrices.Matrix(1..n,1..m)
      := Standard_Random_Series.Random_Series_Matrix(1,n,1,m,degree);
    x : constant Standard_Dense_Series_Vectors.Vector(1..m)
      := Standard_Random_Series.Random_Series_Vector(1,m,degree);
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
  end Standard_Random_Least_Squares;

  procedure DoblDobl_Random_Least_Squares ( n,m,degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-m matrix A with series of the given degree.
  --   A random solution x is generated and then the right hand side
  --   vector b is A*x.  for n > m this is an overdetermined linear
  --   system A*x = b, where the solution x is the generated vector.
  --   Computations are done in double double precision.

    use DoblDobl_Dense_Series_Matrices;

    A : constant DoblDobl_Dense_Series_Matrices.Matrix(1..n,1..m)
      := DoblDobl_Random_Series.Random_Series_Matrix(1,n,1,m,degree);
    x : constant DoblDobl_Dense_Series_Vectors.Vector(1..m)
      := DoblDobl_Random_Series.Random_Series_Vector(1,m,degree);
    b : constant DoblDobl_Dense_Series_Vectors.Vector(1..n) := A*x;
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
  end DoblDobl_Random_Least_Squares;

  procedure QuadDobl_Random_Least_Squares ( n,m,degree : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-m matrix A with series of the given degree.
  --   A random solution x is generated and then the right hand side
  --   vector b is A*x.  for n > m this is an overdetermined linear
  --   system A*x = b, where the solution x is the generated vector.
  --   Computations are done in quad double precision.

    use QuadDobl_Dense_Series_Matrices;

    A : constant QuadDobl_Dense_Series_Matrices.Matrix(1..n,1..m)
      := QuadDobl_Random_Series.Random_Series_Matrix(1,n,1,m,degree);
    x : constant QuadDobl_Dense_Series_Vectors.Vector(1..m)
      := QuadDobl_Random_Series.Random_Series_Vector(1,m,degree);
    b : constant QuadDobl_Dense_Series_Vectors.Vector(1..n) := A*x;
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
  end QuadDobl_Random_Least_Squares;

  function Prompt_for_Precision return character is

  -- DESCRIPTION :
  --   Displays the menu for the working precision,
  --   prompts for '0', '1', or '2', depending whether double,
  --   double double, or quad double precision is selected.
  --   Returns '0', '1', or '2'.

    ans : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    return ans;
  end Prompt_for_Precision;

  procedure Test_Solving is

  -- DESCRIPTION :
  --   Prompts the user for the dimension of the linear system
  --   and the degree of the series, solves a random system.

    dim,ord,nrows,ncols : integer32 := 0;
    ans,prc : character;

  begin
    new_line;
    put("Give the degree of the series : "); get(ord);
    new_line;
    put("Square system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give the dimension : "); get(dim);
      prc := Prompt_for_Precision;
      case prc is
        when '0' => Standard_Random_Linear_Solve(dim,ord);
        when '1' => DoblDobl_Random_Linear_Solve(dim,ord);
        when '2' => QuadDobl_Random_Linear_Solve(dim,ord);
        when others => null;
      end case;
    else
      put("Give number of rows : "); get(nrows);
      put("Give number of colums : "); get(ncols);
      prc := Prompt_for_Precision;
      case prc is
        when '0' => Standard_Random_Least_Squares(nrows,ncols,ord);
        when '1' => DoblDobl_Random_Least_Squares(nrows,ncols,ord);
        when '2' => QuadDobl_Random_Least_Squares(nrows,ncols,ord);
        when others => null;
      end case;
    end if;
  end Test_Solving;

  procedure Standard_Test_Conversion ( n,m,d : in integer32 ) is

  -- DESCRIPTION :
  --   Converts an n-by-m matrix of series of degree d with standard
  --   double precision complex coefficients into a matrix series.

    sA : constant Standard_Dense_Series_Matrices.Matrix(1..n,1..m)
       := Standard_Random_Series.Random_Series_Matrix(1,n,1,m,d);
    As : constant Standard_Dense_Matrix_Series.Matrix 
       := Standard_Dense_Matrix_Series.Create(sA); 

  begin
    Write(sA);
    put_line("The coefficients of the matrix series :");
    put(As);
  end Standard_Test_Conversion;

  procedure DoblDobl_Test_Conversion ( n,m,d : in integer32 ) is

  -- DESCRIPTION :
  --   Converts an n-by-m matrix of series of degree d with double
  --   double precision complex coefficients into a matrix series.

    sA : constant DoblDobl_Dense_Series_Matrices.Matrix(1..n,1..m)
       := DoblDobl_Random_Series.Random_Series_Matrix(1,n,1,m,d);
    As : constant DoblDobl_Dense_Matrix_Series.Matrix 
       := DoblDobl_Dense_Matrix_Series.Create(sA); 

  begin
    Write(sA);
    put_line("The coefficients of the matrix series :");
    put(As);
  end DoblDobl_Test_Conversion;

  procedure QuadDobl_Test_Conversion ( n,m,d : in integer32 ) is

  -- DESCRIPTION :
  --   Converts an n-by-m matrix of series of degree d with double
  --   double precision complex coefficients into a matrix series.

    sA : constant QuadDobl_Dense_Series_Matrices.Matrix(1..n,1..m)
       := QuadDobl_Random_Series.Random_Series_Matrix(1,n,1,m,d);
    As : constant QuadDobl_Dense_Matrix_Series.Matrix 
       := QuadDobl_Dense_Matrix_Series.Create(sA); 

  begin
    Write(sA);
    put_line("The coefficients of the matrix series :");
    put(As);
  end QuadDobl_Test_Conversion;

  procedure Test_Conversion is

  -- DESCRIPTION :
  --   Prompts for the number of rows, the number of columns,
  --   and the degree of the series.

    nbrows,nbcols,deg : integer32 := 0;
    prc : character;

  begin
    new_line;
    put_line("Converting random series matrix into matrix series ...");
    put("  Give the number of rows : "); get(nbrows);
    put("  Give the number of columns : "); get(nbcols);
    put("  Give the degree of the series : "); get(deg);
    prc := Prompt_for_Precision;
    case prc is
      when '0' => Standard_Test_Conversion(nbrows,nbcols,deg);
      when '1' => DoblDobl_Test_Conversion(nbrows,nbcols,deg);
      when '2' => QuadDobl_Test_Conversion(nbrows,nbcols,deg);
      when others => null;
    end case;
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
end ts_sermat;
