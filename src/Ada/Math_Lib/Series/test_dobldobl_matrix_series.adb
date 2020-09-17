with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Double_Double_Numbers;               use Double_Double_Numbers;
with Double_Double_Numbers_io;            use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers_io;         use DoblDobl_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;         use Standard_Integer_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_QR_Least_Squares;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_io;
with DoblDobl_Complex_Matrix_Series;
with DoblDobl_Complex_Matrix_Series_io;   use DoblDobl_Complex_Matrix_Series_io;
with DoblDobl_Random_Series_Vectors;
with DoblDobl_Random_Series_Matrices;
with DoblDobl_CSeries_Vector_Norms;
with DoblDobl_Series_Linear_Solvers;
with DoblDobl_Series_Least_Squares;

package body Test_DoblDobl_Matrix_Series is

  procedure Write ( v : in DoblDobl_Complex_Series_Vectors.Vector ) is
  begin
    for i in v'range loop
      put("Component "); put(i,1); put_line(" :");
      DoblDobl_Complex_Series_io.put(v(i));
    end loop;
  end Write;

  procedure Write ( A : in DoblDobl_Complex_Series_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("Component ");
        put(i,1); put(", "); put(j,1); put_line(" :");
        DoblDobl_Complex_Series_io.put(A(i,j));
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

  procedure DoblDobl_Random_Linear_Solve ( n,degree : in integer32 ) is

    use DoblDobl_Complex_Series_Vectors;
    use DoblDobl_Complex_Series_Matrices;
    use DoblDobl_Series_Linear_Solvers;

    A : constant DoblDobl_Complex_Series_Matrices.Matrix(1..n,1..n)
      := DoblDobl_Random_Series_Matrices.Random_Series_Matrix(1,n,1,n,degree);
    wrk : DoblDobl_Complex_Series_Matrices.Matrix(1..n,1..n);
    x : constant DoblDobl_Complex_Series_Vectors.Vector(1..n)
      := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,n,degree);
    b : constant DoblDobl_Complex_Series_Vectors.Vector(1..n) := A*x;
    y : DoblDobl_Complex_Series_Vectors.Vector(1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    rsd,err : DoblDobl_Complex_Series_Vectors.Vector(1..n);
    nrm_rsd,nrm_err : double_double;

  begin
    put_line("The coefficient matrix A :"); Write(A);
    put_line("The right hand side vector b :"); Write(b);
    DoblDobl_Complex_Series_Matrices.Copy(A,wrk); -- deep copy
    LUfac(wrk,n,ipvt,info);
    put("The pivots : "); put(ipvt); new_line;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      DoblDobl_Complex_Series_Vectors.Copy(b,y); -- deep copy
      LUsolve(wrk,n,ipvt,y);
      put_line("The generated solution x to A*x = b :"); Write(x);
      put_line("The computed solution x to A*x = b :"); Write(y);
      rsd := b - A*y; -- backward error
      put_line("The residual b - A*x :"); Write(rsd);
      err := x - y; -- forward error
      put_line("The difference between generated and computed :"); Write(err);
      nrm_rsd := DoblDobl_CSeries_Vector_Norms.Max_Norm(rsd);
      nrm_err := DoblDobl_CSeries_Vector_Norms.Max_Norm(err);
      put("Max norm of the backward error : "); put(nrm_rsd,3); new_line;
      put("Max norm of the forward error  : "); put(nrm_err,3); new_line;
    end if;
  end DoblDobl_Random_Linear_Solve;

  function Trunc ( A : DoblDobl_Complex_Series_Matrices.Matrix )
                 return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := A(i,j).cff(0);
      end loop;
    end loop;
    return res;
  end Trunc;

  procedure Zero_QRD ( A : DoblDobl_Complex_Series_Matrices.Matrix ) is

    B : DoblDobl_Complex_Matrices.Matrix(A'range(1),A'range(2)) := Trunc(A);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    qraux : DoblDobl_Complex_Vectors.Vector(A'range(2));

  begin
    DoblDobl_Complex_QR_Least_Squares.QRD(B,qraux,ipvt,false);
    put_line("QRD on zero degree terms:"); Write(B);
  end Zero_QRD;

  procedure Solve_Normal_Equations
              ( A : in DoblDobl_Complex_Series_Matrices.Matrix;
                x,b : in DoblDobl_Complex_Series_Vectors.Vector ) is

    use DoblDobl_Complex_Series_Vectors;
    use DoblDobl_Complex_Series_Matrices;
    use DoblDobl_Series_Linear_Solvers;

    T : constant Matrix(A'range(2),A'range(1)) := Transpose(A);
    TA : constant Matrix(A'range(2),A'range(2)) := T*A;
    wrk : Matrix(A'range(2),A'range(2));
    Tb : constant DoblDobl_Complex_Series_Vectors.Vector(A'range(2)) := T*b;
    n : constant integer32 := A'last(2);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    sol,rsd : DoblDobl_Complex_Series_Vectors.Vector(1..n);

  begin
    DoblDobl_Complex_Series_Matrices.Copy(TA,wrk);
    LUfac(wrk,n,ipvt,info);
    put("The pivots : "); put(ipvt); new_line;
    if info /= 0 then
      put("info : "); put(info,1); new_line;
    else
      DoblDobl_Complex_Series_Vectors.Copy(Tb,sol);
      LUsolve(wrk,n,ipvt,sol);
      put_line("The constructed solution x to A*x = b :"); Write(x);
      put_line("The computed solution x to A*x = b :"); Write(sol);
      rsd := Tb - TA*sol;
      put_line("The residual b - A*x :"); Write(rsd);
    end if;
  end Solve_Normal_Equations;

  procedure Test_Normality
              ( qr : in DoblDobl_Complex_Series_Matrices.Matrix ) is

    use DoblDobl_Complex_Series;
    use DoblDobl_Complex_Series_io;
    use DoblDobl_Complex_Series_Vectors;

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
      nrm := DoblDobl_CSeries_Vector_Norms.Norm(col);
      put_line("The norm :"); put(nrm);
    end loop;
  end Test_Normality;

  procedure Test_Orthogonality
              ( qr : in DoblDobl_Complex_Series_Matrices.Matrix ) is

    u,v : DoblDobl_Complex_Series_Vectors.Vector(qr'range(1));
    deg : constant integer32 := qr(qr'first(1),qr'first(2)).deg;
    ipr : DoblDobl_Complex_Series.Series(deg);

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
        ipr := DoblDobl_CSeries_Vector_Norms.Inner_Product(u,v);
        put("Inner product of "); put(i,1); put(" with "); put(j,1);
        put_line(" :"); DoblDobl_Complex_Series_io.put(ipr);
      end loop;
    end loop;
  end Test_Orthogonality;

  procedure Test_Basis
              ( wrk,A : in DoblDobl_Complex_Series_Matrices.Matrix ) is

    use DoblDobl_Complex_Series_Matrices;

    qr : Matrix(wrk'range(1),wrk'range(2));
    ans : character;

  begin
    DoblDobl_Complex_Series_Matrices.Copy(wrk,qr);
    new_line;
    put_line("Computing the orthogonal part of the QR decomposition ...");
    DoblDobl_Series_Least_Squares.Basis(qr,A);
    Test_Normality(qr);   
    new_line;
    put("Continue to the orthogonality test ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Test_Orthogonality(qr);
    end if;
  end Test_Basis;

  procedure QR_Solve_Least_Squares
              ( A : in DoblDobl_Complex_Series_Matrices.Matrix;
                x,b : in DoblDobl_Complex_Series_Vectors.Vector ) is

    use DoblDobl_Complex_Series;
    use DoblDobl_Complex_Series_Vectors;
    use DoblDobl_Complex_Series_Matrices;
    use DoblDobl_Series_Least_Squares;

    deg : constant integer32 := A(A'first(1),A'first(2)).deg;
    wrk : DoblDobl_Complex_Series_Matrices.Matrix(A'range(1),A'range(2));
    ipvt : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    qraux,sol : DoblDobl_Complex_Series_Vectors.Vector(A'range(2));
    info : integer32;
    n : constant integer32 := A'last(1);
    m : constant integer32 := A'last(2);
    rsd,dum,dum2,dum3 : DoblDobl_Complex_Series_Vectors.Vector(1..n);
    err : DoblDobl_Complex_Series_Vectors.Vector(1..m);
    ans : character;
    nrm_rsd,nrm_err : double_double;

  begin
    DoblDobl_Complex_Series_Matrices.Copy(A,wrk);
    for i in qraux'range loop
      qraux(i) := new Series'(DoblDobl_Complex_Series.Create(0,deg));
      sol(i) := new Series'(DoblDobl_Complex_Series.Create(0,deg));
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
        nrm_rsd := DoblDobl_CSeries_Vector_Norms.Max_Norm(rsd);
        nrm_err := DoblDobl_CSeries_Vector_Norms.Max_Norm(err);
        put("Max norm of backward error : "); put(nrm_rsd,3); new_line;
        put("Max norm of forward error  : "); put(nrm_err,3); new_line;
      end if;
    end if;
  end QR_Solve_Least_Squares;

  procedure DoblDobl_Random_Least_Squares ( n,m,degree : in integer32 ) is

    use DoblDobl_Complex_Series_Matrices;

    A : constant DoblDobl_Complex_Series_Matrices.Matrix(1..n,1..m)
      := DoblDobl_Random_Series_Matrices.Random_Series_Matrix(1,n,1,m,degree);
    x : constant DoblDobl_Complex_Series_Vectors.Vector(1..m)
      := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,m,degree);
    b : constant DoblDobl_Complex_Series_Vectors.Vector(1..n) := A*x;
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

  procedure Test_Solving is

    dim,deg,nrows,ncols : integer32 := 0;
    ans : character;

  begin
    put("Give the degree of the series : "); get(deg);
    new_line;
    put("Square system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give the dimension : "); get(dim);
      DoblDobl_Random_Linear_Solve(dim,deg);
    else
      put("Give number of rows : "); get(nrows);
      put("Give number of colums : "); get(ncols);
      DoblDobl_Random_Least_Squares(nrows,ncols,deg);
    end if;
  end Test_Solving;

  procedure DoblDobl_Test_Conversion ( n,m,d : in integer32 ) is

    sA : constant DoblDobl_Complex_Series_Matrices.Matrix(1..n,1..m)
       := DoblDobl_Random_Series_Matrices.Random_Series_Matrix(1,n,1,m,d);
    As : constant DoblDobl_Complex_Matrix_Series.Matrix 
       := DoblDobl_Complex_Matrix_Series.Create(sA); 

  begin
    Write(sA);
    put_line("The coefficients of the matrix series :");
    put(As);
  end DoblDobl_Test_Conversion;

  procedure Test_Conversion is

    nbrows,nbcols,deg : integer32 := 0;

  begin
    new_line;
    put_line("Converting random series matrix into matrix series ...");
    put("  Give the number of rows : "); get(nbrows);
    put("  Give the number of columns : "); get(nbcols);
    put("  Give the degree of the series : "); get(deg);
    DoblDobl_Test_Conversion(nbrows,nbcols,deg);
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
    new_line;
    case ans is
      when '0' => Test_Solving;
      when '1' => Test_Conversion;
      when others => null;
    end case;
  end Main;

end Test_DoblDobl_Matrix_Series;
