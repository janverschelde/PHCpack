with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Hexa_Double_Numbers;                 use Hexa_Double_Numbers;
with Hexa_Double_Numbers_io;              use Hexa_Double_Numbers_io;
with HexaDobl_Complex_Numbers_io;         use HexaDobl_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;         use Standard_Integer_Vectors_io;
with HexaDobl_Complex_Vectors;
with HexaDobl_Complex_QR_Least_Squares;
with HexaDobl_Complex_Series;
with HexaDobl_Complex_Series_io;
with HexaDobl_Complex_Matrix_Series;
with HexaDobl_Complex_Matrix_Series_io;   use HexaDobl_Complex_Matrix_Series_io;
with HexaDobl_Random_Series_Vectors;
with HexaDobl_Random_Series_Matrices;
with HexaDobl_CSeries_Vector_Norms;
with HexaDobl_Series_Linear_Solvers;
with HexaDobl_Series_Least_Squares;

package body Test_HexaDobl_Matrix_Series is

  procedure Write ( v : in HexaDobl_Complex_Series_Vectors.Vector ) is
  begin
    for i in v'range loop
      put("Component "); put(i,1); put_line(" :");
      HexaDobl_Complex_Series_io.put(v(i));
    end loop;
  end Write;

  procedure Write ( A : in HexaDobl_Complex_Series_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("Component ");
        put(i,1); put(", "); put(j,1); put_line(" :");
        HexaDobl_Complex_Series_io.put(A(i,j));
      end loop;
    end loop;
  end Write;

  procedure Write ( A : in HexaDobl_Complex_Matrices.Matrix ) is
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("Component ");
        put(i,1); put(", "); put(j,1); put_line(" :");
        put(A(i,j)); new_line;
      end loop;
    end loop;
  end Write;

  procedure HexaDobl_Random_Linear_Solve ( n,degree : in integer32 ) is

    use HexaDobl_Complex_Series_Vectors;
    use HexaDobl_Complex_Series_Matrices;
    use HexaDobl_Series_Linear_Solvers;

    A : constant HexaDobl_Complex_Series_Matrices.Matrix(1..n,1..n)
      := HexaDobl_Random_Series_Matrices.Random_Series_Matrix(1,n,1,n,degree);
    wrk : HexaDobl_Complex_Series_Matrices.Matrix(1..n,1..n);
    x : constant HexaDobl_Complex_Series_Vectors.Vector(1..n)
      := HexaDobl_Random_Series_Vectors.Random_Series_Vector(1,n,degree);
    b : constant HexaDobl_Complex_Series_Vectors.Vector(1..n) := A*x;
    y : HexaDobl_Complex_Series_Vectors.Vector(1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    rsd,err : HexaDobl_Complex_Series_Vectors.Vector(1..n);
    nrm_rsd,nrm_err : hexa_double;

  begin
    put_line("The coefficient matrix A :"); Write(A);
    put_line("The right hand side vector b :"); Write(b);
    HexaDobl_Complex_Series_Matrices.Copy(A,wrk); -- deep copy
    LUfac(wrk,n,ipvt,info);
    put("The pivots : "); put(ipvt); new_line;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      HexaDobl_Complex_Series_Vectors.Copy(b,y); -- deep copy
      LUsolve(wrk,n,ipvt,y);
      put_line("The generated solution x to A*x = b :"); Write(x);
      put_line("The computed solution x to A*x = b :"); Write(y);
      rsd := b - A*y; -- backward error
      put_line("The residual b - A*x :"); Write(rsd);
      err := x - y; -- forward error
      put_line("The difference between generated and computed :"); Write(err);
      nrm_rsd := HexaDobl_CSeries_Vector_Norms.Max_Norm(rsd);
      nrm_err := HexaDobl_CSeries_Vector_Norms.Max_Norm(err);
      put("Max norm of the backward error : "); put(nrm_rsd,3); new_line;
      put("Max norm of the forward error  : "); put(nrm_err,3); new_line;
    end if;
  end HexaDobl_Random_Linear_Solve;

  function Trunc ( A : HexaDobl_Complex_Series_Matrices.Matrix )
                 return HexaDobl_Complex_Matrices.Matrix is

    res : HexaDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := A(i,j).cff(0);
      end loop;
    end loop;
    return res;
  end Trunc;

  procedure Zero_QRD ( A : HexaDobl_Complex_Series_Matrices.Matrix ) is

    B : HexaDobl_Complex_Matrices.Matrix(A'range(1),A'range(2)) := Trunc(A);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    qraux : HexaDobl_Complex_Vectors.Vector(A'range(2));

  begin
    HexaDobl_Complex_QR_Least_Squares.QRD(B,qraux,ipvt,false);
    put_line("QRD on zero degree terms:"); Write(B);
  end Zero_QRD;

  procedure Solve_Normal_Equations
              ( A : in HexaDobl_Complex_Series_Matrices.Matrix;
                x,b : in HexaDobl_Complex_Series_Vectors.Vector ) is

    use HexaDobl_Complex_Series_Vectors;
    use HexaDobl_Complex_Series_Matrices;
    use HexaDobl_Series_Linear_Solvers;

    T : constant Matrix(A'range(2),A'range(1)) := Transpose(A);
    TA : constant Matrix(A'range(2),A'range(2)) := T*A;
    wrk : Matrix(A'range(2),A'range(2));
    Tb : constant HexaDobl_Complex_Series_Vectors.Vector(A'range(2)) := T*b;
    n : constant integer32 := A'last(2);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    sol,rsd : HexaDobl_Complex_Series_Vectors.Vector(1..n);

  begin
    HexaDobl_Complex_Series_Matrices.Copy(TA,wrk);
    LUfac(wrk,n,ipvt,info);
    put("The pivots : "); put(ipvt); new_line;
    if info /= 0 then
      put("info : "); put(info,1); new_line;
    else
      HexaDobl_Complex_Series_Vectors.Copy(Tb,sol);
      LUsolve(wrk,n,ipvt,sol);
      put_line("The constructed solution x to A*x = b :"); Write(x);
      put_line("The computed solution x to A*x = b :"); Write(sol);
      rsd := Tb - TA*sol;
      put_line("The residual b - A*x :"); Write(rsd);
    end if;
  end Solve_Normal_Equations;

  procedure Test_Normality
              ( qr : in HexaDobl_Complex_Series_Matrices.Matrix ) is

    use HexaDobl_Complex_Series;
    use HexaDobl_Complex_Series_io;
    use HexaDobl_Complex_Series_Vectors;

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
      nrm := HexaDobl_CSeries_Vector_Norms.Norm(col);
      put_line("The norm :"); put(nrm);
    end loop;
  end Test_Normality;

  procedure Test_Orthogonality
              ( qr : in HexaDobl_Complex_Series_Matrices.Matrix ) is

    u,v : HexaDobl_Complex_Series_Vectors.Vector(qr'range(1));
    deg : constant integer32 := qr(qr'first(1),qr'first(2)).deg;
    ipr : HexaDobl_Complex_Series.Series(deg);

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
        ipr := HexaDobl_CSeries_Vector_Norms.Inner_Product(u,v);
        put("Inner product of "); put(i,1); put(" with "); put(j,1);
        put_line(" :"); HexaDobl_Complex_Series_io.put(ipr);
      end loop;
    end loop;
  end Test_Orthogonality;

  procedure Test_Basis
              ( wrk,A : in HexaDobl_Complex_Series_Matrices.Matrix ) is

    use HexaDobl_Complex_Series_Matrices;

    qr : Matrix(wrk'range(1),wrk'range(2));
    ans : character;

  begin
    HexaDobl_Complex_Series_Matrices.Copy(wrk,qr);
    new_line;
    put_line("Computing the orthogonal part of the QR decomposition ...");
    HexaDobl_Series_Least_Squares.Basis(qr,A);
    Test_Normality(qr);   
    new_line;
    put("Continue to the orthogonality test ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Test_Orthogonality(qr);
    end if;
  end Test_Basis;

  procedure QR_Solve_Least_Squares
              ( A : in HexaDobl_Complex_Series_Matrices.Matrix;
                x,b : in HexaDobl_Complex_Series_Vectors.Vector ) is

    use HexaDobl_Complex_Series;
    use HexaDobl_Complex_Series_Vectors;
    use HexaDobl_Complex_Series_Matrices;
    use HexaDobl_Series_Least_Squares;

    deg : constant integer32 := A(A'first(1),A'first(2)).deg;
    wrk : HexaDobl_Complex_Series_Matrices.Matrix(A'range(1),A'range(2));
    ipvt : Standard_Integer_Vectors.Vector(A'range(2)) := (A'range(2) => 0);
    qraux,sol : HexaDobl_Complex_Series_Vectors.Vector(A'range(2));
    info : integer32;
    n : constant integer32 := A'last(1);
    m : constant integer32 := A'last(2);
    rsd,dum,dum2,dum3 : HexaDobl_Complex_Series_Vectors.Vector(1..n);
    err : HexaDobl_Complex_Series_Vectors.Vector(1..m);
    ans : character;
    nrm_rsd,nrm_err : hexa_double;

  begin
    HexaDobl_Complex_Series_Matrices.Copy(A,wrk);
    for i in qraux'range loop
      qraux(i) := new Series'(HexaDobl_Complex_Series.Create(0,deg));
      sol(i) := new Series'(HexaDobl_Complex_Series.Create(0,deg));
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
        nrm_rsd := HexaDobl_CSeries_Vector_Norms.Max_Norm(rsd);
        nrm_err := HexaDobl_CSeries_Vector_Norms.Max_Norm(err);
        put("Max norm of backward error : "); put(nrm_rsd,3); new_line;
        put("Max norm of forward error  : "); put(nrm_err,3); new_line;
      end if;
    end if;
  end QR_Solve_Least_Squares;

  procedure HexaDobl_Random_Least_Squares ( n,m,degree : in integer32 ) is

    use HexaDobl_Complex_Series_Matrices;

    A : constant HexaDobl_Complex_Series_Matrices.Matrix(1..n,1..m)
      := HexaDobl_Random_Series_Matrices.Random_Series_Matrix(1,n,1,m,degree);
    x : constant HexaDobl_Complex_Series_Vectors.Vector(1..m)
      := HexaDobl_Random_Series_Vectors.Random_Series_Vector(1,m,degree);
    b : constant HexaDobl_Complex_Series_Vectors.Vector(1..n) := A*x;
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
  end HexaDobl_Random_Least_Squares;

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
      HexaDobl_Random_Linear_Solve(dim,deg);
    else
      put("Give number of rows : "); get(nrows);
      put("Give number of colums : "); get(ncols);
      HexaDobl_Random_Least_Squares(nrows,ncols,deg);
    end if;
  end Test_Solving;

  procedure HexaDobl_Test_Conversion ( n,m,d : in integer32 ) is

    sA : constant HexaDobl_Complex_Series_Matrices.Matrix(1..n,1..m)
       := HexaDobl_Random_Series_Matrices.Random_Series_Matrix(1,n,1,m,d);
    As : constant HexaDobl_Complex_Matrix_Series.Matrix 
       := HexaDobl_Complex_Matrix_Series.Create(sA); 

  begin
    Write(sA);
    put_line("The coefficients of the matrix series :");
    put(As);
  end HexaDobl_Test_Conversion;

  procedure Test_Conversion is

    nbrows,nbcols,deg : integer32 := 0;

  begin
    put_line("Converting random series matrix into matrix series ...");
    put("  Give the number of rows : "); get(nbrows);
    put("  Give the number of columns : "); get(nbcols);
    put("  Give the degree of the series : "); get(deg);
    HexaDobl_Test_Conversion(nbrows,nbcols,deg);
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

end Test_HexaDobl_Matrix_Series;
