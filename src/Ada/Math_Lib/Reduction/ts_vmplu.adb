with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with DoblDobl_Mathematical_Functions;   use DoblDobl_Mathematical_Functions;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;     use Standard_Floating_Matrices_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Double_Double_Vectors;
with Double_Double_Vectors_io;          use Double_Double_Vectors_io;
with Double_Double_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Matrices;
with Standard_Random_Matrices;          use Standard_Random_Matrices;
with Standard_Floating_Linear_Solvers;
with Standard_Complex_Linear_Solvers;
with Double_Double_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;

procedure ts_vmplu is

-- DESCRIPTION :
--   Test on variable precision linear system solving with LU decomposition,
--   with real and complex standard double and double double numbers.

  function d2dd ( mat : Standard_Floating_Matrices.Matrix )
                return Double_Double_Matrices.Matrix is

    res : Double_Double_Matrices.Matrix(mat'range(1),mat'range(2));

  begin
    for i in mat'range(1) loop
      for j in mat'range(2) loop
        res(i,j) := create(mat(i,j));
      end loop;
    end loop;
    return res;
  end d2dd;

  function d2dd ( mat : Standard_Complex_Matrices.Matrix )
                return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix(mat'range(1),mat'range(2));

  begin
    for i in mat'range(1) loop
      for j in mat'range(2) loop
        declare
          strp : constant double_float
               := Standard_Complex_Numbers.REAL_PART(mat(i,j));
          ddrp : constant double_double := create(strp);
          stip : constant double_float
               := Standard_Complex_Numbers.IMAG_PART(mat(i,j));
          ddip : constant double_double := create(stip);
        begin
          res(i,j) := DoblDobl_Complex_Numbers.create(ddrp,ddip);
        end;
      end loop;
    end loop;
    return res;
  end d2dd;

  function Singular_Value_Matrix
             ( n : integer32; c : double_float )
             return Standard_Floating_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns an n-dimensional diagonal matrix that has on its diagonal 
  --   the values 1.0 + (c - 1.0)/(n-1)*(i - 1), for i in 1..n. 

   res : Standard_Floating_Matrices.Matrix(1..n,1..n);
   dm1 : constant double_float := double_float(n) - 1.0;

  begin
    for i in 1..n loop
      for j in 1..n loop
        res(i,j) := 0.0;
      end loop;
      res(i,i) := 1.0 + (c - 1.0)/dm1*(double_float(i) - 1.0);
    end loop;
   -- put_line("The diagonal matrix :"); put(res);
    return res;
  end Singular_Value_Matrix;

  function Singular_Value_Matrix
             ( n : integer32; c : double_float )
             return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns an n-dimensional diagonal matrix that has on its diagonal 
  --   the values 1.0 + (c - 1.0)/(n-1)*(i - 1), for i in 1..n. 

   res : Standard_Complex_Matrices.Matrix(1..n,1..n);
   dm1 : constant double_float := double_float(n) - 1.0;
   rdi : double_float;

  begin
    for i in 1..n loop
      for j in 1..n loop
        res(i,j) := Standard_Complex_Numbers.Create(0.0);
      end loop;
      rdi := 1.0 + (c - 1.0)/dm1*(double_float(i) - 1.0);
      res(i,i) := Standard_Complex_Numbers.Create(rdi);
    end loop;
    return res;
  end Singular_Value_Matrix;

  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return Standard_Floating_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns a random n-dimensional matrix with condition number c.

    res : Standard_Floating_Matrices.Matrix(1..n,1..n);
    svm : constant Standard_Floating_Matrices.Matrix(1..n,1..n)
        := Singular_Value_Matrix(n,c);
    rq1 : constant Standard_Floating_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 
    rq2 : constant Standard_Floating_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 

    use Standard_Floating_Matrices;

  begin
    res := rq1*svm*rq2;
    return res;
  end Random_Conditioned_Matrix;

  function Random_Conditioned_Matrix
             ( n : integer32; c : double_float )
             return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns a random n-dimensional matrix with condition number c.

    res : Standard_Complex_Matrices.Matrix(1..n,1..n);
    svm : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
        := Singular_Value_Matrix(n,c);
    rq1 : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 
    rq2 : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
        := Random_Orthogonal_Matrix(natural32(n),natural32(n)); 

    use Standard_Complex_Matrices;

  begin
    res := rq1*svm*rq2;
    return res;
  end Random_Conditioned_Matrix;

  function Estimated_Loss_in_Fixed_Precision
             ( mat : in Standard_Floating_Matrices.Matrix )
             return integer32 is

  -- DESCRIPTION :
  --   Returns the estimated loss of decimal places that may occur
  --   when solving a linear system with coefficient matrix mat.
  --   As the estimate is computed in fixed double precision,
  --   the result is not reliable for matrices with condition numbers
  --   that are larger than 1.0E+15.

  -- REQUIRED : mat'range(1) = mat'range(2).

    res : integer32;
    dim : constant integer32 := mat'last(1);
    wrk : Standard_Floating_Matrices.Matrix(mat'range(1),mat'range(2)) := mat;
    piv : Standard_Integer_Vectors.Vector(1..dim);
    rco : double_float;

  begin
    Standard_Floating_Linear_Solvers.lufco(wrk,dim,piv,rco);
    put("The value of rco : "); put(rco); new_line;
    put("Estimated condition number :");
    Standard_Floating_Numbers_io.put(1.0/rco,3); new_line;
    res := integer32(log10(rco));
    return res;
  end Estimated_Loss_in_Fixed_Precision;

  function Estimated_Loss_in_Fixed_Precision
             ( mat : in Standard_Complex_Matrices.Matrix )
             return integer32 is

  -- DESCRIPTION :
  --   Returns the estimated loss of decimal places that may occur
  --   when solving a linear system with coefficient matrix mat.
  --   As the estimate is computed in fixed double precision,
  --   the result is not reliable for matrices with condition numbers
  --   that are larger than 1.0E+15.

  -- REQUIRED : mat'range(1) = mat'range(2).

    res : integer32;
    dim : constant integer32 := mat'last(1);
    wrk : Standard_Complex_Matrices.Matrix(mat'range(1),mat'range(2)) := mat;
    piv : Standard_Integer_Vectors.Vector(1..dim);
    rco : double_float;

  begin
    Standard_Complex_Linear_Solvers.lufco(wrk,dim,piv,rco);
    put("The value of rco : "); put(rco); new_line;
    put("Estimated condition number :");
    Standard_Floating_Numbers_io.put(1.0/rco,3); new_line;
    res := integer32(log10(rco));
    return res;
  end Estimated_Loss_in_Fixed_Precision;

  function Estimated_Loss_in_Fixed_Precision
             ( mat : in Double_Double_Matrices.Matrix )
             return integer32 is

  -- DESCRIPTION :
  --   Returns the estimated loss of decimal places that may occur
  --   when solving a linear system with coefficient matrix mat.

  -- REQUIRED : mat'range(1) = mat'range(2).

    res : integer32;
    dim : constant integer32 := mat'last(1);
    wrk : Double_Double_Matrices.Matrix(mat'range(1),mat'range(2)) := mat;
    piv : Standard_Integer_Vectors.Vector(1..dim);
    one : constant double_double := create(1.0);
    rco,condnumb : double_double;

  begin
    Double_Double_Linear_Solvers.lufco(wrk,dim,piv,rco);
    put("The value of rco : "); put(rco); new_line;
    condnumb := one/rco;
    put("Estimated condition number :  "); put(condnumb,3); new_line;
    res := integer32(to_double(log10(rco)));
    return res;
  end Estimated_Loss_in_Fixed_Precision;

  function Estimated_Loss_in_Fixed_Precision
             ( mat : in DoblDobl_Complex_Matrices.Matrix )
             return integer32 is

  -- DESCRIPTION :
  --   Returns the estimated loss of decimal places that may occur
  --   when solving a linear system with coefficient matrix mat.
  --   As the estimate is computed in fixed double precision,
  --   the result is not reliable for matrices with condition numbers
  --   that are larger than 1.0E+31.

  -- REQUIRED : mat'range(1) = mat'range(2).

    res : integer32;
    dim : constant integer32 := mat'last(1);
    wrk : DoblDobl_Complex_Matrices.Matrix(mat'range(1),mat'range(2)) := mat;
    piv : Standard_Integer_Vectors.Vector(1..dim);
    one : constant double_double := create(1.0);
    rco,condnumb : double_double;

  begin
    DoblDobl_Complex_Linear_Solvers.lufco(wrk,dim,piv,rco);
    put("The value of rco : "); put(rco); new_line;
    condnumb := one/rco;
    put("Estimated condition number :  "); put(condnumb,3); new_line;
    res := integer32(to_double(log10(rco)));
    return res;
  end Estimated_Loss_in_Fixed_Precision;

  function Estimated_Loss
             ( mat : in Standard_Floating_Matrices.Matrix )
             return integer32 is

  -- DESCRIPTION :
  --   Returns the estimated loss of decimal places that may occur
  --   when solving a linear system with coefficient matrix mat.

    res : integer32 := Estimated_Loss_in_Fixed_Precision(mat);

  begin
     if res < -15 then
       declare
         ddm : constant Double_Double_Matrices.Matrix
                 (mat'range(1),mat'range(2)) := d2dd(mat);
       begin
         res := Estimated_Loss_in_Fixed_Precision(ddm);
       end;
     end if;
    return res;
  end Estimated_Loss;

  function Estimated_Loss
             ( mat : in Standard_Complex_Matrices.Matrix )
             return integer32 is

  -- DESCRIPTION :
  --   Returns the estimated loss of decimal places that may occur
  --   when solving a linear system with coefficient matrix mat.

    res : integer32 := Estimated_Loss_in_Fixed_Precision(mat);

  begin
     if res < -15 then
       declare
         ddm : constant DoblDobl_Complex_Matrices.Matrix
                 (mat'range(1),mat'range(2)) := d2dd(mat);
       begin
         put_line("estimating with double double arithmetic ...");
         res := Estimated_Loss_in_Fixed_Precision(ddm);
       end;
     end if;
    return res;
  end Estimated_Loss;

  procedure Standard_Real_Solve
              ( dim,want_dcp : in integer32;
                mat : in Standard_Floating_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Makes a righthandside vector b so x = (1,1,..,1) is the exact solution.
  --   Solves the linear system in standard double floating-point arithmetic.
  --   Verifies whether the computed solution has indeed at least as many 
  --   decimal places correct as the value of want_dcp.

    x : constant Standard_Floating_Vectors.Vector(1..dim) := (1..dim => 1.0);
    b : Standard_Floating_Vectors.Vector(1..dim);
    wrk : Standard_Floating_Matrices.Matrix(1..dim,1..dim) := mat;
    piv : Standard_Integer_Vectors.Vector(1..dim);
    inf : integer32;
    tol : constant double_float := 10.0**(integer(-want_dcp - 1));
    okay : boolean := true;
    adf : double_float;
  
    use Standard_Floating_Matrices;

  begin
    b := mat*x;
    Standard_Floating_Linear_Solvers.lufac(wrk,dim,piv,inf);
    Standard_Floating_Linear_Solvers.lusolve(wrk,dim,piv,b);
    put_line("The computed solution :"); put_line(b);
    for i in b'range loop
      adf := abs(b(i) - x(i));
      if adf > tol then
        put("x("); put(i,1); put(") = "); put(b(i));
        put_line(" is off!"); okay := false;
      end if;
      exit when not okay;
    end loop;
    if okay then
      put("Solution accurate with at least ");
      put(want_dcp,1); put_line(" decimal places.");
    else
      put("Accuracy of "); put(want_dcp,1); 
      put_line(" decimal places is not obtained.");
    end if;
  end Standard_Real_Solve;

  procedure Standard_Complex_Solve
              ( dim,want_dcp : in integer32;
                mat : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Makes a righthandside vector b so x = (1,1,..,1) is the exact solution.
  --   Solves the linear system in standard double floating-point arithmetic.
  --   Verifies whether the computed solution has indeed at least as many 
  --   decimal places correct as the value of want_dcp.

    use Standard_Complex_Numbers;

    one : constant Complex_Number := create(1.0);
    x : constant Standard_Complex_Vectors.Vector(1..dim) := (1..dim => one);
    b : Standard_Complex_Vectors.Vector(1..dim);
    wrk : Standard_Complex_Matrices.Matrix(1..dim,1..dim) := mat;
    piv : Standard_Integer_Vectors.Vector(1..dim);
    inf : integer32;
    tol : constant double_float := 10.0**(integer(-want_dcp - 1));
    okay : boolean := true;
    adf : double_float;
  
    use Standard_Complex_Matrices;

  begin
    b := mat*x;
    Standard_Complex_Linear_Solvers.lufac(wrk,dim,piv,inf);
    Standard_Complex_Linear_Solvers.lusolve(wrk,dim,piv,b);
    put_line("The computed solution :"); put_line(b);
    for i in b'range loop
      adf := AbsVal(b(i) - x(i));
      if adf > tol then
        put("x("); put(i,1); put(") = "); put(b(i));
        put_line(" is off!"); okay := false;
      end if;
      exit when not okay;
    end loop;
    if okay then
      put("Solution accurate with at least ");
      put(want_dcp,1); put_line(" decimal places.");
    else
      put("Accuracy of "); put(want_dcp,1); 
      put_line(" decimal places is not obtained.");
    end if;
  end Standard_Complex_Solve;

  procedure DoblDobl_Real_Solve
              ( dim,want_dcp : in integer32;
                mat : in Double_Double_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Makes a righthandside vector b so x = (1,1,..,1) is the exact solution.
  --   Solves the linear system in double double floating-point arithmetic.
  --   Verifies whether the computed solution has indeed at least as many 
  --   decimal places correct as the value of want_dcp.

    one : constant double_double := create(1.0);
    x : constant Double_Double_Vectors.Vector(1..dim) := (1..dim => one);
    b : Double_Double_Vectors.Vector(1..dim);
    wrk : Double_Double_Matrices.Matrix(1..dim,1..dim) := mat;
    piv : Standard_Integer_Vectors.Vector(1..dim);
    inf : integer32;
    tol : constant double_float := 10.0**(integer(-want_dcp - 1));
    okay : boolean := true;
    adf : double_double;
  
    use Double_Double_Matrices;

  begin
    b := mat*x;
    Double_Double_Linear_Solvers.lufac(wrk,dim,piv,inf);
    Double_Double_Linear_Solvers.lusolve(wrk,dim,piv,b);
    put_line("The computed solution :"); put_line(b);
    for i in b'range loop
      adf := abs(b(i) - x(i));
      if adf > tol then
        put("x("); put(i,1); put(") = "); put(b(i));
        put_line(" is off!"); okay := false;
      end if;
      exit when not okay;
    end loop;
    if okay then
      put("Solution accurate with at least ");
      put(want_dcp,1); put_line(" decimal places.");
    else
      put("Accuracy of "); put(want_dcp,1); 
      put_line(" decimal places is not obtained.");
    end if;
  end DoblDobl_Real_Solve;

  procedure DoblDobl_Complex_Solve
              ( dim,want_dcp : in integer32;
                mat : in DoblDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Makes a righthandside vector b so x = (1,1,..,1) is the exact solution.
  --   Solves the linear system in double double floating-point arithmetic.
  --   Verifies whether the computed solution has indeed at least as many 
  --   decimal places correct as the value of want_dcp.

    use DoblDobl_Complex_Numbers;

    dd_one : constant double_double := create(1.0);
    one : constant Complex_Number := create(dd_one);
    x : constant DoblDobl_Complex_Vectors.Vector(1..dim) := (1..dim => one);
    b : DoblDobl_Complex_Vectors.Vector(1..dim);
    wrk : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim) := mat;
    piv : Standard_Integer_Vectors.Vector(1..dim);
    inf : integer32;
    tol : constant double_float := 10.0**(integer(-want_dcp - 1));
    okay : boolean := true;
    adf : double_double;
  
    use DoblDobl_Complex_Matrices;

  begin
    b := mat*x;
    DoblDobl_Complex_Linear_Solvers.lufac(wrk,dim,piv,inf);
    DoblDobl_Complex_Linear_Solvers.lusolve(wrk,dim,piv,b);
    put_line("The computed solution :"); put_line(b);
    for i in b'range loop
      adf := AbsVal(b(i) - x(i));
      if adf > tol then
        put("x("); put(i,1); put(") = "); put(b(i));
        put_line(" is off!"); okay := false;
      end if;
      exit when not okay;
    end loop;
    if okay then
      put("Solution accurate with at least ");
      put(want_dcp,1); put_line(" decimal places.");
    else
      put("Accuracy of "); put(want_dcp,1); 
      put_line(" decimal places is not obtained.");
    end if;
  end DoblDobl_Complex_Solve;

  procedure Real_Test ( n : in integer32; c : in double_float ) is

  -- DESCRIPTION :
  --   Generates a random n-dimensional matrix with condition number c
  --   and calls the LU factorization with condition number estimator.

    mat : Standard_Floating_Matrices.Matrix(1..n,1..n)
        := Random_Conditioned_Matrix(n,c);
    loss_dcp,want_dcp : integer32 := 0;
    precision : integer32 := 15;

  begin
    new_line;
    put("The given condition number :"); put(c,3); new_line;
    loss_dcp := Estimated_Loss(mat);
    put("Estimated loss of decimal places : "); put(loss_dcp,1); new_line;
    new_line;
    put("Give the wanted number of decimal places : "); get(want_dcp);
    precision := precision + loss_dcp;
    put("Number of decimal places left in precision : ");
    put(precision,1); new_line;
    if precision >= want_dcp then
      put_line("Standard double precision will suffice, solving ...");
      Standard_Real_Solve(n,want_dcp,mat);
    else
      put_line("Standard double precision will not suffice.");
      precision := 31 + loss_dcp;
      if precision >= want_dcp then
        declare
          dd_mat : constant Double_Double_Matrices.Matrix(1..n,1..n)
                 := d2dd(mat);
        begin
          put_line("Double double precision will suffice.");
          DoblDobl_Real_Solve(n,want_dcp,dd_mat);
        end;
      else
        put_line("Double double precision will not suffice.");
      end if;
    end if;
  end Real_Test;

  procedure Complex_Test ( n : in integer32; c : in double_float ) is

  -- DESCRIPTION :
  --   Generates a random n-dimensional matrix with condition number c
  --   and calls the LU factorization with condition number estimator.

    mat : Standard_Complex_Matrices.Matrix(1..n,1..n)
        := Random_Conditioned_Matrix(n,c);
    loss_dcp,want_dcp : integer32 := 0;
    precision : integer32 := 15;

  begin
    new_line;
    put("The given condition number :"); put(c,3); new_line;
    loss_dcp := Estimated_Loss(mat);
    put("Estimated loss of decimal places : "); put(loss_dcp,1); new_line;
    new_line;
    put("Give the wanted number of decimal places : "); get(want_dcp);
    precision := precision + loss_dcp;
    put("Number of decimal places left in precision : ");
    put(precision,1); new_line;
    if precision >= want_dcp then
      put_line("Standard double precision will suffice, solving ...");
      Standard_Complex_Solve(n,want_dcp,mat);
    else
      put_line("Standard double precision will not suffice.");
      precision := 31 + loss_dcp;
      if precision >= want_dcp then
        declare
          dd_mat : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
                 := d2dd(mat);
        begin
          put_line("Double double precision will suffice.");
          DoblDobl_Complex_Solve(n,want_dcp,dd_mat);
        end;
      else
        put_line("Double double precision will not suffice.");
      end if;
    end if;
  end Complex_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension, condition number,
  --   and ask whether real or complex arithmetic is needed.

    dim : integer32 := 0;
    cnd : double_float := 1.0;
    ans : character;

  begin
    new_line;
    put_line("Interactive test on variable precision LU decomposition ...");
    put("Give the dimension : "); get(dim);
    put("Give the condition number : "); get(cnd);
    put("Real or complex arithmetic ? (r/c) ");
    Ask_Alternative(ans,"rc");
    if ans = 'r'
     then Real_Test(dim,cnd);
     else Complex_Test(dim,cnd);
    end if;
  end Main;

begin
  Main;
end ts_vmplu;
