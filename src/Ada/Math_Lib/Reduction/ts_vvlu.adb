with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Multprec_Complex_Matrices;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with DoblDobl_Random_Matrices;           use DoblDobl_Random_Matrices;
with QuadDobl_Random_Matrices;           use QuadDobl_Random_Matrices;
with Multprec_Random_Matrices;           use Multprec_Random_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with Multprec_Complex_Linear_Solvers;    use Multprec_Complex_Linear_Solvers;

procedure ts_vvlu is

-- DESCRIPTION :
--   Test on LU factorization on matrices given as vectors of columns.

  function mat2vv ( A : Standard_Complex_Matrices.Matrix )
                  return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns a vector of vectors, with as content in the vectors
  --   the columns of the matrix A.

    res : Standard_Complex_VecVecs.VecVec(A'range(2));
    col : Standard_Complex_Vectors.Vector(A'range(1));

  begin
    for k in A'range(2) loop
      for i in col'range loop
        col(i) := A(i,k); 
      end loop;
      res(k) := new Standard_Complex_Vectors.Vector'(col);
    end loop;
    return res;
  end mat2vv;

  function mat2vv ( A : DoblDobl_Complex_Matrices.Matrix )
                  return DoblDobl_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns a vector of vectors, with as content in the vectors
  --   the columns of the matrix A.

    res : DoblDobl_Complex_VecVecs.VecVec(A'range(2));
    col : DoblDobl_Complex_Vectors.Vector(A'range(1));

  begin
    for k in A'range(2) loop
      for i in col'range loop
        col(i) := A(i,k); 
      end loop;
      res(k) := new DoblDobl_Complex_Vectors.Vector'(col);
    end loop;
    return res;
  end mat2vv;

  function mat2vv ( A : QuadDobl_Complex_Matrices.Matrix )
                  return QuadDobl_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns a vector of vectors, with as content in the vectors
  --   the columns of the matrix A.

    res : QuadDobl_Complex_VecVecs.VecVec(A'range(2));
    col : QuadDobl_Complex_Vectors.Vector(A'range(1));

  begin
    for k in A'range(2) loop
      for i in col'range loop
        col(i) := A(i,k); 
      end loop;
      res(k) := new QuadDobl_Complex_Vectors.Vector'(col);
    end loop;
    return res;
  end mat2vv;

  function mat2vv ( A : Multprec_Complex_Matrices.Matrix )
                  return Multprec_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns a vector of vectors, with as content in the vectors
  --   the columns of the matrix A.

    res : Multprec_Complex_VecVecs.VecVec(A'range(2));

  begin
    for k in A'range(2) loop
      declare
        col : Multprec_Complex_Vectors.Vector(A'range(1));
      begin
        for i in col'range loop
          Multprec_Complex_Numbers.Copy(A(i,k),col(i)); 
        end loop;
        res(k) := new Multprec_Complex_Vectors.Vector'(col);
      end;
    end loop;
    return res;
  end mat2vv;

  procedure Compare ( A : in Standard_Complex_Matrices.Matrix;
                      B : in Standard_Complex_VecVecs.VecVec;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the entries in A with the corresponding entries in B.
  --   Writes the results to screen if output is true or if the 
  --   difference between the two entries exceeds the tolerance tol.

    use Standard_Complex_Numbers;
    dff : Complex_Number;
    val : double_float;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        dff := A(i,j) - B(j)(i);
        val := AbsVal(dff);
        if output or val > tol then
          put(" A("); put(i,1); put(","); put(j,1); put(") : ");
          put(A(i,j)); new_line;
          put("B("); put(j,1); put(")("); put(i,1); put(") : ");
          put(B(j)(i)); new_line;
          put("error : "); put(val,3); new_line;
        end if;
      end loop;
    end loop;
  end Compare;

  procedure Compare ( A : in DoblDobl_Complex_Matrices.Matrix;
                      B : in DoblDobl_Complex_VecVecs.VecVec;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the entries in A with the corresponding entries in B.
  --   Writes the results to screen if output is true or if the 
  --   difference between the two entries exceeds the tolerance tol.

    use DoblDobl_Complex_Numbers;
    dff : Complex_Number;
    val : double_double;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        dff := A(i,j) - B(j)(i);
        val := AbsVal(dff);
        if output or val > tol then
          put(" A("); put(i,1); put(","); put(j,1); put(") : ");
          put(A(i,j)); new_line;
          put("B("); put(j,1); put(")("); put(i,1); put(") : ");
          put(B(j)(i)); new_line;
          put("error : "); put(val,3); new_line;
        end if;
      end loop;
    end loop;
  end Compare;

  procedure Compare ( A : in QuadDobl_Complex_Matrices.Matrix;
                      B : in QuadDobl_Complex_VecVecs.VecVec;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the entries in A with the corresponding entries in B.
  --   Writes the results to screen if output is true or if the 
  --   difference between the two entries exceeds the tolerance tol.

    use QuadDobl_Complex_Numbers;
    dff : Complex_Number;
    val : quad_double;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        dff := A(i,j) - B(j)(i);
        val := AbsVal(dff);
        if output or val > tol then
          put(" A("); put(i,1); put(","); put(j,1); put(") : ");
          put(A(i,j)); new_line;
          put("B("); put(j,1); put(")("); put(i,1); put(") : ");
          put(B(j)(i)); new_line;
          put("error : "); put(val,3); new_line;
        end if;
      end loop;
    end loop;
  end Compare;

  procedure Compare ( A : in Multprec_Complex_Matrices.Matrix;
                      B : in Multprec_Complex_VecVecs.VecVec;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the entries in A with the corresponding entries in B.
  --   Writes the results to screen if output is true or if the 
  --   difference between the two entries exceeds the tolerance tol.

    use Multprec_Complex_Numbers;
    dff : Complex_Number;
    val : Floating_Number;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        dff := A(i,j) - B(j)(i);
        val := AbsVal(dff);
        if output or val > tol then
          put(" A("); put(i,1); put(","); put(j,1); put(") : ");
          put(A(i,j)); new_line;
          put("B("); put(j,1); put(")("); put(i,1); put(") : ");
          put(B(j)(i)); new_line;
          put("error : "); put(val,3); new_line;
        end if;
        Clear(dff); Clear(val);
      end loop;
    end loop;
  end Compare;

  procedure Standard_Test ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Performs a test in standard double precision,
  --   on a problem of the given dimension dim.

    use Standard_Complex_VecVecs;
    use Standard_Complex_Matrices;

    A : Matrix(1..dim,1..dim) := Random_Matrix(1,dim,1,dim);
    B : VecVec(1..dim) := mat2vv(A);
    Apiv,Bpiv : Standard_Integer_Vectors.Vector(1..dim);
    Ainfo,Binfo : integer32;
    otp : boolean;
    ans : character;

  begin
    new_line;
    put("Intermediate output wanted ? (y/n) ");
    Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    new_line;
    put_line("Checking the matrices on input ...");
    Compare(A,B,1.0E-8,otp);
    lufac(A,dim,Apiv,Ainfo);
    put("pivots on A : "); put(Apiv); new_line;
    lufac(B,dim,Bpiv,Binfo);
    put("pivots on B : "); put(Bpiv); new_line;
    put_line("Checking the matrices on output ...");
    Compare(A,B,1.0E-8,otp);
  end Standard_Test;

  procedure DoblDobl_Test ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Performs a test in double double precision,
  --   on a problem of the given dimension dim.

    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Complex_Matrices;
   
    n : constant natural32 := natural32(dim);
    A : Matrix(1..dim,1..dim) := Random_Matrix(n,n);
    B : VecVec(1..dim) := mat2vv(A);
    Apiv,Bpiv : Standard_Integer_Vectors.Vector(1..dim);
    Ainfo,Binfo : integer32;
    otp : boolean;
    ans : character;

  begin
    new_line;
    put("Intermediate output wanted ? (y/n) ");
    Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    new_line;
    put_line("Checking the matrices on input ...");
    Compare(A,B,1.0E-8,otp);
    lufac(A,dim,Apiv,Ainfo);
    put("pivots on A : "); put(Apiv); new_line;
    lufac(B,dim,Bpiv,Binfo);
    put("pivots on B : "); put(Bpiv); new_line;
    put_line("Checking the matrices on output ...");
    Compare(A,B,1.0E-8,otp);
  end DoblDobl_Test;

  procedure QuadDobl_Test ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Performs a test in double double precision,
  --   on a problem of the given dimension dim.

    use QuadDobl_Complex_VecVecs;
    use QuadDobl_Complex_Matrices;
   
    n : constant natural32 := natural32(dim);
    A : Matrix(1..dim,1..dim) := Random_Matrix(n,n);
    B : VecVec(1..dim) := mat2vv(A);
    Apiv,Bpiv : Standard_Integer_Vectors.Vector(1..dim);
    Ainfo,Binfo : integer32;
    otp : boolean;
    ans : character;

  begin
    new_line;
    put("Intermediate output wanted ? (y/n) ");
    Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    new_line;
    put_line("Checking the matrices on input ...");
    Compare(A,B,1.0E-8,otp);
    lufac(A,dim,Apiv,Ainfo);
    put("pivots on A : "); put(Apiv); new_line;
    lufac(B,dim,Bpiv,Binfo);
    put("pivots on B : "); put(Bpiv); new_line;
    put_line("Checking the matrices on output ...");
    Compare(A,B,1.0E-8,otp);
  end QuadDobl_Test;

  procedure Multprec_Test ( dim : in integer32; size : in natural32 ) is

  -- DESCRIPTION :
  --   Performs a test in double double precision,
  --   on a problem of the given dimension dim.
  --   The size of the numbers is defined by size.

    use Multprec_Complex_VecVecs;
    use Multprec_Complex_Matrices;
   
    n : constant natural32 := natural32(dim);
    A : Matrix(1..dim,1..dim) := Random_Matrix(n,n,size);
    B : VecVec(1..dim) := mat2vv(A);
    Apiv,Bpiv : Standard_Integer_Vectors.Vector(1..dim);
    Ainfo,Binfo : integer32;
    otp : boolean;
    ans : character;

  begin
    new_line;
    put("Intermediate output wanted ? (y/n) ");
    Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    new_line;
    put_line("Checking the matrices on input ...");
    Compare(A,B,1.0E-8,otp);
    lufac(A,dim,Apiv,Ainfo);
    put("pivots on A : "); put(Apiv); new_line;
    lufac(B,dim,Bpiv,Binfo);
    put("pivots on B : "); put(Bpiv); new_line;
    put_line("Checking the matrices on output ...");
    Compare(A,B,1.0E-8,otp);
  end Multprec_Test;

  procedure Main is

    dim : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("LU factorization of matrices as vectors of columns ...");
    put("Give the dimension : "); get(dim);
    new_line;
    put_line("MENU to select the precision : ");
    put_line("  0. standard double precision; or");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision; or");
    put_line("  3. arbitrary multiprecision.");
    put("Type 0, 1, 2, or 3 to select the precision : ");
    Ask_Alternative(ans,"0123");
    case ans is
      when '0' => Standard_Test(dim);
      when '1' => DoblDobl_Test(dim);
      when '2' => QuadDobl_Test(dim);
      when '3' =>
        declare
          deci,size : natural32 := 0;
        begin
          new_line;
          put("Give the number of decimal places : "); get(deci);
          size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
          Multprec_Test(dim,size);
        end;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_vvlu;
