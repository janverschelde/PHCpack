with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;

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

  procedure Standard_Test ( dim : integer32 ) is

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
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_vvlu;
