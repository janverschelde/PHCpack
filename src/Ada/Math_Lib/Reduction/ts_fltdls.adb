with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Floating_Linear_Solvers;   use Standard_Floating_LInear_Solvers;
with Standard_Random_Matrices;           use Standard_Random_Matrices;

procedure ts_fltdls is

-- DESCRIPTION :
--   Test the dynamic triangulators of floating-point matrices.

  tol : constant double_float := 10.0**(-8);

  procedure Write_Matrix ( mat : in Matrix ) is
  begin
    for i in mat'range(1) loop
      for j in mat'range(2) loop
        put(mat(i,j),3,3,3);
      end loop;
      new_line;
    end loop;
  end Write_Matrix;

  procedure Write_Vector ( v : in Standard_Floating_Vectors.Vector ) is
  begin
    for i in v'range loop
      put(v(i),3,3,3);
    end loop;
    new_line;
  end Write_Vector;

  function Product ( n,col : integer32; mat : Matrix;
                     sol : Standard_Floating_Vectors.Vector )
                   return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the product of the solution with the columns of the matrix.

    res : Standard_Floating_Vectors.Vector(1..n);

  begin
    for i in 1..n loop
      res(i) := mat(i,col)*sol(n+1);
      for j in 1..n loop
        res(i) := res(i) + mat(i,j)*sol(j);
      end loop;
    end loop;
    return res;
  end Product;

  function Residual ( n,col : integer32; mat : Matrix;
                      ipvt : Standard_Integer_Vectors.Vector;
                      sol : Standard_Floating_Vectors.Vector )
                    return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the residual of the solution.

    res : Standard_Floating_Vectors.Vector(1..n);

  begin
    for i in 1..n loop
      res(i) := mat(i,ipvt(col))*sol(n+1);
      for j in 1..n loop
        res(i) := res(i) + mat(i,ipvt(j))*sol(j);
      end loop;
    end loop;
    return res;
  end Residual;

  procedure Upper_Triangulate ( n,m : in integer32; mat : in Matrix ) is

  -- DESCRIPTION :
  --   Triangulates the matrix, computes a solution and checks residual.

    wrk : Matrix(1..n,1..m) := mat;
    ipvt : Standard_Integer_Vectors.Vector(1..m);
    pivot : integer32;
    sol : Standard_Floating_Vectors.Vector(1..n+1);
    res,prod : Standard_Floating_Vectors.Vector(1..n);

  begin
    put_line("The matrix on entry : "); Write_Matrix(mat);
    for i in 1..m loop
      ipvt(i) := i;
    end loop;
    for i in 1..n loop
      Upper_Triangulate(i,wrk,tol,ipvt,pivot);
      exit when (pivot = 0);
    end loop;
    put_line("The triangulated matrix : "); Write_Matrix(wrk);
    put("Pivots : "); put(ipvt); new_line;
    for i in n+1..m loop
      sol := Solve(n,i,wrk);
      put("Solution : "); Write_Vector(sol);
      prod := Product(n,i,wrk,sol);
      put("Product : "); Write_Vector(prod);
      res := Residual(n,i,mat,ipvt,sol);
      put("Residual : "); Write_Vector(res);
    end loop;
  end Upper_Triangulate;

  procedure Interactive_Testing is

    n,m : integer32 := 0;

  begin
    put("Give number of rows : "); get(n);
    put("Give number of columns : "); get(m);
    put("Give "); put(n,1); put("x"); put(m,1); put_line(" matrix : ");
    declare
      mat : Matrix(1..n,1..m);
    begin
      get(mat);
      Upper_Triangulate(n,m,mat);
    end;
  end Interactive_Testing;

  procedure Random_Testing is

    n,m,nb : integer32 := 0;

  begin 
    put("Give number of rows : "); get(n);
    put("Give number of columns : "); get(m);
    put("Give number of tests : "); get(nb);
    for i in 1..nb loop
      declare
        mat : constant Matrix(1..n,1..m)
            := Random_Matrix(natural32(n),natural32(m));
      begin
        Upper_Triangulate(n,m,mat);
      end;
    end loop;
  end Random_Testing;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing the dynamic triangulators of floating-point matrices.");
    loop
      new_line;
      put_line("Choose one of the following :                           ");
      put_line("  0. Exit this program.                                 ");
      put_line("  1. Interactive testing of dynamic triangulators.      ");
      put_line("  2. Random testing of dynamic triangulators.           ");
      put("Type 0, 1, or 2 to select : "); get(ans);
      exit when (ans = '0');
      case ans is
        when '1' => Interactive_Testing;
        when '2' => Random_Testing;
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_fltdls;
