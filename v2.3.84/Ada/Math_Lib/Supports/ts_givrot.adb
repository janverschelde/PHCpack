with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Givens_Rotations;                   use Givens_Rotations;

procedure ts_givrot is

-- DESCRIPTION :
--   Interactive testing of Givens rotations.

  n,m : integer32 := 0;
  ans : character;
  tol : constant double_float := 10.0**(-12);

  function Compute_Residual ( mat : in matrix; rhs,x : in vector;
                              ipvt : in Standard_Integer_Vectors.Vector )
                            return Vector is

    res : Vector(rhs'range) := rhs;

  begin
    for i in mat'range(1) loop
      for j in mat'range(2) loop
        res(i) := res(i) - mat(i,ipvt(j))*x(j);
        exit when j > mat'last(2);
      end loop;
    end loop;
    return res;
  end Compute_Residual;

  procedure Row_Wise_Triangulation ( mat : in out Matrix ) is

    pivots : Standard_Integer_Vectors.Vector(mat'range(2));
    ipiv : integer32;
 
  begin
    for i in pivots'range loop
      pivots(i) := i;
    end loop;
    for i in mat'range(1) loop
      Upper_Triangulate(i,mat,tol,pivots,ipiv);
      put("At row "); put(i,1); put(" we have pivot : "); put(ipiv,1);
      new_line;
      put_line("with matrix : "); put(mat);
      put("and pivot vector "); put(pivots); new_line;
      exit when (ipiv = 0);
    end loop;
  end Row_Wise_Triangulation;

begin
  loop
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    declare
      mat,wrkmat : matrix(1..n,1..m);
      rhs,wrkrhs : vector(1..n);
      ipvt : Standard_Integer_Vectors.Vector(1..m);
    begin
      put("Give the elements of the ");
      put(n,1); put("x"); put(m,1); put_line("-matrix :"); get(mat);
      put_line("The matrix :"); put(mat);
      wrkmat := mat;
      put("Do you have a right hand side ? (y/n) "); get(ans);
      if ans = 'y' then
        put("Give "); put(n,1); put(" floating point numbers : ");
        rhs := (1..n => 0.0); get(rhs);
        wrkrhs := rhs;
        Upper_Triangulate(wrkmat,wrkrhs,tol,ipvt);
      else
        Row_Wise_Triangulation(wrkmat);
        wrkmat := mat;
        Upper_Triangulate(wrkmat,tol,ipvt);
      end if;
      put_line("The matrix in upper triangular form : "); put(wrkmat);
      put(" with pivoting information : "); put(ipvt); new_line;
      if ans = 'y' then
        put(" and transformed right hand side : ");
        put(wrkrhs); new_line;
        declare
          sol,res : Vector(1..n);
        begin
          Solve(wrkmat,wrkrhs,tol,sol);
          put("The solution : "); put(sol); new_line;
          res := Compute_Residual(mat,rhs,sol,ipvt);
          put(" with residual : "); put(res); new_line;
        end;
      end if;
    end;
    put("Do you want more tests ? (y/n) "); get(ans);
    exit when ans /= 'y';
  end loop;
end ts_givrot;
