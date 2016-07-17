with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;         use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Standard_Complex_Linear_Solvers;     use Standard_Complex_Linear_Solvers;
with Standard_Random_Series;
with Standard_Dense_Series_Vectors;
with Standard_Dense_Series_Vectors_io;    use Standard_Dense_Series_Vectors_io;
with Standard_Dense_Series_Matrices;
with Standard_Dense_Vector_Series;
with Standard_Dense_Vector_Series_io;     use Standard_Dense_Vector_Series_io;
with Standard_Dense_Matrix_Series;
with Standard_Dense_Matrix_Series_io;     use Standard_Dense_Matrix_Series_io;
with Standard_Matrix_Series_Solvers;      use Standard_Matrix_Series_Solvers;

procedure ts_serlin is

-- DESCRIPTION :
--   Tests the linearization of solving linear systems of truncated series.

  procedure Solve ( A : in Standard_Dense_Matrix_Series.Matrix;
                    b : in Standard_Dense_Vector_Series.Vector;
                    x : out Standard_Dense_Vector_Series.Vector ) is

  -- DESCRIPTION :
  --   Solves the linear system A*x = b.

  -- REQUIRED : A.deg >= 0 and b.deg >= 0.

    dim : constant integer32 := A.cff(0)'last;
    lwrk : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    info : integer32;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);

  begin
    Solve_Lead_by_lufac(A,b,lwrk,ipvt,info,x);
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      put_line("The leading vector series of the solution :");
      put_line(x.cff(0));
      Solve_Next_by_lusolve(A,b,lwrk,ipvt,x);
    end if;
  end Solve;

  procedure Standard_Test ( n,m,d : in integer32 ) is

  -- DESCRIPTION :
  --   Converts an n-by-m matrix of series of degree d with standard
  --   double precision complex coefficients into a matrix series.

    use Standard_Dense_Series_Matrices;

    sA : constant Standard_Dense_Series_Matrices.Matrix(1..n,1..m)
       := Standard_Random_Series.Random_Series_Matrix(1,n,1,m,d);
    As : constant Standard_Dense_Matrix_Series.Matrix 
       := Standard_Dense_Matrix_Series.Create(sA); 
    sx : constant Standard_Dense_Series_Vectors.Vector(1..m)
       := Standard_Random_Series.Random_Series_Vector(1,m,d);
    xs : Standard_Dense_Vector_Series.Vector
       := Standard_Dense_Vector_Series.Create(sx);
    sb : constant Standard_Dense_Series_Vectors.Vector(1..n) := sA*sx;
    bs : Standard_Dense_Vector_Series.Vector
       := Standard_Dense_Vector_Series.Create(sb);
    ys : Standard_Dense_Vector_Series.Vector;

  begin
    put_line("The coefficients of the matrix series :"); put(As);
    put_line("The exact solution x :"); put_line(sx);
    put_line("The coefficients of the vector series x :"); put(xs);
    put_line("The right hand side vector b :"); put_line(sb);
    put_line("The coefficients of the vector series b :"); put(bs);
    Solve(As,bs,ys);
    put_line("The generated leading vector series of the solution :");
    put_line(xs.cff(0));
    put_line("The computed leading vector series of the solution :");
    put_line(ys.cff(0));
    put_line("The generated 2nd term of the vector series of the solution :");
    put_line(xs.cff(1));
    put_line("The computed 2nd term of the vector series of the solution :");
    put_line(ys.cff(1));
  end Standard_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension of the linear system
  --   and the degrees of the series in the system.

    dim,deg : integer32 := 0;

  begin
    new_line;
    put_line("Testing the linearization of systems of power series ...");
    put("  Give the dimension of the system : "); get(dim);
    put("  Give the degree of the series : "); get(deg);
    Standard_Test(dim,dim,deg);
  end Main;

begin
  Main;
end ts_serlin;
