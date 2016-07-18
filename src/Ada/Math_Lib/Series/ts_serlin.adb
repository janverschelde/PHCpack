with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;         use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;         use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;         use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Matrices;
with Standard_Random_Series;
with Standard_Dense_Series_Vectors;
with Standard_Dense_Series_Vectors_io;    use Standard_Dense_Series_Vectors_io;
with Standard_Dense_Series_Matrices;
with Standard_Dense_Vector_Series;
with Standard_Dense_Vector_Series_io;     use Standard_Dense_Vector_Series_io;
with Standard_Dense_Matrix_Series;
with Standard_Dense_Matrix_Series_io;     use Standard_Dense_Matrix_Series_io;
with Standard_Matrix_Series_Solvers;
with DoblDobl_Random_Series;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Dense_Series_Vectors_io;    use DoblDobl_Dense_Series_Vectors_io;
with DoblDobl_Dense_Series_Matrices;
with DoblDobl_Dense_vector_Series;
with DoblDobl_Dense_vector_Series_io;     use DoblDobl_Dense_vector_Series_io;
with DoblDobl_Dense_Matrix_Series;
with DoblDobl_Dense_Matrix_Series_io;     use DoblDobl_Dense_Matrix_Series_io;
with DoblDobl_Matrix_Series_Solvers;
with QuadDobl_Random_Series;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series_Vectors_io;    use QuadDobl_Dense_Series_Vectors_io;
with QuadDobl_Dense_Series_Matrices;
with QuadDobl_Dense_vector_Series;
with QuadDobl_Dense_vector_Series_io;     use QuadDobl_Dense_vector_Series_io;
with QuadDobl_Dense_Matrix_Series;
with QuadDobl_Dense_Matrix_Series_io;     use QuadDobl_Dense_Matrix_Series_io;
with QuadDobl_Matrix_Series_Solvers;

procedure ts_serlin is

-- DESCRIPTION :
--   Tests the linearization of solving linear systems of truncated series.

  procedure Solve ( A : in Standard_Dense_Matrix_Series.Matrix;
                    b : in Standard_Dense_Vector_Series.Vector;
                    x : out Standard_Dense_Vector_Series.Vector ) is

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, in standard double precision.

  -- REQUIRED : A.deg >= 0 and b.deg >= 0.

    dim : constant integer32 := A.cff(0)'last;
    lwrk : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    info : integer32;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);

    use Standard_Matrix_Series_Solvers;
   
  begin
    Solve_Lead_by_lufac(A,b,lwrk,ipvt,info,x);
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      put_line("The leading vector series of the solution :");
      put_line(x.cff(0));
      for k in 1..b.deg loop
        Solve_Next_by_lusolve(A,b,lwrk,ipvt,x);
      end loop;
    end if;
  end Solve;

  procedure Solve ( A : in DoblDobl_Dense_Matrix_Series.Matrix;
                    b : in DoblDobl_Dense_Vector_Series.Vector;
                    x : out DoblDobl_Dense_Vector_Series.Vector ) is

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, in double double precision.

  -- REQUIRED : A.deg >= 0 and b.deg >= 0.

    dim : constant integer32 := A.cff(0)'last;
    lwrk : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    info : integer32;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);

    use DoblDobl_Matrix_Series_Solvers;

  begin
    Solve_Lead_by_lufac(A,b,lwrk,ipvt,info,x);
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      put_line("The leading vector series of the solution :");
      put_line(x.cff(0));
      for k in 1..b.deg loop
        Solve_Next_by_lusolve(A,b,lwrk,ipvt,x);
      end loop;
    end if;
  end Solve;

  procedure Solve ( A : in QuadDobl_Dense_Matrix_Series.Matrix;
                    b : in QuadDobl_Dense_Vector_Series.Vector;
                    x : out QuadDobl_Dense_Vector_Series.Vector ) is

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, in double double precision.

  -- REQUIRED : A.deg >= 0 and b.deg >= 0.

    dim : constant integer32 := A.cff(0)'last;
    lwrk : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    info : integer32;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);

    use QuadDobl_Matrix_Series_Solvers;

  begin
    Solve_Lead_by_lufac(A,b,lwrk,ipvt,info,x);
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      put_line("The leading vector series of the solution :");
      put_line(x.cff(0));
      for k in 1..b.deg loop
        Solve_Next_by_lusolve(A,b,lwrk,ipvt,x);
      end loop;
    end if;
  end Solve;

  procedure Standard_Test ( n,m,d : in integer32 ) is

  -- DESCRIPTION :
  --   Generates an n-by-m matrix of series of degree d,
  --   with complex coefficients in standard double precision.
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
    for k in 1..bs.deg loop
      put("The generated term "); put(k,1);
      put_line(" of the vector series of the solution :");
      put_line(xs.cff(k));
      put("The computed term "); put(k,1);
      put_line(" of the vector series of the solution :");
      put_line(ys.cff(k));
    end loop;
  end Standard_Test;

  procedure DoblDobl_Test ( n,m,d : in integer32 ) is

  -- DESCRIPTION :
  --   Generates and n-by-m matrix of series of degree d,
  --   with complex coefficients in double double precision.
  --   Converts an n-by-m matrix of series of degree d with double
  --   double precision complex coefficients into a matrix series.

    use DoblDobl_Dense_Series_Matrices;

    sA : constant DoblDobl_Dense_Series_Matrices.Matrix(1..n,1..m)
       := DoblDobl_Random_Series.Random_Series_Matrix(1,n,1,m,d);
    As : constant DoblDobl_Dense_Matrix_Series.Matrix 
       := DoblDobl_Dense_Matrix_Series.Create(sA); 
    sx : constant DoblDobl_Dense_Series_Vectors.Vector(1..m)
       := DoblDobl_Random_Series.Random_Series_Vector(1,m,d);
    xs : DoblDobl_Dense_Vector_Series.Vector
       := DoblDobl_Dense_Vector_Series.Create(sx);
    sb : constant DoblDobl_Dense_Series_Vectors.Vector(1..n) := sA*sx;
    bs : DoblDobl_Dense_Vector_Series.Vector
       := DoblDobl_Dense_Vector_Series.Create(sb);
    ys : DoblDobl_Dense_Vector_Series.Vector;

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
    for k in 1..bs.deg loop
      put("The generated term "); put(k,1);
      put_line(" of the vector series of the solution :");
      put_line(xs.cff(k));
      put("The computed term "); put(k,1);
      put_line(" of the vector series of the solution :");
      put_line(ys.cff(k));
    end loop;
  end DoblDobl_Test;

  procedure QuadDobl_Test ( n,m,d : in integer32 ) is

  -- DESCRIPTION :
  --   Generates and n-by-m matrix of series of degree d,
  --   with complex coefficients in quad double precision.
  --   Converts an n-by-m matrix of series of degree d with quad
  --   double precision complex coefficients into a matrix series.

    use QuadDobl_Dense_Series_Matrices;

    sA : constant QuadDobl_Dense_Series_Matrices.Matrix(1..n,1..m)
       := QuadDobl_Random_Series.Random_Series_Matrix(1,n,1,m,d);
    As : constant QuadDobl_Dense_Matrix_Series.Matrix 
       := QuadDobl_Dense_Matrix_Series.Create(sA); 
    sx : constant QuadDobl_Dense_Series_Vectors.Vector(1..m)
       := QuadDobl_Random_Series.Random_Series_Vector(1,m,d);
    xs : QuadDobl_Dense_Vector_Series.Vector
       := QuadDobl_Dense_Vector_Series.Create(sx);
    sb : constant QuadDobl_Dense_Series_Vectors.Vector(1..n) := sA*sx;
    bs : QuadDobl_Dense_Vector_Series.Vector
       := QuadDobl_Dense_Vector_Series.Create(sb);
    ys : QuadDobl_Dense_Vector_Series.Vector;

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
    for k in 1..bs.deg loop
      put("The generated term "); put(k,1);
      put_line(" of the vector series of the solution :");
      put_line(xs.cff(k));
      put("The computed term "); put(k,1);
      put_line(" of the vector series of the solution :");
      put_line(ys.cff(k));
    end loop;
  end QuadDobl_Test;

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

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension of the linear system
  --   and the degrees of the series in the system.

    dim,deg : integer32 := 0;
    prc : character;

  begin
    new_line;
    put_line("Testing the linearization of systems of power series ...");
    put("  Give the dimension of the system : "); get(dim);
    put("  Give the degree of the series : "); get(deg);
    prc := Prompt_for_Precision;
    case prc is
      when '0' => Standard_Test(dim,dim,deg);
      when '1' => DoblDobl_Test(dim,dim,deg);
      when '2' => QuadDobl_Test(dim,dim,deg);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_serlin;
