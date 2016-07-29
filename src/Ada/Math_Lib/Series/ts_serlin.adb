with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Double_Double_Numbers;               use Double_Double_Numbers;
with Double_Double_Numbers_io;            use Double_Double_Numbers_io;
with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with Quad_Double_Numbers_io;              use Quad_Double_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;         use Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;         use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;         use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vector_Norms;
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

  procedure Write_Difference
              ( x,y : in Standard_Dense_Vector_Series.Vector ) is

  -- DESCRIPTION :
  --   Writes the max norm of the difference of each coefficient vector
  --   between x and y.  At the end, writes the largest max norm, as an
  --   upper bound on the error.

  -- REQUIRED : x.deg = y.deg >= 0.

    dim : constant integer32 := x.cff(0)'last;
    dif : Standard_Complex_Vectors.Vector(1..dim);
    nrm,err : double_float := 0.0;

    use Standard_Complex_Vectors;

  begin
    for k in 0..x.deg loop
      dif := x.cff(k).all - y.cff(k).all;
      nrm := Standard_Complex_Vector_Norms.Max_Norm(dif);
      put("Max norm of error at component "); put(k,1);
      put(" :"); put(nrm,3); new_line;
      if nrm > err
       then err := nrm;
      end if;
    end loop;
    put("Max norm of the error :"); put(err,3); new_line;
  end Write_Difference;

  procedure Write_Difference
              ( x,y : in DoblDobl_Dense_Vector_Series.Vector ) is

  -- DESCRIPTION :
  --   Writes the max norm of the difference of each coefficient vector
  --   between x and y.  At the end, writes the largest max norm, as an
  --   upper bound on the error.

  -- REQUIRED : x.deg = y.deg >= 0.

    dim : constant integer32 := x.cff(0)'last;
    dif : DoblDobl_Complex_Vectors.Vector(1..dim);
    nrm,err : double_double := create(0.0);

    use DoblDobl_Complex_Vectors;

  begin
    for k in 0..x.deg loop
      dif := x.cff(k).all - y.cff(k).all;
      nrm := DoblDobl_Complex_Vector_Norms.Max_Norm(dif);
      put("Max norm of error at component "); put(k,1);
      put(" : "); put(nrm,3); new_line;
      if nrm > err
       then err := nrm;
      end if;
    end loop;
    put("Max norm of the error : "); put(err,3); new_line;
  end Write_Difference;

  procedure Write_Difference
              ( x,y : in QuadDobl_Dense_Vector_Series.Vector ) is

  -- DESCRIPTION :
  --   Writes the max norm of the difference of each coefficient vector
  --   between x and y.  At the end, writes the largest max norm, as an
  --   upper bound on the error.

  -- REQUIRED : x.deg = y.deg >= 0.

    dim : constant integer32 := x.cff(0)'last;
    dif : QuadDobl_Complex_Vectors.Vector(1..dim);
    nrm,err : quad_double := create(0.0);

    use QuadDobl_Complex_Vectors;

  begin
    for k in 0..x.deg loop
      dif := x.cff(k).all - y.cff(k).all;
      nrm := QuadDobl_Complex_Vector_Norms.Max_Norm(dif);
      put("Max norm of error at component "); put(k,1);
      put(" : "); put(nrm,3); new_line;
      if nrm > err
       then err := nrm;
      end if;
    end loop;
    put("Max norm of the error : "); put(err,3); new_line;
  end Write_Difference;

  procedure Standard_Test ( n,m,d : in integer32 ) is

  -- DESCRIPTION :
  --   Generates an n-by-m matrix of series of degree d,
  --   with complex coefficients in standard double precision.
  --   Converts an n-by-m matrix of series of degree d with standard
  --   double precision complex coefficients into a matrix series.

    use Standard_Dense_Series_Matrices;
    use Standard_Matrix_Series_Solvers;

    sA : constant Standard_Dense_Series_Matrices.Matrix(1..n,1..m)
       := Standard_Random_Series.Random_Series_Matrix(1,n,1,m,d);
    As : constant Standard_Dense_Matrix_Series.Matrix 
       := Standard_Dense_Matrix_Series.Create(sA); 
    sx : constant Standard_Dense_Series_Vectors.Vector(1..m)
       := Standard_Random_Series.Random_Series_Vector(1,m,d);
    xs : constant Standard_Dense_Vector_Series.Vector
       := Standard_Dense_Vector_Series.Create(sx);
    sb : constant Standard_Dense_Series_Vectors.Vector(1..n) := sA*sx;
    bs : constant Standard_Dense_Vector_Series.Vector
       := Standard_Dense_Vector_Series.Create(sb);
    ys : Standard_Dense_Vector_Series.Vector;
    ans : character;
    rcond : double_float;
    info : integer32;

  begin
    put_line("The coefficients of the matrix series :"); put(As);
    put_line("The exact solution x :"); put_line(sx);
    put_line("The coefficients of the vector series x :"); put(xs);
    put_line("The right hand side vector b :"); put_line(sb);
    put_line("The coefficients of the vector series b :"); put(bs);
    if n > m then
      new_line;
      put("Solve with SVD ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y' then
        Solve_by_SVD(As,bs,info,rcond,ys);
        put("rcond : "); put(rcond,3); new_line;
      else
        Solve_by_QRLS(As,bs,info,ys);
        put("info : "); put(info,1); new_line;
      end if;
    else
      new_line;
      put("Condition number wanted ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y' then
        Solve_by_lufco(As,bs,rcond,ys);
        put("rcond : "); put(rcond,3); new_line;
      else
        Solve_by_lufac(As,bs,info,ys);
        put("info : "); put(info,1); new_line;
      end if;
    end if;
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
    Write_Difference(xs,ys);
  end Standard_Test;

  procedure DoblDobl_Test ( n,m,d : in integer32 ) is

  -- DESCRIPTION :
  --   Generates and n-by-m matrix of series of degree d,
  --   with complex coefficients in double double precision.
  --   Converts an n-by-m matrix of series of degree d with double
  --   double precision complex coefficients into a matrix series.

    use DoblDobl_Dense_Series_Matrices;
    use DoblDobl_Matrix_Series_Solvers;

    sA : constant DoblDobl_Dense_Series_Matrices.Matrix(1..n,1..m)
       := DoblDobl_Random_Series.Random_Series_Matrix(1,n,1,m,d);
    As : constant DoblDobl_Dense_Matrix_Series.Matrix 
       := DoblDobl_Dense_Matrix_Series.Create(sA); 
    sx : constant DoblDobl_Dense_Series_Vectors.Vector(1..m)
       := DoblDobl_Random_Series.Random_Series_Vector(1,m,d);
    xs : constant DoblDobl_Dense_Vector_Series.Vector
       := DoblDobl_Dense_Vector_Series.Create(sx);
    sb : constant DoblDobl_Dense_Series_Vectors.Vector(1..n) := sA*sx;
    bs : constant DoblDobl_Dense_Vector_Series.Vector
       := DoblDobl_Dense_Vector_Series.Create(sb);
    ys : DoblDobl_Dense_Vector_Series.Vector;
    ans : character;
    rcond : double_double;
    info : integer32;

  begin
    put_line("The coefficients of the matrix series :"); put(As);
    put_line("The exact solution x :"); put_line(sx);
    put_line("The coefficients of the vector series x :"); put(xs);
    put_line("The right hand side vector b :"); put_line(sb);
    put_line("The coefficients of the vector series b :"); put(bs);
    if n > m then
      Solve_by_QRLS(As,bs,info,ys);
      put("info : "); put(info,1); new_line;
    else
      new_line;
      put("Condition number wanted ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        Solve_by_lufco(As,bs,rcond,ys);
        put("rcond : "); put(rcond,3); new_line;
      else
        Solve_by_lufac(As,bs,info,ys);
        put("info : "); put(info,1); new_line;
      end if;
    end if;
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
    Write_Difference(xs,ys);
  end DoblDobl_Test;

  procedure QuadDobl_Test ( n,m,d : in integer32 ) is

  -- DESCRIPTION :
  --   Generates and n-by-m matrix of series of degree d,
  --   with complex coefficients in quad double precision.
  --   Converts an n-by-m matrix of series of degree d with quad
  --   double precision complex coefficients into a matrix series.

    use QuadDobl_Dense_Series_Matrices;
    use QuadDobl_Matrix_Series_Solvers;

    sA : constant QuadDobl_Dense_Series_Matrices.Matrix(1..n,1..m)
       := QuadDobl_Random_Series.Random_Series_Matrix(1,n,1,m,d);
    As : constant QuadDobl_Dense_Matrix_Series.Matrix 
       := QuadDobl_Dense_Matrix_Series.Create(sA); 
    sx : constant QuadDobl_Dense_Series_Vectors.Vector(1..m)
       := QuadDobl_Random_Series.Random_Series_Vector(1,m,d);
    xs : constant QuadDobl_Dense_Vector_Series.Vector
       := QuadDobl_Dense_Vector_Series.Create(sx);
    sb : constant QuadDobl_Dense_Series_Vectors.Vector(1..n) := sA*sx;
    bs : constant QuadDobl_Dense_Vector_Series.Vector
       := QuadDobl_Dense_Vector_Series.Create(sb);
    ys : QuadDobl_Dense_Vector_Series.Vector;
    ans : character;
    rcond : quad_double;
    info : integer32;

  begin
    put_line("The coefficients of the matrix series :"); put(As);
    put_line("The exact solution x :"); put_line(sx);
    put_line("The coefficients of the vector series x :"); put(xs);
    put_line("The right hand side vector b :"); put_line(sb);
    put_line("The coefficients of the vector series b :"); put(bs);
    if n > m then
      Solve_by_QRLS(As,bs,info,ys);
      put("info : "); put(info,1); new_line;
    else
      new_line;
      put("Condition number wanted ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        Solve_by_lufco(As,bs,rcond,ys);
        put("rcond : "); put(rcond,3); new_line;
      else
        Solve_by_lufac(As,bs,info,ys);
        put("info : "); put(info,1); new_line;
      end if;
    end if;
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
    Write_Difference(xs,ys);
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

    neq,nvr,deg : integer32 := 0;
    prc : character;

  begin
    new_line;
    put_line("Testing the linearization of systems of power series ...");
    put("  Give the number of equations in the system : "); get(neq);
    put("  Give the number of variables in the system : "); get(nvr);
    put("  Give the degree of the series : "); get(deg);
    prc := Prompt_for_Precision;
    case prc is
      when '0' => Standard_Test(neq,nvr,deg);
      when '1' => DoblDobl_Test(neq,nvr,deg);
      when '2' => QuadDobl_Test(neq,nvr,deg);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_serlin;
