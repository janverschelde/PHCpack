with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Standard_Complex_Vector_Norms;
with Standard_Dense_Vector_Series2;
with Standard_Dense_Vector_Series2_io;   use Standard_Dense_Vector_Series2_io;
with Standard_Dense_Matrix_Series2;
with Standard_Dense_Matrix_Series2_io;   use Standard_Dense_Matrix_Series2_io;
with Standard_Interpolating_Series2;
with Random_Series_Vectors;
with Random_Series_Matrices;

procedure ts_seritp2 is

-- DESCRIPTION :
--   Development of solving linear systems of series with interpolation.

  procedure Standard_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension.
  --   Solves the problem in standard double precision.

    mat : constant Standard_Dense_Matrix_Series2.Matrix
        := Random_Series_Matrices.Random_Matrix_Series(deg,dim,0,1);
    sol : constant Standard_Dense_Vector_Series2.Vector
        := Random_Series_Vectors.Random_Vector_Series(1,dim,deg);
    rhs : constant Standard_Dense_Vector_Series2.Vector
        := Standard_Dense_Matrix_Series2.Multiply(mat,sol);
    rnk : integer32;
    cff : Standard_Dense_Vector_Series2.Vector(deg);

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
    rnk := Standard_Interpolating_Series2.Full_Rank(mat);
    if rnk = -1 then
      put_line("The matrix series does not have full rank.");
    else
      put_line("The matrix series has full rank.");
      put("The smallest degree for full rank : "); put(rnk,1); new_line;
      cff := Standard_Interpolating_Series2.Interpolate(mat,rhs);
      put_line("The computed solution :"); put(cff);
      put_line("The constructed solution :"); put(sol);
    end if;
  end Standard_Test;

  function Standard_Special_Matrix_Series
             ( deg : integer32 ) 
             return Standard_Dense_Matrix_Series2.Matrix is

  -- DESCRIPTION :
  --   Returns a special matrix series of degree equal to deg,
  --   which has an infinite series as solution when the
  --   right hand side vector is (1, 0).
  --   The degree should be at least one.

    res : Standard_Dense_Matrix_Series2.Matrix(deg);
    wrk : Standard_Complex_Matrices.Matrix(1..2,1..2);

  begin
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        wrk(i,j) := Standard_Complex_Numbers.Create(0.0);
      end loop;
    end loop;
    for k in 0..deg loop
      res.cff(k) := new Standard_Complex_Matrices.Matrix'(wrk);
    end loop;
    res.cff(0)(1,1) := Standard_Complex_Numbers.Create(1.0);
    res.cff(0)(1,2) := Standard_Complex_Numbers.Create(1.0);
    res.cff(0)(2,2) := Standard_Complex_Numbers.Create(1.0);
    res.cff(1)(2,1) := Standard_Complex_Numbers.Create(1.0);
    return res;
  end Standard_Special_Matrix_Series;

  function Standard_Special_Vector_Series
             ( deg : integer32 )
             return Standard_Dense_Vector_Series2.Vector is

  -- DESCRIPTION :
  --   Returns a special right hand side series with degree equal to deg.

    res : Standard_Dense_Vector_Series2.Vector(deg);
    wrk : Standard_Complex_Vectors.Vector(1..2);

  begin
    for i in wrk'range loop
      wrk(i) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    for k in 0..deg loop
      res.cff(k) := new Standard_Complex_Vectors.Vector'(wrk);
    end loop;
    res.cff(0)(1) := Standard_Complex_Numbers.Create(1.0);
    return res;
  end Standard_Special_Vector_Series;

  function Standard_Singular_Matrix_Series
             return Standard_Dense_Matrix_Series2.Matrix is

  -- DESCRIPTION :
  --   Returns a singular matrix series of degree of degree one,
  --   in standard double precision.

    res : Standard_Dense_Matrix_Series2.Matrix(1);
    wrk : Standard_Complex_Matrices.Matrix(1..2,1..2);

  begin
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        wrk(i,j) := Standard_Complex_Numbers.Create(0.0);
      end loop;
    end loop;
    for k in 0..res.deg loop
      res.cff(k) := new Standard_Complex_Matrices.Matrix'(wrk);
    end loop;
    res.cff(0)(1,2) := Standard_Complex_Numbers.Create(4.0);
    res.cff(1)(1,1) := Standard_Complex_Numbers.Create(4.0);
    res.cff(1)(2,1) := Standard_Complex_Numbers.Create(4.0);
    return res;
  end Standard_Singular_Matrix_Series;

  function Standard_Singular_Vector_Series
             return Standard_Dense_Vector_Series2.Vector is

  -- DESCRIPTION :
  --   Returns a special right hand side series with degree equal to one,
  --   in standard double precision.

    res : Standard_Dense_Vector_Series2.Vector(1);
    wrk : Standard_Complex_Vectors.Vector(1..2);

  begin
    for i in wrk'range loop
      wrk(i) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    for k in 0..res.deg loop
      res.cff(k) := new Standard_Complex_Vectors.Vector'(wrk);
    end loop;
    res.cff(0)(1) := Standard_Complex_Numbers.Create(1.0);
    res.cff(0)(2) := Standard_Complex_Numbers.Create(0.0);
    res.cff(1)(1) := Standard_Complex_Numbers.Create(0.0);
    res.cff(1)(2) := Standard_Complex_Numbers.Create(0.0);
    return res;
  end Standard_Singular_Vector_Series;

  procedure Write_Differences
              ( x,y : in Standard_Dense_Vector_Series2.Vector ) is

  -- DESCRIPTION :
  --  Writes the differences between the series x and y.

    use Standard_Complex_Vectors;

    deg : integer32;
    dff : Standard_Complex_Vectors.Vector(x.cff(0)'range);
    nrm : double_float;

  begin
    if x.deg <= y.deg
     then deg := x.deg;
     else deg := y.deg;
    end if;
    for k in 0..deg loop
      put("x coefficient "); put(k,1); put_line(" :");
      put_line(x.cff(k));
      put("y coefficient "); put(k,1); put_line(" :");
      put_line(y.cff(k));
      dff := x.cff(k).all - y.cff(k).all;
      nrm := Standard_Complex_Vector_Norms.Max_Norm(dff);
      put("Max norm of the difference : "); put(nrm,3); new_line;
    end loop;
  end Write_Differences;

  procedure Compare_Values
              ( x,y : in Standard_Dense_Vector_Series2.Vector;
                t : in double_float ) is

  -- DESCRIPTION :
  --   Evaluates x and y at t and writes their values.

    xt : constant Standard_Complex_Vectors.Vector
       := Standard_Dense_Vector_Series2.Eval(x,t);
    yt : constant Standard_Complex_Vectors.Vector
       := Standard_Dense_Vector_Series2.Eval(y,t);

  begin
    put("x evaluated at "); put(t); put_line(" :");
    put_line(xt);
    put("y evaluated at "); put(t); put_line(" :");
    put_line(yt);
  end Compare_Values;

  procedure Standard_Random_Hermite_Test
              ( deg,dim : in integer32; laurent : in boolean ) is

  -- DESCRIPTION :
  --   Generates random data of a degree larger than deg
  --   for a problem of the given dimension dim.
  --   If laurent, then Hermite-Laurent interpolation is used.

    use Standard_Complex_Numbers;
    use Standard_Interpolating_Series2;

    mat : Standard_Dense_Matrix_Series2.Matrix
        := Random_Series_Matrices.Random_Matrix_Series(deg,dim,0,1);
    sol : constant Standard_Dense_Vector_Series2.Vector
        := Random_Series_Vectors.Random_Vector_Series(1,dim,deg);
    rhs : constant Standard_Dense_Vector_Series2.Vector
        := Standard_Dense_Matrix_Series2.Multiply(mat,sol);
    rnd : constant Standard_Complex_Numbers.Complex_Number
        := Standard_Random_Numbers.Random1;
    t : Standard_Complex_Numbers.Complex_Number
      := Standard_Complex_Numbers.Create(0.01)*rnd;
    x : Standard_Dense_Vector_Series2.Vector(deg);
    rnk : integer32;

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
   -- mat.deg := deg; rhs.deg := deg; -- truncate the degrees
    rnk := Standard_Interpolating_Series2.Full_Rank(mat);
    if rnk = -1 then
      put_line("The matrix series does not have full rank.");
    else
      put_line("The matrix series has full rank.");
      put("The smallest degree for full rank : "); put(rnk,1); new_line;
      if laurent 
       then x := Hermite_Laurent_Interpolate(mat,rhs);
       else x := Hermite_Interpolate(mat,rhs,t);
      end if;
      put_line("The computed solution :"); put(x);
      Write_Differences(x,sol);
      Compare_Values(x,sol,0.01);
    end if;
  end Standard_Random_Hermite_Test;

  procedure Standard_Special_Hermite_Test ( deg : integer32 ) is

  -- DESCRIPTION :
  --   Applies Hermite interpolation on a special case,
  --   of dimension two, in standard double precision.

    mat : constant Standard_Dense_Matrix_Series2.Matrix
        := Standard_Special_Matrix_Series(deg);
    rhs : constant Standard_Dense_Vector_Series2.Vector
        := Standard_Special_Vector_Series(deg);
    sol : Standard_Dense_Vector_Series2.Vector(deg);
    t : Standard_Complex_Numbers.Complex_Number
      := Standard_Complex_Numbers.Create(0.01);

  begin
    put_line("The special matrix series : "); put(mat);
    put_line("The special vector series : "); put(rhs);
    sol := Standard_Interpolating_Series2.Hermite_Interpolate(mat,rhs,t);
    put_line("The solution : "); put(sol);
  end Standard_Special_Hermite_Test;

  procedure Standard_Singular_Hermite_Test is

  -- DESCRIPTION :
  --   Applies Hermite interpolation on a singular case,
  --   of dimension two, in standard double precision.

    mat : constant Standard_Dense_Matrix_Series2.Matrix
        := Standard_Singular_Matrix_Series;
    rhs : constant Standard_Dense_Vector_Series2.Vector
        := Standard_Singular_Vector_Series;
    sol,y : Standard_Dense_Vector_Series2.Vector(mat.deg);

    use Standard_Interpolating_Series2;

  begin
    put_line("The singular matrix series : "); put(mat);
    put_line("The singular vector series : "); put(rhs);
    sol := Hermite_Laurent_Interpolate(mat,rhs);
    put_line("The solution : "); put(sol);
    y := Standard_Dense_Matrix_Series2.Multiply(mat,sol);
    put_line("The matrix multiplied with the solution :"); put(y);
  end Standard_Singular_Hermite_Test;

  procedure Standard_Hermite_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   To test the Hermite interpolation, we have a special case
  --   and the general case, tested on a randomly generated problem.
  --   Hermite interpolation is applied in standard double precision.

    lau,spc : character;

  begin
    put("Apply Hermite-Laurent interpolation ? (y/n) ");
    Ask_Yes_or_No(lau);
    put("Test special or singular case ? (y/n) ");
    Ask_Yes_or_No(spc);
    new_line;
    if spc = 'y' then
      if lau = 'y'
       then Standard_Singular_Hermite_Test;
       else Standard_Special_Hermite_Test(deg);
      end if;
     else Standard_Random_Hermite_Test(deg,dim,lau='y');
    end if;
  end Standard_Hermite_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree and the dimension of the series.
  --   The degree is the highest exponent in the power series.
  --   The dimension is the number of variables in the series.

    deg,dim : integer32 := 0;
    hrm : character;

  begin
    new_line;
    put_line("Testing solving linear systems with interpolation.");
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    put("Hermite interpolation ? (y/n) ");
    Ask_Yes_or_No(hrm);
    if hrm = 'y' then
      Standard_Hermite_Test(deg,dim);
    else
      new_line;
      Standard_Test(deg,dim);
    end if;
  end Main;

begin
  Main;
end ts_seritp2;
