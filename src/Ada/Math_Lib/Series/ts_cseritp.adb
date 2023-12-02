with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Vector_Series;
with Standard_Complex_Vector_Series_io;  use Standard_Complex_Vector_Series_io;
with Standard_Complex_Matrix_Series;
with Standard_Complex_Matrix_Series_io;  use Standard_Complex_Matrix_Series_io;
with Standard_Random_Series_Vectors;
with Standard_Random_Series_Matrices;
with DoblDobl_Complex_Vector_Series;
with DoblDobl_Complex_Vector_Series_io;  use DoblDobl_Complex_Vector_Series_io;
with DoblDobl_Complex_Matrix_Series;
with DoblDobl_Complex_Matrix_Series_io;  use DoblDobl_Complex_Matrix_Series_io;
with DoblDobl_Random_Series_Vectors;
with DoblDobl_Random_Series_Matrices;
with QuadDobl_Complex_Vector_Series;
with QuadDobl_Complex_Vector_Series_io;  use QuadDobl_Complex_Vector_Series_io;
with QuadDobl_Complex_Matrix_Series;
with QuadDobl_Complex_Matrix_Series_io;  use QuadDobl_Complex_Matrix_Series_io;
with QuadDobl_Random_Series_Vectors;
with QuadDobl_Random_Series_Matrices;
with Standard_Interpolating_CSeries;
with DoblDobl_Interpolating_CSeries;
with QuadDobl_Interpolating_CSeries;

procedure ts_cseritp is

-- DESCRIPTION :
--   Development of solving linear systems of series with interpolation.

  procedure Standard_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension.
  --   Solves the problem in standard double precision.

    mat : constant Standard_Complex_Matrix_Series.Matrix
        := Standard_Random_Series_Matrices.Random_Matrix_Series(deg,dim,0,1);
    sol : constant Standard_Complex_Vector_Series.Vector
        := Standard_Random_Series_Vectors.Random_Vector_Series(1,dim,deg);
    rhs : constant Standard_Complex_Vector_Series.Vector
        := Standard_Complex_Matrix_Series.Multiply(mat,sol);
    rnk : integer32;
    cff : Standard_Complex_Vector_Series.Vector(deg);

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
    rnk := Standard_Interpolating_CSeries.Full_Rank(mat);
    if rnk = -1 then
      put_line("The matrix series does not have full rank.");
    else
      put_line("The matrix series has full rank.");
      put("The smallest degree for full rank : "); put(rnk,1); new_line;
      cff := Standard_Interpolating_CSeries.Interpolate(mat,rhs);
      put_line("The computed solution :"); put(cff);
      put_line("The constructed solution :"); put(sol);
    end if;
  end Standard_Test;

  procedure DoblDobl_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension.
  --   Solves the problem in double double precision.

    mat : constant DoblDobl_Complex_Matrix_Series.Matrix
        := DoblDobl_Random_Series_Matrices.Random_Matrix_Series(deg,dim,0,1);
    sol : constant DoblDobl_Complex_Vector_Series.Vector
        := DoblDobl_Random_Series_Vectors.Random_Vector_Series(1,dim,deg);
    rhs : constant DoblDobl_Complex_Vector_Series.Vector
        := DoblDobl_Complex_Matrix_Series.Multiply(mat,sol);
    rnk : integer32;
    cff : DoblDobl_Complex_Vector_Series.Vector(deg);

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
    rnk := DoblDobl_Interpolating_CSeries.Full_Rank(mat);
    if rnk = -1 then
      put_line("The matrix series does not have full rank.");
    else
      put_line("The matrix series has full rank.");
      put("The smallest degree for full rank : "); put(rnk,1); new_line;
      cff := DoblDobl_Interpolating_CSeries.Interpolate(mat,rhs);
      put_line("The computed solution :"); put(cff);
      put_line("The constructed solution :"); put(sol);
    end if;
  end DoblDobl_Test;

  procedure QuadDobl_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension.
  --   Solves the problem in quad double precision.

    mat : constant QuadDobl_Complex_Matrix_Series.Matrix
        := QuadDobl_Random_Series_Matrices.Random_Matrix_Series(deg,dim,0,1);
    sol : constant QuadDobl_Complex_Vector_Series.Vector
        := QuadDobl_Random_Series_Vectors.Random_Vector_Series(1,dim,deg);
    rhs : constant QuadDobl_Complex_Vector_Series.Vector
        := QuadDobl_Complex_Matrix_Series.Multiply(mat,sol);
    rnk : integer32;
    cff : QuadDobl_Complex_Vector_Series.Vector(deg);

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
    rnk := QuadDobl_Interpolating_CSeries.Full_Rank(mat);
    if rnk = -1 then
      put_line("The matrix series does not have full rank.");
    else
      put_line("The matrix series has full rank.");
      put("The smallest degree for full rank : "); put(rnk,1); new_line;
      cff := QuadDobl_Interpolating_CSeries.Interpolate(mat,rhs);
      put_line("The computed solution :"); put(cff);
      put_line("The constructed solution :"); put(sol);
    end if;
  end QuadDobl_Test;

  function Standard_Special_Matrix_Series
             ( deg : integer32 ) 
             return Standard_Complex_Matrix_Series.Matrix is

  -- DESCRIPTION :
  --   Returns a special matrix series of degree equal to deg,
  --   which has an infinite series as solution when the
  --   right hand side vector is (1, 0).
  --   The degree should be at least one.

    res : Standard_Complex_Matrix_Series.Matrix(deg);
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

  function DoblDobl_Special_Matrix_Series
             ( deg : integer32 ) 
             return DoblDobl_Complex_Matrix_Series.Matrix is

  -- DESCRIPTION :
  --   Returns a special matrix series of degree equal to deg,
  --   which has an infinite series as solution when the
  --   right hand side vector is (1, 0).
  --   The degree should be at least one.

    res : DoblDobl_Complex_Matrix_Series.Matrix(deg);
    wrk : DoblDobl_Complex_Matrices.Matrix(1..2,1..2);
    one : constant double_double := create(1.0);
    zero : constant double_double := create(0.0);

  begin
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        wrk(i,j) := DoblDobl_Complex_Numbers.Create(zero);
      end loop;
    end loop;
    for k in 0..deg loop
      res.cff(k) := new DoblDobl_Complex_Matrices.Matrix'(wrk);
    end loop;
    res.cff(0)(1,1) := DoblDobl_Complex_Numbers.Create(one);
    res.cff(0)(1,2) := DoblDobl_Complex_Numbers.Create(one);
    res.cff(0)(2,2) := DoblDobl_Complex_Numbers.Create(one);
    res.cff(1)(2,1) := DoblDobl_Complex_Numbers.Create(one);
    return res;
  end DoblDobl_Special_Matrix_Series;

  function QuadDobl_Special_Matrix_Series
             ( deg : integer32 ) 
             return QuadDobl_Complex_Matrix_Series.Matrix is

  -- DESCRIPTION :
  --   Returns a special matrix series of degree equal to deg,
  --   which has an infinite series as solution when the
  --   right hand side vector is (1, 0).
  --   The degree should be at least one.

    res : QuadDobl_Complex_Matrix_Series.Matrix(deg);
    wrk : QuadDobl_Complex_Matrices.Matrix(1..2,1..2);
    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);

  begin
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        wrk(i,j) := QuadDobl_Complex_Numbers.Create(zero);
      end loop;
    end loop;
    for k in 0..deg loop
      res.cff(k) := new QuadDobl_Complex_Matrices.Matrix'(wrk);
    end loop;
    res.cff(0)(1,1) := QuadDobl_Complex_Numbers.Create(one);
    res.cff(0)(1,2) := QuadDobl_Complex_Numbers.Create(one);
    res.cff(0)(2,2) := QuadDobl_Complex_Numbers.Create(one);
    res.cff(1)(2,1) := QuadDobl_Complex_Numbers.Create(one);
    return res;
  end QuadDobl_Special_Matrix_Series;

  function Standard_Special_Vector_Series
             ( deg : integer32 )
             return Standard_Complex_Vector_Series.Vector is

  -- DESCRIPTION :
  --   Returns a special right hand side series with degree equal to deg.

    res : Standard_Complex_Vector_Series.Vector(deg);
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

  function DoblDobl_Special_Vector_Series
             ( deg : integer32 )
             return DoblDobl_Complex_Vector_Series.Vector is

  -- DESCRIPTION :
  --   Returns a special right hand side series with degree equal to deg.

    res : DoblDobl_Complex_Vector_Series.Vector(deg);
    wrk : DoblDobl_Complex_Vectors.Vector(1..2);
    one : constant double_double := create(0.0);
    zero : constant double_double := create(1.0);

  begin
    for i in wrk'range loop
      wrk(i) := DoblDobl_Complex_Numbers.Create(zero);
    end loop;
    for k in 0..deg loop
      res.cff(k) := new DoblDobl_Complex_Vectors.Vector'(wrk);
    end loop;
    res.cff(0)(1) := DoblDobl_Complex_Numbers.Create(one);
    return res;
  end DoblDobl_Special_Vector_Series;

  function QuadDobl_Special_Vector_Series
             ( deg : integer32 )
             return QuadDobl_Complex_Vector_Series.Vector is

  -- DESCRIPTION :
  --   Returns a special right hand side series with degree equal to deg.

    res : QuadDobl_Complex_Vector_Series.Vector(deg);
    wrk : QuadDobl_Complex_Vectors.Vector(1..2);
    one : constant quad_double := create(0.0);
    zero : constant quad_double := create(1.0);

  begin
    for i in wrk'range loop
      wrk(i) := QuadDobl_Complex_Numbers.Create(zero);
    end loop;
    for k in 0..deg loop
      res.cff(k) := new QuadDobl_Complex_Vectors.Vector'(wrk);
    end loop;
    res.cff(0)(1) := QuadDobl_Complex_Numbers.Create(one);
    return res;
  end QuadDobl_Special_Vector_Series;

  function Standard_Singular_Matrix_Series
             return Standard_Complex_Matrix_Series.Matrix is

  -- DESCRIPTION :
  --   Returns a singular matrix series of degree of degree one,
  --   in standard double precision.

    res : Standard_Complex_Matrix_Series.Matrix(1);
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

  function DoblDobl_Singular_Matrix_Series
             return DoblDobl_Complex_Matrix_Series.Matrix is

  -- DESCRIPTION :
  --   Returns a singular matrix series of degree of degree one,
  --   in double double precision.

    res : DoblDobl_Complex_Matrix_Series.Matrix(1);
    wrk : DoblDobl_Complex_Matrices.Matrix(1..2,1..2);
    zero : constant double_double := create(0.0);
    four : constant double_double := create(4.0);

  begin
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        wrk(i,j) := DoblDobl_Complex_Numbers.Create(zero);
      end loop;
    end loop;
    for k in 0..res.deg loop
      res.cff(k) := new DoblDobl_Complex_Matrices.Matrix'(wrk);
    end loop;
    res.cff(0)(1,2) := DoblDobl_Complex_Numbers.Create(four);
    res.cff(1)(1,1) := DoblDobl_Complex_Numbers.Create(four);
    res.cff(1)(2,1) := DoblDobl_Complex_Numbers.Create(four);
    return res;
  end DoblDobl_Singular_Matrix_Series;

  function QuadDobl_Singular_Matrix_Series
             return QuadDobl_Complex_Matrix_Series.Matrix is

  -- DESCRIPTION :
  --   Returns a singular matrix series of degree of degree one,
  --   in double double precision.

    res : QuadDobl_Complex_Matrix_Series.Matrix(1);
    wrk : QuadDobl_Complex_Matrices.Matrix(1..2,1..2);
    zero : constant quad_double := create(0.0);
    four : constant quad_double := create(4.0);

  begin
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        wrk(i,j) := QuadDobl_Complex_Numbers.Create(zero);
      end loop;
    end loop;
    for k in 0..res.deg loop
      res.cff(k) := new QuadDobl_Complex_Matrices.Matrix'(wrk);
    end loop;
    res.cff(0)(1,2) := QuadDobl_Complex_Numbers.Create(four);
    res.cff(1)(1,1) := QuadDobl_Complex_Numbers.Create(four);
    res.cff(1)(2,1) := QuadDobl_Complex_Numbers.Create(four);
    return res;
  end QuadDobl_Singular_Matrix_Series;

  function Standard_Singular_Vector_Series
             return Standard_Complex_Vector_Series.Vector is

  -- DESCRIPTION :
  --   Returns a special right hand side series with degree equal to one,
  --   in standard double precision.

    res : Standard_Complex_Vector_Series.Vector(1);
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

  function DoblDobl_Singular_Vector_Series
             return DoblDobl_Complex_Vector_Series.Vector is

  -- DESCRIPTION :
  --   Returns a special right hand side series with degree equal to one,
  --   in double double precision.

    res : DoblDobl_Complex_Vector_Series.Vector(1);
    wrk : DoblDobl_Complex_Vectors.Vector(1..2);
    one : constant double_double := create(1.0);
    zero : constant double_double := create(0.0);

  begin
    for i in wrk'range loop
      wrk(i) := DoblDobl_Complex_Numbers.Create(zero);
    end loop;
    for k in 0..res.deg loop
      res.cff(k) := new DoblDobl_Complex_Vectors.Vector'(wrk);
    end loop;
    res.cff(0)(1) := DoblDobl_Complex_Numbers.Create(one);
    res.cff(0)(2) := DoblDobl_Complex_Numbers.Create(zero);
    res.cff(1)(1) := DoblDobl_Complex_Numbers.Create(zero);
    res.cff(1)(2) := DoblDobl_Complex_Numbers.Create(zero);
    return res;
  end DoblDobl_Singular_Vector_Series;

  function QuadDobl_Singular_Vector_Series
             return QuadDobl_Complex_Vector_Series.Vector is

  -- DESCRIPTION :
  --   Returns a special right hand side series with degree equal to one,
  --   in quad double precision.

    res : QuadDobl_Complex_Vector_Series.Vector(1);
    wrk : QuadDobl_Complex_Vectors.Vector(1..2);
    one : constant quad_double := create(1.0);
    zero : constant quad_double := create(0.0);

  begin
    for i in wrk'range loop
      wrk(i) := QuadDobl_Complex_Numbers.Create(zero);
    end loop;
    for k in 0..res.deg loop
      res.cff(k) := new QuadDobl_Complex_Vectors.Vector'(wrk);
    end loop;
    res.cff(0)(1) := QuadDobl_Complex_Numbers.Create(one);
    res.cff(0)(2) := QuadDobl_Complex_Numbers.Create(zero);
    res.cff(1)(1) := QuadDobl_Complex_Numbers.Create(zero);
    res.cff(1)(2) := QuadDobl_Complex_Numbers.Create(zero);
    return res;
  end QuadDobl_Singular_Vector_Series;

  procedure Write_Differences
              ( x,y : in Standard_Complex_Vector_Series.Vector ) is

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

  procedure Write_Differences
              ( x,y : in DoblDobl_Complex_Vector_Series.Vector ) is

  -- DESCRIPTION :
  --  Writes the differences between the series x and y.

    use DoblDobl_Complex_Vectors;

    deg : integer32;
    dff : DoblDobl_Complex_Vectors.Vector(x.cff(0)'range);
    nrm : double_double;

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
      nrm := DoblDobl_Complex_Vector_Norms.Max_Norm(dff);
      put("Max norm of the difference : "); put(nrm,3); new_line;
    end loop;
  end Write_Differences;

  procedure Write_Differences
              ( x,y : in QuadDobl_Complex_Vector_Series.Vector ) is

  -- DESCRIPTION :
  --  Writes the differences between the series x and y.

    use QuadDobl_Complex_Vectors;

    deg : integer32;
    dff : QuadDobl_Complex_Vectors.Vector(x.cff(0)'range);
    nrm : quad_double;

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
      nrm := QuadDobl_Complex_Vector_Norms.Max_Norm(dff);
      put("Max norm of the difference : "); put(nrm,3); new_line;
    end loop;
  end Write_Differences;

  procedure Compare_Values
              ( x,y : in Standard_Complex_Vector_Series.Vector;
                t : in double_float ) is

  -- DESCRIPTION :
  --   Evaluates x and y at t and writes their values.

    xt : constant Standard_Complex_Vectors.Vector
       := Standard_Complex_Vector_Series.Eval(x,t);
    yt : constant Standard_Complex_Vectors.Vector
       := Standard_Complex_Vector_Series.Eval(y,t);

  begin
    put("x evaluated at "); put(t); put_line(" :"); put_line(xt);
    put("y evaluated at "); put(t); put_line(" :"); put_line(yt);
  end Compare_Values;

  procedure Compare_Values
              ( x,y : in DoblDobl_Complex_Vector_Series.Vector;
                t : in double_double ) is

  -- DESCRIPTION :
  --   Evaluates x and y at t and writes their values.

    xt : constant DoblDobl_Complex_Vectors.Vector
       := DoblDobl_Complex_Vector_Series.Eval(x,t);
    yt : constant DoblDobl_Complex_Vectors.Vector
       := DoblDobl_Complex_Vector_Series.Eval(y,t);

  begin
    put("x evaluated at "); put(t); put_line(" :"); put_line(xt);
    put("y evaluated at "); put(t); put_line(" :"); put_line(yt);
  end Compare_Values;

  procedure Compare_Values
              ( x,y : in QuadDobl_Complex_Vector_Series.Vector;
                t : in quad_double ) is

  -- DESCRIPTION :
  --   Evaluates x and y at t and writes their values.

    xt : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Complex_Vector_Series.Eval(x,t);
    yt : constant QuadDobl_Complex_Vectors.Vector
       := QuadDobl_Complex_Vector_Series.Eval(y,t);

  begin
    put("x evaluated at "); put(t); put_line(" :"); put_line(xt);
    put("y evaluated at "); put(t); put_line(" :"); put_line(yt);
  end Compare_Values;

  procedure Standard_Random_Hermite_Test
              ( deg,dim : in integer32; laurent : in boolean ) is

  -- DESCRIPTION :
  --   Generates random data of a degree larger than deg
  --   for a problem of the given dimension dim.
  --   If laurent, then Hermite-Laurent interpolation is used.

    use Standard_Complex_Numbers;
    use Standard_Interpolating_CSeries;

    mat : constant Standard_Complex_Matrix_Series.Matrix
        := Standard_Random_Series_Matrices.Random_Matrix_Series(deg,dim,0,1);
    sol : constant Standard_Complex_Vector_Series.Vector
        := Standard_Random_Series_Vectors.Random_Vector_Series(1,dim,deg);
    rhs : constant Standard_Complex_Vector_Series.Vector
        := Standard_Complex_Matrix_Series.Multiply(mat,sol);
    rnd : constant Standard_Complex_Numbers.Complex_Number
        := Standard_Random_Numbers.Random1;
    t : constant Standard_Complex_Numbers.Complex_Number
      := Standard_Complex_Numbers.Create(0.01)*rnd;
    x : Standard_Complex_Vector_Series.Vector(deg);
    rnk : integer32;

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
   -- mat.deg := deg; rhs.deg := deg; -- truncate the degrees
    rnk := Standard_Interpolating_CSeries.Full_Rank(mat);
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

  procedure DoblDobl_Random_Hermite_Test
              ( deg,dim : in integer32; laurent : in boolean ) is

  -- DESCRIPTION :
  --   Generates random data of a degree larger than deg
  --   for a problem of the given dimension dim.
  --   If laurent, then Hermite-Laurent interpolation is used.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Interpolating_CSeries;

    mat : constant DoblDobl_Complex_Matrix_Series.Matrix
        := DoblDobl_Random_Series_Matrices.Random_Matrix_Series(deg,dim,0,1);
    sol : constant DoblDobl_Complex_Vector_Series.Vector
        := DoblDobl_Random_Series_Vectors.Random_Vector_Series(1,dim,deg);
    rhs : constant DoblDobl_Complex_Vector_Series.Vector
        := DoblDobl_Complex_Matrix_Series.Multiply(mat,sol);
    rnd : constant DoblDobl_Complex_Numbers.Complex_Number
        := DoblDobl_Random_Numbers.Random1;
    val : constant double_double := create(0.01);
    t : constant DoblDobl_Complex_Numbers.Complex_Number
      := DoblDobl_Complex_Numbers.Create(val)*rnd;
    x : DoblDobl_Complex_Vector_Series.Vector(deg);
    rnk : integer32;

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
   -- mat.deg := deg; rhs.deg := deg; -- truncate the degrees
    rnk := DoblDobl_Interpolating_CSeries.Full_Rank(mat);
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
      Compare_Values(x,sol,val);
    end if;
  end DoblDobl_Random_Hermite_Test;

  procedure QuadDobl_Random_Hermite_Test
              ( deg,dim : in integer32; laurent : in boolean ) is

  -- DESCRIPTION :
  --   Generates random data of a degree larger than deg
  --   for a problem of the given dimension dim.
  --   If laurent, then Hermite-Laurent interpolation is used.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Interpolating_CSeries;

    mat : constant QuadDobl_Complex_Matrix_Series.Matrix
        := QuadDobl_Random_Series_Matrices.Random_Matrix_Series(deg,dim,0,1);
    sol : constant QuadDobl_Complex_Vector_Series.Vector
        := QuadDobl_Random_Series_Vectors.Random_Vector_Series(1,dim,deg);
    rhs : constant QuadDobl_Complex_Vector_Series.Vector
        := QuadDobl_Complex_Matrix_Series.Multiply(mat,sol);
    rnd : constant QuadDobl_Complex_Numbers.Complex_Number
        := QuadDobl_Random_Numbers.Random1;
    val : constant quad_double := create(0.01);
    t : constant QuadDobl_Complex_Numbers.Complex_Number
      := QuadDobl_Complex_Numbers.Create(val)*rnd;
    x : QuadDobl_Complex_Vector_Series.Vector(deg);
    rnk : integer32;

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
   -- mat.deg := deg; rhs.deg := deg; -- truncate the degrees
    rnk := QuadDobl_Interpolating_CSeries.Full_Rank(mat);
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
      Compare_Values(x,sol,val);
    end if;
  end QuadDobl_Random_Hermite_Test;

  procedure Standard_Special_Hermite_Test ( deg : integer32 ) is

  -- DESCRIPTION :
  --   Applies Hermite interpolation on a special case,
  --   of dimension two, in standard double precision.

    mat : constant Standard_Complex_Matrix_Series.Matrix
        := Standard_Special_Matrix_Series(deg);
    rhs : constant Standard_Complex_Vector_Series.Vector
        := Standard_Special_Vector_Series(deg);
    sol : Standard_Complex_Vector_Series.Vector(deg);
    t : constant Standard_Complex_Numbers.Complex_Number
      := Standard_Complex_Numbers.Create(0.01);

  begin
    put_line("The special matrix series : "); put(mat);
    put_line("The special vector series : "); put(rhs);
    sol := Standard_Interpolating_CSeries.Hermite_Interpolate(mat,rhs,t);
    put_line("The solution : "); put(sol);
  end Standard_Special_Hermite_Test;

  procedure DoblDobl_Special_Hermite_Test ( deg : integer32 ) is

  -- DESCRIPTION :
  --   Applies Hermite interpolation on a special case,
  --   of dimension two, in double double precision.

    mat : constant DoblDobl_Complex_Matrix_Series.Matrix
        := DoblDobl_Special_Matrix_Series(deg);
    rhs : constant DoblDobl_Complex_Vector_Series.Vector
        := DoblDobl_Special_Vector_Series(deg);
    sol : DoblDobl_Complex_Vector_Series.Vector(deg);
    val : constant double_double := create(0.01);
    t : constant DoblDobl_Complex_Numbers.Complex_Number
      := DoblDobl_Complex_Numbers.Create(val);

  begin
    put_line("The special matrix series : "); put(mat);
    put_line("The special vector series : "); put(rhs);
    sol := DoblDobl_Interpolating_CSeries.Hermite_Interpolate(mat,rhs,t);
    put_line("The solution : "); put(sol);
  end DoblDobl_Special_Hermite_Test;

  procedure QuadDobl_Special_Hermite_Test ( deg : integer32 ) is

  -- DESCRIPTION :
  --   Applies Hermite interpolation on a special case,
  --   of dimension two, in double double precision.

    mat : constant QuadDobl_Complex_Matrix_Series.Matrix
        := QuadDobl_Special_Matrix_Series(deg);
    rhs : constant QuadDobl_Complex_Vector_Series.Vector
        := QuadDobl_Special_Vector_Series(deg);
    sol : QuadDobl_Complex_Vector_Series.Vector(deg);
    val : constant quad_double := create(0.01);
    t : constant QuadDobl_Complex_Numbers.Complex_Number
      := QuadDobl_Complex_Numbers.Create(val);

  begin
    put_line("The special matrix series : "); put(mat);
    put_line("The special vector series : "); put(rhs);
    sol := QuadDobl_Interpolating_CSeries.Hermite_Interpolate(mat,rhs,t);
    put_line("The solution : "); put(sol);
  end QuadDobl_Special_Hermite_Test;

  procedure Standard_Singular_Hermite_Test is

  -- DESCRIPTION :
  --   Applies Hermite interpolation on a singular case,
  --   of dimension two, in standard double precision.

    mat : constant Standard_Complex_Matrix_Series.Matrix
        := Standard_Singular_Matrix_Series;
    rhs : constant Standard_Complex_Vector_Series.Vector
        := Standard_Singular_Vector_Series;
    sol,y : Standard_Complex_Vector_Series.Vector(mat.deg);

    use Standard_Interpolating_CSeries;

  begin
    put_line("The singular matrix series : "); put(mat);
    put_line("The singular vector series : "); put(rhs);
    sol := Hermite_Laurent_Interpolate(mat,rhs);
    put_line("The solution : "); put(sol);
    y := Standard_Complex_Matrix_Series.Multiply(mat,sol);
    put_line("The matrix multiplied with the solution :"); put(y);
  end Standard_Singular_Hermite_Test;

  procedure DoblDobl_Singular_Hermite_Test is

  -- DESCRIPTION :
  --   Applies Hermite interpolation on a singular case,
  --   of dimension two, in double double precision.

    mat : constant DoblDobl_Complex_Matrix_Series.Matrix
        := DoblDobl_Singular_Matrix_Series;
    rhs : constant DoblDobl_Complex_Vector_Series.Vector
        := DoblDobl_Singular_Vector_Series;
    sol,y : DoblDobl_Complex_Vector_Series.Vector(mat.deg);

    use DoblDobl_Interpolating_CSeries;

  begin
    put_line("The singular matrix series : "); put(mat);
    put_line("The singular vector series : "); put(rhs);
    sol := Hermite_Laurent_Interpolate(mat,rhs);
    put_line("The solution : "); put(sol);
    y := DoblDobl_Complex_Matrix_Series.Multiply(mat,sol);
    put_line("The matrix multiplied with the solution :"); put(y);
  end DoblDobl_Singular_Hermite_Test;

  procedure QuadDobl_Singular_Hermite_Test is

  -- DESCRIPTION :
  --   Applies Hermite interpolation on a singular case,
  --   of dimension two, in quad double precision.

    mat : constant QuadDobl_Complex_Matrix_Series.Matrix
        := QuadDobl_Singular_Matrix_Series;
    rhs : constant QuadDobl_Complex_Vector_Series.Vector
        := QuadDobl_Singular_Vector_Series;
    sol,y : QuadDobl_Complex_Vector_Series.Vector(mat.deg);

    use QuadDobl_Interpolating_CSeries;

  begin
    put_line("The singular matrix series : "); put(mat);
    put_line("The singular vector series : "); put(rhs);
    sol := Hermite_Laurent_Interpolate(mat,rhs);
    put_line("The solution : "); put(sol);
    y := QuadDobl_Complex_Matrix_Series.Multiply(mat,sol);
    put_line("The matrix multiplied with the solution :"); put(y);
  end QuadDobl_Singular_Hermite_Test;

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

  procedure DoblDobl_Hermite_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   To test the Hermite interpolation, we have a special case
  --   and the general case, tested on a randomly generated problem.
  --   Hermite interpolation is applied in double double precision.

    lau,spc : character;

  begin
    put("Apply Hermite-Laurent interpolation ? (y/n) ");
    Ask_Yes_or_No(lau);
    put("Test special or singular case ? (y/n) ");
    Ask_Yes_or_No(spc);
    new_line;
    if spc = 'y' then
      if lau = 'y'
       then DoblDobl_Singular_Hermite_Test;
       else DoblDobl_Special_Hermite_Test(deg);
      end if;
     else DoblDobl_Random_Hermite_Test(deg,dim,lau='y');
    end if;
  end DoblDobl_Hermite_Test;

  procedure QuadDobl_Hermite_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   To test the Hermite interpolation, we have a special case
  --   and the general case, tested on a randomly generated problem.
  --   Hermite interpolation is applied in quad double precision.

    lau,spc : character;

  begin
    put("Apply Hermite-Laurent interpolation ? (y/n) ");
    Ask_Yes_or_No(lau);
    put("Test special or singular case ? (y/n) ");
    Ask_Yes_or_No(spc);
    new_line;
    if spc = 'y' then
      if lau = 'y'
       then QuadDobl_Singular_Hermite_Test;
       else QuadDobl_Special_Hermite_Test(deg);
      end if;
     else QuadDobl_Random_Hermite_Test(deg,dim,lau='y');
    end if;
  end QuadDobl_Hermite_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree and the dimension of the series.
  --   The degree is the highest exponent in the power series.
  --   The dimension is the number of variables in the series.

    deg,dim : integer32 := 0;
    hrm,ans : character;

  begin
    new_line;
    put_line("Testing solving linear systems with interpolation.");
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    put("Hermite interpolation ? (y/n) ");
    Ask_Yes_or_No(hrm);
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    if hrm = 'y' then
      case ans is
        when '0' => Standard_Hermite_Test(deg,dim);
        when '1' => DoblDobl_Hermite_Test(deg,dim);
        when '2' => QuadDobl_Hermite_Test(deg,dim);
        when others => null;
      end case;
    else
      new_line;
      case ans is
        when '0' => Standard_Test(deg,dim);
        when '1' => DoblDobl_Test(deg,dim);
        when '2' => QuadDobl_Test(deg,dim);
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_cseritp;
