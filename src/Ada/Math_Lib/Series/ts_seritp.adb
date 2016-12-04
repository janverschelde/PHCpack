with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Random_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Matrices;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Standard_Complex_VecMats;
with Standard_Dense_Series;
with Standard_Dense_Vector_Series;
with Standard_Dense_Vector_Series_io;    use Standard_Dense_Vector_Series_io;
with Standard_Random_Series;
with DoblDobl_Random_Series;
with QuadDobl_Random_Series;
with Standard_Dense_Matrix_Series;
with Standard_Dense_Matrix_Series_io;    use Standard_Dense_Matrix_Series_io;
with DoblDobl_Dense_Series;
with DoblDobl_Dense_Vector_Series;
with DoblDobl_Dense_Vector_Series_io;    use DoblDobl_Dense_Vector_Series_io;
with DoblDobl_Dense_Matrix_Series;
with DoblDobl_Dense_Matrix_Series_io;    use DoblDobl_Dense_Matrix_Series_io;
with QuadDobl_Dense_Series;
with QuadDobl_Dense_Vector_Series;
with QuadDobl_Dense_Vector_Series_io;    use QuadDobl_Dense_Vector_Series_io;
with QuadDobl_Dense_Matrix_Series;
with QuadDobl_Dense_Matrix_Series_io;    use QuadDobl_Dense_Matrix_Series_io;
with Standard_Interpolating_Series;
with DoblDobl_Interpolating_Series;
with QuadDobl_Interpolating_Series;

procedure ts_seritp is

-- DESCRIPTION :
--   Development of solving linear systems of series with interpolation.

  function Standard_Random_Matrix_Series
             ( deg,dim : integer32 )
             return Standard_Dense_Matrix_Series.Matrix is

  -- DESCRIPTION :
  --   Returns a random matrix series, with coefficients zero or one,
  --   of the given degree deg and dimension dim.
  --   The coefficients are stored in standard double precision.

    res : Standard_Dense_Matrix_Series.Matrix;

  begin
    if deg > Standard_Dense_Series.max_deg
     then res.deg := Standard_Dense_Series.max_deg;
     else res.deg := deg;
    end if;
    for d in 0..res.deg loop
      declare
        mat : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
        rnd : integer32;
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            rnd := Standard_Random_Numbers.Random(0,1);
            mat(i,j) := Standard_Complex_Numbers.Create(double_float(rnd));
          end loop;
        end loop;
        res.cff(d) := new Standard_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end Standard_Random_Matrix_Series;

  function DoblDobl_Random_Matrix_Series
             ( deg,dim : integer32 )
             return DoblDobl_Dense_Matrix_Series.Matrix is

  -- DESCRIPTION :
  --   Returns a random matrix series, with coefficients zero or one,
  --   of the given degree deg and dimension dim.

    res : DoblDobl_Dense_Matrix_Series.Matrix;

  begin
    if deg > DoblDobl_Dense_Series.max_deg
     then res.deg := DoblDobl_Dense_Series.max_deg;
     else res.deg := deg;
    end if;
    for d in 0..res.deg loop
      declare
        mat : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
        rnd : integer32;
        ddr : double_double;
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            rnd := Standard_Random_Numbers.Random(0,1);
            ddr := Double_Double_Numbers.Create(rnd);
            mat(i,j) := DoblDobl_Complex_Numbers.Create(ddr);
          end loop;
        end loop;
        res.cff(d) := new DoblDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end DoblDobl_Random_Matrix_Series;

  function QuadDobl_Random_Matrix_Series
             ( deg,dim : integer32 )
             return QuadDobl_Dense_Matrix_Series.Matrix is

  -- DESCRIPTION :
  --   Returns a random matrix series, with coefficients zero or one,
  --   of the given degree deg and dimension dim.

    res : QuadDobl_Dense_Matrix_Series.Matrix;

  begin
    if deg > QuadDobl_Dense_Series.max_deg
     then res.deg := QuadDobl_Dense_Series.max_deg;
     else res.deg := deg;
    end if;
    for d in 0..res.deg loop
      declare
        mat : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
        rnd : integer32;
        qdr : quad_double;
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            rnd := Standard_Random_Numbers.Random(0,1);
            qdr := Quad_Double_Numbers.Create(rnd);
            mat(i,j) := QuadDobl_Complex_Numbers.Create(qdr);
          end loop;
        end loop;
        res.cff(d) := new QuadDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end QuadDobl_Random_Matrix_Series;

  function Standard_Multiply
             ( mat : Standard_Dense_Matrix_Series.Matrix;
               vec : Standard_Dense_Vector_Series.Vector )
             return Standard_Dense_Vector_Series.Vector is

  -- DESCRIPTION :
  --   Multiplies the matrix series with the vector series,
  --   also convoluting the additional terms of higher degrees,
  --   up to the maximum degree.

  -- REQUIRED : mat.deg = vec.deg.

    res : Standard_Dense_Vector_Series.Vector;
    deg : constant integer32 := mat.deg + vec.deg;
    dim : constant integer32 := vec.cff(0)'last;
    xdg : integer32;

    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;

  begin
    if deg > Standard_Dense_Series.max_deg
     then res.deg := Standard_Dense_Series.max_deg;
     else res.deg := deg;
    end if;
    for k in 0..mat.deg loop
      declare
        acc : Standard_Complex_Vectors.Vector(1..dim)
            := mat.cff(0).all*vec.cff(k).all;
      begin
        for i in 1..k loop
          acc := acc + mat.cff(i).all*vec.cff(k-i).all;
        end loop;
        res.cff(k) := new Standard_Complex_Vectors.Vector'(acc);
      end;
    end loop;
    xdg := mat.deg+1;
    for k in 1..mat.deg loop -- computed extended degree terms
      declare
        acc : Standard_Complex_Vectors.Vector(1..dim)
            := mat.cff(k).all*vec.cff(xdg-k).all;
      begin
        for i in (k+1)..mat.deg loop
          acc := acc + mat.cff(i).all*vec.cff(xdg-i).all;
        end loop;
        res.cff(xdg) := new Standard_Complex_Vectors.Vector'(acc);
      end;
      xdg := xdg + 1;
    end loop;
    return res;
  end Standard_Multiply;

  function DoblDobl_Multiply
             ( mat : DoblDobl_Dense_Matrix_Series.Matrix;
               vec : DoblDobl_Dense_Vector_Series.Vector )
             return DoblDobl_Dense_Vector_Series.Vector is

  -- DESCRIPTION :
  --   Multiplies the matrix series with the vector series,
  --   also convoluting the additional terms of higher degrees,
  --   up to the maximum degree.

  -- REQUIRED : mat.deg = vec.deg.

    res : DoblDobl_Dense_Vector_Series.Vector;
    deg : constant integer32 := mat.deg + vec.deg;
    dim : constant integer32 := vec.cff(0)'last;
    xdg : integer32;

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;

  begin
    if deg > DoblDobl_Dense_Series.max_deg
     then res.deg := DoblDobl_Dense_Series.max_deg;
     else res.deg := deg;
    end if;
    for k in 0..mat.deg loop
      declare
        acc : DoblDobl_Complex_Vectors.Vector(1..dim)
            := mat.cff(0).all*vec.cff(k).all;
      begin
        for i in 1..k loop
          acc := acc + mat.cff(i).all*vec.cff(k-i).all;
        end loop;
        res.cff(k) := new DoblDobl_Complex_Vectors.Vector'(acc);
      end;
    end loop;
    xdg := mat.deg+1;
    for k in 1..mat.deg loop -- computed extended degree terms
      declare
        acc : DoblDobl_Complex_Vectors.Vector(1..dim)
            := mat.cff(k).all*vec.cff(xdg-k).all;
      begin
        for i in (k+1)..mat.deg loop
          acc := acc + mat.cff(i).all*vec.cff(xdg-i).all;
        end loop;
        res.cff(xdg) := new DoblDobl_Complex_Vectors.Vector'(acc);
      end;
      xdg := xdg + 1;
    end loop;
    return res;
  end DoblDobl_Multiply;

  function QuadDobl_Multiply
             ( mat : QuadDobl_Dense_Matrix_Series.Matrix;
               vec : QuadDobl_Dense_Vector_Series.Vector )
             return QuadDobl_Dense_Vector_Series.Vector is

  -- DESCRIPTION :
  --   Multiplies the matrix series with the vector series,
  --   also convoluting the additional terms of higher degrees,
  --   up to the maximum degree.

  -- REQUIRED : mat.deg = vec.deg.

    res : QuadDobl_Dense_Vector_Series.Vector;
    deg : constant integer32 := mat.deg + vec.deg;
    dim : constant integer32 := vec.cff(0)'last;
    xdg : integer32;

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Matrices;

  begin
    if deg > QuadDobl_Dense_Series.max_deg
     then res.deg := QuadDobl_Dense_Series.max_deg;
     else res.deg := deg;
    end if;
    for k in 0..mat.deg loop
      declare
        acc : QuadDobl_Complex_Vectors.Vector(1..dim)
            := mat.cff(0).all*vec.cff(k).all;
      begin
        for i in 1..k loop
          acc := acc + mat.cff(i).all*vec.cff(k-i).all;
        end loop;
        res.cff(k) := new QuadDobl_Complex_Vectors.Vector'(acc);
      end;
    end loop;
    xdg := mat.deg+1;
    for k in 1..mat.deg loop -- computed extended degree terms
      declare
        acc : QuadDobl_Complex_Vectors.Vector(1..dim)
            := mat.cff(k).all*vec.cff(xdg-k).all;
      begin
        for i in (k+1)..mat.deg loop
          acc := acc + mat.cff(i).all*vec.cff(xdg-i).all;
        end loop;
        res.cff(xdg) := new QuadDobl_Complex_Vectors.Vector'(acc);
      end;
      xdg := xdg + 1;
    end loop;
    return res;
  end QuadDobl_Multiply;

  procedure Standard_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension.
  --   Solves the problem in standard double precision.

    mat : constant Standard_Dense_Matrix_Series.Matrix
        := Standard_Random_Matrix_Series(deg,dim);
    sol : constant Standard_Dense_Vector_Series.Vector
        := Standard_Random_Series.Random_Vector_Series(1,dim,deg);
    rhs : constant Standard_Dense_Vector_Series.Vector
        := Standard_Multiply(mat,sol);
    rnk : integer32;
    cff : Standard_Dense_Vector_Series.Vector;

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
    rnk := Standard_Interpolating_Series.Full_Rank(mat);
    if rnk = -1 then
      put_line("The matrix series does not have full rank.");
    else
      put_line("The matrix series has full rank.");
      put("The smallest degree for full rank : "); put(rnk,1); new_line;
      cff := Standard_Interpolating_Series.Interpolate(mat,rhs);
      put_line("The computed solution :"); put(cff);
      put_line("The constructed solution :"); put(sol);
    end if;
  end Standard_Test;

  procedure DoblDobl_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension.
  --   Solves the problem in standard double precision.

    mat : constant DoblDobl_Dense_Matrix_Series.Matrix
        := DoblDobl_Random_Matrix_Series(deg,dim);
    sol : constant DoblDobl_Dense_Vector_Series.Vector
        := DoblDobl_Random_Series.Random_Vector_Series(1,dim,deg);
    rhs : constant DoblDobl_Dense_Vector_Series.Vector
        := DoblDobl_Multiply(mat,sol);
    rnk : integer32;
    cff : DoblDobl_Dense_Vector_Series.Vector;

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
    rnk := DoblDobl_Interpolating_Series.Full_Rank(mat);
    if rnk = -1 then
      put_line("The matrix series does not have full rank.");
    else
      put_line("The matrix series has full rank.");
      put("The smallest degree for full rank : "); put(rnk,1); new_line;
      cff := DoblDobl_Interpolating_Series.Interpolate(mat,rhs);
      put_line("The computed solution :"); put(cff);
      put_line("The constructed solution :"); put(sol);
    end if;
  end DoblDobl_Test;

  procedure QuadDobl_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension.
  --   Solves the problem in standard double precision.

    mat : constant QuadDobl_Dense_Matrix_Series.Matrix
        := QuadDobl_Random_Matrix_Series(deg,dim);
    sol : constant QuadDobl_Dense_Vector_Series.Vector
        := QuadDobl_Random_Series.Random_Vector_Series(1,dim,deg);
    rhs : constant QuadDobl_Dense_Vector_Series.Vector
        := QuadDobl_Multiply(mat,sol);
    rnk : integer32;
    cff : QuadDobl_Dense_Vector_Series.Vector;

  begin
    put_line("The generated matrix series :"); put(mat);
    put_line("The generated solution :"); put(sol);
    put_line("The multiplied right hand side vector : "); put(rhs);
    rnk := QuadDobl_Interpolating_Series.Full_Rank(mat);
    if rnk = -1 then
      put_line("The matrix series does not have full rank.");
    else
      put_line("The matrix series has full rank.");
      put("The smallest degree for full rank : "); put(rnk,1); new_line;
      cff := QuadDobl_Interpolating_Series.Interpolate(mat,rhs);
      put_line("The computed solution :"); put(cff);
      put_line("The constructed solution :"); put(sol);
    end if;
  end QuadDobl_Test;

  function Standard_Special_Matrix_Series
             ( deg : integer32 ) 
             return Standard_Dense_Matrix_Series.Matrix is

  -- DESCRIPTION :
  --   Returns a special matrix series of degree equal to deg,
  --   which has an infinite series as solution when the
  --   right hand side vector is (1, 0).
  --   The degree should be at least one.

    res : Standard_Dense_Matrix_Series.Matrix;
    wrk : Standard_Complex_Matrices.Matrix(1..2,1..2);

  begin
    res.deg := deg;
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
             return DoblDobl_Dense_Matrix_Series.Matrix is

  -- DESCRIPTION :
  --   Returns a special matrix series of degree equal to deg,
  --   which has an infinite series as solution when the
  --   right hand side vector is (1, 0).
  --   The degree should be at least one.

    res : DoblDobl_Dense_Matrix_Series.Matrix;
    wrk : DoblDobl_Complex_Matrices.Matrix(1..2,1..2);

  begin
    res.deg := deg;
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        wrk(i,j) := DoblDobl_Complex_Numbers.Create(integer32(0));
      end loop;
    end loop;
    for k in 0..deg loop
      res.cff(k) := new DoblDobl_Complex_Matrices.Matrix'(wrk);
    end loop;
    res.cff(0)(1,1) := DoblDobl_Complex_Numbers.Create(integer32(1));
    res.cff(0)(1,2) := DoblDobl_Complex_Numbers.Create(integer32(1));
    res.cff(0)(2,2) := DoblDobl_Complex_Numbers.Create(integer32(1));
    res.cff(1)(2,1) := DoblDobl_Complex_Numbers.Create(integer32(1));
    return res;
  end DoblDobl_Special_Matrix_Series;

  function QuadDobl_Special_Matrix_Series
             ( deg : integer32 ) 
             return QuadDobl_Dense_Matrix_Series.Matrix is

  -- DESCRIPTION :
  --   Returns a special matrix series of degree equal to deg,
  --   which has an infinite series as solution when the
  --   right hand side vector is (1, 0).
  --   The degree should be at least one.

    res : QuadDobl_Dense_Matrix_Series.Matrix;
    wrk : QuadDobl_Complex_Matrices.Matrix(1..2,1..2);

  begin
    res.deg := deg;
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        wrk(i,j) := QuadDobl_Complex_Numbers.Create(integer32(0));
      end loop;
    end loop;
    for k in 0..deg loop
      res.cff(k) := new QuadDobl_Complex_Matrices.Matrix'(wrk);
    end loop;
    res.cff(0)(1,1) := QuadDobl_Complex_Numbers.Create(integer32(1));
    res.cff(0)(1,2) := QuadDobl_Complex_Numbers.Create(integer32(1));
    res.cff(0)(2,2) := QuadDobl_Complex_Numbers.Create(integer32(1));
    res.cff(1)(2,1) := QuadDobl_Complex_Numbers.Create(integer32(1));
    return res;
  end QuadDobl_Special_Matrix_Series;

  function Standard_Special_Vector_Series
             ( deg : integer32 )
             return Standard_Dense_Vector_Series.Vector is

  -- DESCRIPTION :
  --   Returns a special right hand side series with degree equal to deg.

    res : Standard_Dense_Vector_Series.Vector;
    wrk : Standard_Complex_Vectors.Vector(1..2);

  begin
    for i in wrk'range loop
      wrk(i) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    res.deg := deg;
    for k in 0..deg loop
      res.cff(k) := new Standard_Complex_Vectors.Vector'(wrk);
    end loop;
    res.cff(0)(1) := Standard_Complex_Numbers.Create(1.0);
    return res;
  end Standard_Special_Vector_Series;

  function DoblDobl_Special_Vector_Series
             ( deg : integer32 )
             return DoblDobl_Dense_Vector_Series.Vector is

  -- DESCRIPTION :
  --   Returns a special right hand side series with degree equal to deg.

    res : DoblDobl_Dense_Vector_Series.Vector;
    wrk : DoblDobl_Complex_Vectors.Vector(1..2);

  begin
    for i in wrk'range loop
      wrk(i) := DoblDobl_Complex_Numbers.Create(integer32(0));
    end loop;
    res.deg := deg;
    for k in 0..deg loop
      res.cff(k) := new DoblDobl_Complex_Vectors.Vector'(wrk);
    end loop;
    res.cff(0)(1) := DoblDobl_Complex_Numbers.Create(integer32(1));
    return res;
  end DoblDobl_Special_Vector_Series;

  function QuadDobl_Special_Vector_Series
             ( deg : integer32 )
             return QuadDobl_Dense_Vector_Series.Vector is

  -- DESCRIPTION :
  --   Returns a special right hand side series with degree equal to deg.

    res : QuadDobl_Dense_Vector_Series.Vector;
    wrk : QuadDobl_Complex_Vectors.Vector(1..2);

  begin
    for i in wrk'range loop
      wrk(i) := QuadDobl_Complex_Numbers.Create(integer32(0));
    end loop;
    res.deg := deg;
    for k in 0..deg loop
      res.cff(k) := new QuadDobl_Complex_Vectors.Vector'(wrk);
    end loop;
    res.cff(0)(1) := QuadDobl_Complex_Numbers.Create(integer32(1));
    return res;
  end QuadDobl_Special_Vector_Series;

  procedure Standard_Hermite_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Applies Hermite interpolation on a special case,
  --   in standard double precision.
  --   The dimension will be assumed to be 2, no matter what dim is.

    mat : constant Standard_Dense_Matrix_Series.Matrix
        := Standard_Special_Matrix_Series(deg);
    rhs : constant Standard_Dense_Vector_Series.Vector
        := Standard_Special_Vector_Series(deg);
    sol : Standard_Dense_Vector_Series.Vector;
    t : Standard_Complex_Numbers.Complex_Number
      := Standard_Complex_Numbers.Create(0.01);

  begin
    put_line("The special matrix series : "); put(mat);
    put_line("The special vector series : "); put(rhs);
    sol := Standard_Interpolating_Series.Hermite_Interpolate(mat,rhs,t);
    put_line("The solution : "); put(sol);
  end Standard_Hermite_Test;

  procedure DoblDobl_Hermite_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Applies Hermite interpolation on a special case,
  --   in double double precision.
  --   The dimension will be assumed to be 2, no matter what dim is.

    mat : constant DoblDobl_Dense_Matrix_Series.Matrix
        := DoblDobl_Special_Matrix_Series(deg);
    rhs : constant DoblDobl_Dense_Vector_Series.Vector
        := DoblDobl_Special_Vector_Series(deg);
    sol : DoblDobl_Dense_Vector_Series.Vector;
    ddt : double_double := Double_Double_Numbers.Create(0.01);
    t : DoblDobl_Complex_Numbers.Complex_Number
      := DoblDobl_Complex_Numbers.Create(ddt);

  begin
    put_line("The special matrix series : "); put(mat);
    put_line("The special vector series : "); put(rhs);
    sol := DoblDobl_Interpolating_Series.Hermite_Interpolate(mat,rhs,t);
    put_line("The solution : "); put(sol);
  end DoblDobl_Hermite_Test;

  procedure QuadDobl_Hermite_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Applies Hermite interpolation on a special case,
  --   in quad double precision.
  --   The dimension will be assumed to be 2, no matter what dim is.

    mat : constant QuadDobl_Dense_Matrix_Series.Matrix
        := QuadDobl_Special_Matrix_Series(deg);
    rhs : constant QuadDobl_Dense_Vector_Series.Vector
        := QuadDobl_Special_Vector_Series(deg);
    sol : QuadDobl_Dense_Vector_Series.Vector;
    qdt : quad_double := Quad_Double_Numbers.Create(0.01);
    t : QuadDobl_Complex_Numbers.Complex_Number
      := QuadDobl_Complex_Numbers.Create(qdt);

  begin
    put_line("The special matrix series : "); put(mat);
    put_line("The special vector series : "); put(rhs);
    sol := QuadDobl_Interpolating_Series.Hermite_Interpolate(mat,rhs,t);
    put_line("The solution : "); put(sol);
  end QuadDobl_Hermite_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree and the dimension of the series.
  --   The degree is the highest exponent in the power series.
  --   The dimension is the number of variables in the series.

    deg,dim : integer32 := 0;
    ans,hrm : character;

  begin
    new_line;
    put_line("Testing solving linear systems with interpolation.");
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    put("Hermite interpolation ? (y/n) ");
    Ask_Yes_or_No(hrm);
    new_line;
    if hrm = 'y' then
      case ans is
        when '0' => Standard_Hermite_Test(deg,dim);
        when '1' => DoblDobl_Hermite_Test(deg,dim);
        when '2' => QuadDobl_Hermite_Test(deg,dim);
        when others => null;
      end case;
    else
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
end ts_seritp;
