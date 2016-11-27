with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Random_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Standard_Complex_VecMats;
with Standard_Dense_Series;
with Standard_Dense_Vector_Series;
with Standard_Dense_Vector_Series_io;    use Standard_Dense_Vector_Series_io;
with Standard_Random_Series;
with Standard_Dense_Matrix_Series;
with Standard_Dense_Matrix_Series_io;    use Standard_Dense_Matrix_Series_io;
with Standard_Interpolating_Series;

procedure ts_seritp is

-- DESCRIPTION :
--   Development of solving linear systems of series with interpolation.

  function Random_Matrix_Series
             ( deg,dim : integer32 )
             return Standard_Dense_Matrix_Series.Matrix is

  -- DESCRIPTION :
  --   Returns a random matrix series, with coefficients zero or one,
  --   of the given degree deg and dimension dim.

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
  end Random_Matrix_Series;

  function Multiply
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
  end Multiply;

  function Interpolate
             ( mat : in Standard_Dense_Matrix_Series.Matrix;
               rhs : in Standard_Dense_Vector_Series.Vector;
               verbose : in boolean := true )
             return Standard_Dense_Vector_Series.Vector is

  -- DESCRIPTION :
  --   Samples the matrix and vector series at random points
  --   and solves the linear systems defined by the evaluated
  --   coefficient matrices and right hand side vectors.

    res : Standard_Dense_Vector_Series.Vector;
    dim : constant integer32 := rhs.cff(0)'last;
    t : constant Standard_Complex_Vectors.Vector(0..mat.deg)
      := Standard_Random_Vectors.Random_Vector(0,mat.deg);
    m : Standard_Complex_VecMats.VecMat(t'range)
      := Standard_Interpolating_Series.Sample(mat,t);
    v : Standard_Complex_VecVecs.VecVec(t'range)
      := Standard_Interpolating_Series.Sample(rhs,t);
    x : Standard_Complex_VecVecs.VecVec(t'range);
    r : Standard_Floating_Vectors.Vector(t'range);
    xt : Standard_Complex_VecVecs.VecVec(1..dim);
    vdm : Standard_Complex_Matrices.Matrix(1..mat.deg+1,1..mat.deg+1);
    cff : Standard_Complex_VecVecs.VecVec(1..dim);

  begin
    if verbose then
      put_line("The sample points :");
      put_line(t);
    end if;
    x := Standard_Interpolating_Series.Solve_Linear_Systems(m,v);
    r := Standard_Interpolating_Series.Residuals(m,v,x);
    if verbose then
      put_line("The solutions to the interpolated linear systems :");
      put_line(x);
      put_line("The residuals of the solved sampled linear systems :");
      put_line(r);
    end if;
    xt := Standard_Interpolating_Series.Transpose(x);
    if verbose then
      put_line("The transposed solution vectors :");
      put_line(xt);
    end if;
    vdm := Standard_Interpolating_Series.Vandermonde_Matrix(t);
    cff := Standard_Interpolating_Series.Solve_Interpolation_Systems(vdm,xt);
    if verbose then
      put_line("The coefficients computed via interpolation :");
      put_line(cff);
    end if;
    res := Standard_Interpolating_Series.Construct(cff);
    return res;
  end Interpolate;

  procedure Standard_Test ( deg,dim : integer32 ) is

  -- DESCRIPTION :
  --   Generates a linear system of the given degree and dimension.

    mat : constant Standard_Dense_Matrix_Series.Matrix
        := Random_Matrix_Series(deg,dim);
    sol : constant Standard_Dense_Vector_Series.Vector
        := Standard_Random_Series.Random_Vector_Series(1,dim,deg);
    rhs : constant Standard_Dense_Vector_Series.Vector
        := Multiply(mat,sol);
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
      cff := Interpolate(mat,rhs);
      put_line("The computed solution :"); put(cff);
      put_line("The constructed solution :"); put(sol);
    end if;
  end Standard_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree and the dimension of the series.
  --   The degree is the highest exponent in the power series.
  --   The dimension is the number of variables in the series.

    deg,dim : integer32 := 0;

  begin
    new_line;
    put_line("Testing solving linear systems with interpolation.");
    put("Give the dimension : "); get(dim);
    put("Give the degree : "); get(deg);
    Standard_Test(deg,dim);
  end Main;

begin
  Main;
end ts_seritp;
