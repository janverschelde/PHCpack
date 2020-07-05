with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with Standard_Complex_Vector_Series;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Matrix_Series;
with Standard_Complex_Series_Matrices;
with Standard_Random_Series_Vectors;
with Standard_Random_Series_Matrices;
with Standard_Series_Matrix_Solvers;
with Series_Coefficient_Vectors;
with Standard_Inlined_Linearization;

procedure ts_perfserlin is

-- DESCRIPTION :
--   Test the performance on solving linear systems of power series with
--   linearization, flat data structures, and inlined linear solvers.

  procedure Standard_Test ( n,d : in integer32 ) is

  -- DESCRIPTION :
  --   Generates an n-by-n matrix of series of degree d,
  --   with complex coefficients in standard double precision.
  --   Converts an n-by-n matrix of series of degree d with standard
  --   double precision complex coefficients into a matrix series.

    use Standard_Complex_Series_Matrices;

    sA : constant Standard_Complex_Series_Matrices.Matrix(1..n,1..n)
       := Standard_Random_Series_Matrices.Random_Series_Matrix(1,n,1,n,d);
    As : constant Standard_Complex_Matrix_Series.Matrix 
       := Standard_Complex_Matrix_Series.Create(sA); 
    vm1 : constant Standard_Complex_VecMats.VecMat(0..As.deg)
        := Series_Coefficient_Vectors.Standard_Series_Coefficients(As);
    vm2 : constant Standard_Complex_VecMats.VecMat(0..As.deg)
        := Series_Coefficient_Vectors.Standard_Series_Coefficients(As);
    sx : constant Standard_Complex_Series_Vectors.Vector(1..n)
       := Standard_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    xs : constant Standard_Complex_Vector_Series.Vector(d)
       := Standard_Complex_Vector_Series.Create(sx);
    sb : constant Standard_Complex_Series_Vectors.Vector(1..n) := sA*sx;
    bs : constant Standard_Complex_Vector_Series.Vector(d)
       := Standard_Complex_Vector_Series.Create(sb);
    bcf1 : constant Standard_Complex_VecVecs.VecVec(0..bs.deg)
         := Series_Coefficient_Vectors.Standard_Series_Coefficients(bs);
    bcf2 : constant Standard_Complex_VecVecs.VecVec(0..bs.deg)
         := Series_Coefficient_Vectors.Standard_Series_Coefficients(bs);
    info : integer32;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    wrk : constant Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..n);

  begin
    Standard_Series_Matrix_Solvers.Solve_by_lufac(vm1,bcf1,ipvt,info,wrk);
    new_line;
    put_line("The generated leading vector series of the solution :");
    put_line(xs.cff(0));
    put_line("The computed leading vector series of the solution :");
    put_line(bcf1(0));
    for k in 1..bs.deg loop
      put("The generated term "); put(k,1);
      put_line(" of the vector series of the solution :"); put_line(xs.cff(k));
      put("The computed term "); put(k,1);
      put_line(" of the vector series of the solution :"); put_line(bcf1(k));
    end loop;
    Standard_Inlined_Linearization.Inlined_Solve_by_lufac(vm2,bcf2,ipvt,info);
    new_line;
    put_line("The generated leading vector series of the solution :");
    put_line(xs.cff(0));
    put_line("The recomputed leading vector series of the solution :");
    put_line(bcf2(0));
    for k in 1..bs.deg loop
      put("The term "); put(k,1); put_line(" in the solution :");
      put_line(xs.cff(k));
      put("Recomputed term "); put(k,1); put_line(" in the solution :");
      put_line(bcf2(k));
    end loop;
  end Standard_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension of the linear system
  --   and the degrees of the series in the system.

    dim,deg : integer32 := 0;

  begin
    new_line;
    put_line("Performance of solving linear systems of power series ...");
    put("  Give the dimension of the system : "); get(dim);
    put("  Give the degree of the series : "); get(deg);
    Standard_Test(dim,deg);
  end Main;

begin
  Main;
end ts_perfserlin;
