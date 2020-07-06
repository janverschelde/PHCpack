with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with Standard_Vector_Splitters;
with Standard_Matrix_Splitters;
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
  --   double precision complex coefficients into a matrix series,
  --   and then solves the linear system of power series twice,
  --   once with complex, and once with inlined LU factorization.

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

  procedure Standard_Timing_Test ( n,d,f : in integer32 ) is

  -- DESCRIPTION :
  --   Generates an n-by-n matrix of series of degree d,
  --   with complex coefficients in standard double precision.
  --   Converts an n-by-n matrix of series of degree d with standard
  --   double precision complex coefficients into a matrix series.
  --   There are three operations that are timed:
  --   1) the complex linearization to solve a linear series system;
  --   2) wrapped inlined linearization, with allocation/deallocation;
  --   3) inlined linearization on allocated floating-point vectors.

    use Standard_Complex_Series_Matrices;

    timer : Timing_Widget;
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
    sb : constant Standard_Complex_Series_Vectors.Vector(1..n) := sA*sx;
    bs : constant Standard_Complex_Vector_Series.Vector(d)
       := Standard_Complex_Vector_Series.Create(sb);
    bcf1 : Standard_Complex_VecVecs.VecVec(0..bs.deg)
         := Series_Coefficient_Vectors.Standard_Series_Coefficients(bs);
    bcf2 : Standard_Complex_VecVecs.VecVec(0..bs.deg)
         := Series_Coefficient_Vectors.Standard_Series_Coefficients(bs);
    info : integer32;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    wrk : constant Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..n);
    lead : Standard_Complex_Matrices.Link_to_Matrix;
    A0 : Standard_Complex_Matrices.Matrix(1..n,1..n);
    rcols,icols,rb,ib : Standard_Floating_VecVecs.Link_to_VecVec;
    ry,iy : Standard_Floating_Vectors.Link_to_Vector;
    rv,iv : Standard_Floating_VecVecVecs.Link_to_VecVecVec;

  begin
    lead := vm1(0);
    for i in 1..n loop        -- copy lead matrix into the backup A0
      for j in 1..n loop
         A0(i,j) := lead(i,j);
      end loop;
    end loop;
    tstart(timer);
    for k in 1..f loop
      for i in 1..n loop      -- restore the lead matrix from backup A0
        for j in 1..n loop
          lead(i,j) := A0(i,j);
        end loop;
      end loop;
      bcf1 := Series_Coefficient_Vectors.Standard_Series_Coefficients(bs);
      Standard_Series_Matrix_Solvers.Solve_by_lufac(vm1,bcf1,ipvt,info,wrk);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex linearization");
    lead := vm2(0);
    for i in 1..n loop       -- copy lead matrix into the backup A0
      for j in 1..n loop
         A0(i,j) := lead(i,j);
      end loop;
    end loop;
    tstart(timer);
    for k in 1..f loop
      for i in 1..n loop      -- restore the lead matrix from backup A0
        for j in 1..n loop
          lead(i,j) := A0(i,j);
        end loop;
      end loop;
      bcf2 := Series_Coefficient_Vectors.Standard_Series_Coefficients(bs);
      Standard_Inlined_Linearization.Inlined_Solve_by_lufac
        (vm2,bcf2,ipvt,info);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"wrapped inlined linearization");
    Standard_Floating_VecVecVecs.Allocate(rv,1,d,1,n,1,n);
    Standard_Floating_VecVecVecs.Allocate(iv,1,d,1,n,1,n);
    Standard_Matrix_Splitters.Split_Rows(vm2,rv,iv);
    rcols := Standard_Vector_Splitters.Allocate(n,n,1,1);
    icols := Standard_Vector_Splitters.Allocate(n,n,1,1);
    rb := Standard_Vector_Splitters.Allocate(n,n,0,1);
    ib := Standard_Vector_Splitters.Allocate(n,n,0,1);
    ry := new Standard_Floating_Vectors.Vector'(1..n => 0.0);
    iy := new Standard_Floating_Vectors.Vector'(1..n => 0.0);
    tstart(timer);
    for k in 1..f loop
      Standard_Matrix_Splitters.Complex_Parts(A0,rcols,icols);
      bcf2 := Series_Coefficient_Vectors.Standard_Series_Coefficients(bs);
      Standard_Inlined_Linearization.Inlined_Solve_by_lufac
        (n,bcf2,ipvt,info,rv,iv,rcols,icols,rb,ib,ry,iy);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"real inlined linearization");
    Standard_Floating_VecVecs.Deep_Clear(rcols);
    Standard_Floating_VecVecs.Deep_Clear(icols);
    Standard_Floating_VecVecs.Deep_Clear(rb);
    Standard_Floating_VecVecs.Deep_Clear(ib);
    Standard_Floating_Vectors.Clear(ry);
    Standard_Floating_Vectors.Clear(iy);
    Standard_Floating_VecVecVecs.Clear(rv);
    Standard_Floating_VecVecVecs.Clear(iv);
  end Standard_Timing_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension of the linear system
  --   and the degrees of the series in the system.

    dim,deg,frq : integer32 := 0;

  begin
    new_line;
    put_line("Performance of solving linear systems of power series ...");
    put("  Give the dimension of the system : "); get(dim);
    put("  Give the degree of the series : "); get(deg);
    put("  Give the frequency (0 for interactive testing) : "); get(frq);
    if frq = 0
     then Standard_Test(dim,deg);
     else Standard_Timing_Test(dim,deg,frq);
    end if;
  end Main;

begin
  Main;
end ts_perfserlin;
