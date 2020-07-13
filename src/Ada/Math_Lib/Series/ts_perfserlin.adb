with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
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

  procedure Standard_Indexed_Test
              ( needrco : in boolean; n,d : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                info : out integer32; rcond : out double_float ) is

  -- DESCRIPTION :
  --   An indexed solver starts with the head and the computes
  --   the solution term after term.
  --   This procedure tests the Solve_Head and Solve_Tail procedures,
  --   solving the linear system in a staggered manner.

  -- REQUIRED : A'range = 0..d = b'range = x'range,
  --   and for all k in 0..d: b(k)'range = 1..n,
  --   and A(k)'range(1) = A(k)'range(2) = 1..n.

  -- ON ENTRY :
  --   needrco  flag to indicate if condition estimate is needed;
  --   n        dimension of the linear system;
  --   d        degree of the power series;
  --   A        matrix coefficients of the system;
  --   b        right hand side vector of the system.

  -- ON RETURN :
  --   b        vector coefficients of the solution series;
  --   info     pivoting information of lufac if no needrco;
  --   rcond    estimated for inverse condition number if needrco.

    rcols,icols,rb,ib : Standard_Floating_VecVecs.Link_to_VecVec;
    ry,iy : Standard_Floating_Vectors.Link_to_Vector;
    rv,iv : Standard_Floating_VecVecVecs.Link_to_VecVecVec;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    wrkidx,wrkdeg : integer32;

    use Standard_Inlined_Linearization;

  begin
    Standard_Floating_VecVecVecs.Allocate(rv,1,d,1,n,1,n);
    Standard_Floating_VecVecVecs.Allocate(iv,1,d,1,n,1,n);
    Standard_Matrix_Splitters.Split_Rows(A,rv,iv);
    rcols := Standard_Vector_Splitters.Allocate(n,n,1,1);
    icols := Standard_Vector_Splitters.Allocate(n,n,1,1);
    rb := Standard_Vector_Splitters.Allocate(d,n,0,1);
    ib := Standard_Vector_Splitters.Allocate(d,n,0,1);
    ry := new Standard_Floating_Vectors.Vector'(1..n => 0.0);
    iy := new Standard_Floating_Vectors.Vector'(1..n => 0.0);
    Standard_Matrix_Splitters.Complex_Parts(A(0).all,rcols,icols);
    Standard_Vector_Splitters.Complex_Parts(b,rb,ib);
    if needrco then
      Inlined_Solve_Head_by_lufco(n,rcols,icols,rb(0),ib(0),ipvt,rcond,ry,iy);
    else
      Inlined_Solve_Head_by_lufac(n,rcols,icols,rb(0),ib(0),ipvt,info);
    end if;
    Standard_Vector_Splitters.Complex_Merge(rb(0),ib(0),b(0));
    wrkdeg := 1; wrkidx := 1;
    loop
      Inlined_Solve_Tail_by_lusolve
        (wrkdeg,n,rcols,icols,rv,iv,rb,ib,ipvt,ry,iy,wrkidx);
      for k in wrkidx..wrkdeg loop
        Standard_Vector_Splitters.Complex_Merge(rb(k),ib(k),b(k));
      end loop;
      exit when (wrkdeg = d);
      wrkidx := wrkdeg+1;
      wrkdeg := 2*wrkdeg;
      if wrkdeg > d
       then wrkdeg := d;
      end if;
    end loop;
  end Standard_Indexed_Test;

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
    vm3 : Standard_Complex_VecMats.VecMat(0..As.deg);
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
    bcf3 : Standard_Complex_VecVecs.VecVec(0..bs.deg);
    info : integer32;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    wrk : constant Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..n);
    rcond : double_float;
    ans : character;
    lurcond,indexed : boolean := false;

  begin
    new_line;
    put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
    lurcond := (ans = 'y');
    new_line;
    put("Test indexed solver ? (y/n) "); Ask_Yes_or_No(ans);
    indexed := (ans = 'y');
    if indexed then
      vm3 := Series_Coefficient_Vectors.Standard_Series_Coefficients(As);
      bcf3 := Series_Coefficient_Vectors.Standard_Series_Coefficients(bs);
    end if;
    new_line;
    if lurcond then
      Standard_Series_Matrix_Solvers.Solve_by_lufco(vm1,bcf1,ipvt,rcond,wrk);
      put("rcond : "); put(rcond); new_line;
    else
      Standard_Series_Matrix_Solvers.Solve_by_lufac(vm1,bcf1,ipvt,info,wrk);
      put("info : "); put(info,1); new_line;
    end if;
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
    new_line;
    if lurcond then
      Standard_Inlined_Linearization.Inlined_Solve_by_lufco
        (vm2,bcf2,ipvt,rcond);
      put("rcond : "); put(rcond); new_line;
    else
      Standard_Inlined_Linearization.Inlined_Solve_by_lufac
        (vm2,bcf2,ipvt,info);
      put("info : "); put(info,1); new_line;
    end if;
    if indexed then
      Standard_Indexed_Test(lurcond,n,d,vm3,bcf3,info,rcond);
      if lurcond
       then put("rcond : "); put(rcond); new_line;
       else put("info : "); put(info,1); new_line;
      end if;
    end if;
    put_line("The generated leading vector series of the solution :");
    put_line(xs.cff(0));
    put_line("The recomputed leading vector series of the solution :");
    put_line(bcf2(0));
    if indexed then
      put_line("The indexed recomputed head of the solution :");
      put_line(bcf3(0));
    end if;
    for k in 1..bs.deg loop
      put("The term "); put(k,1); put_line(" in the solution :");
      put_line(xs.cff(k));
      put("Recomputed term "); put(k,1); put_line(" in the solution :");
      put_line(bcf2(k));
      put("Indexed recomputed term "); put(k,1);
      put_line(" in the solution :"); put_line(bcf3(k));
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
    rcond : double_float;
    ans : character;

  begin
    new_line;
    put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
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
      if ans = 'y' then
        Standard_Series_Matrix_Solvers.Solve_by_lufco(vm1,bcf1,ipvt,rcond,wrk);
      else
        Standard_Series_Matrix_Solvers.Solve_by_lufac(vm1,bcf1,ipvt,info,wrk);
      end if;
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
      if ans = 'y' then
        Standard_Inlined_Linearization.Inlined_Solve_by_lufco
          (vm2,bcf2,ipvt,rcond);
      else
        Standard_Inlined_Linearization.Inlined_Solve_by_lufac
          (vm2,bcf2,ipvt,info);
      end if;
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"wrapped inlined linearization");
    Standard_Floating_VecVecVecs.Allocate(rv,1,d,1,n,1,n);
    Standard_Floating_VecVecVecs.Allocate(iv,1,d,1,n,1,n);
    Standard_Matrix_Splitters.Split_Rows(vm2,rv,iv);
    rcols := Standard_Vector_Splitters.Allocate(n,n,1,1);
    icols := Standard_Vector_Splitters.Allocate(n,n,1,1);
    rb := Standard_Vector_Splitters.Allocate(d,n,0,1);
    ib := Standard_Vector_Splitters.Allocate(d,n,0,1);
    ry := new Standard_Floating_Vectors.Vector'(1..n => 0.0);
    iy := new Standard_Floating_Vectors.Vector'(1..n => 0.0);
    tstart(timer);
    for k in 1..f loop
      Standard_Matrix_Splitters.Complex_Parts(A0,rcols,icols);
      Standard_Vector_Splitters.Complex_Parts(bcf2,rb,ib);
      if ans = 'y' then
        Standard_Inlined_Linearization.Inlined_Solve_by_lufco
          (n,rcols,icols,rv,iv,rb,ib,ipvt,rcond,ry,iy);
      else
        Standard_Inlined_Linearization.Inlined_Solve_by_lufac
          (n,rcols,icols,rv,iv,rb,ib,ipvt,info,ry,iy);
      end if;
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
