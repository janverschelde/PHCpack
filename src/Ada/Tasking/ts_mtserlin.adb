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
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Standard_Complex_VecMats;
with Standard_Complex_VecMats_io;        use Standard_Complex_VecMats_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with DoblDobl_Complex_VecMats;
with DoblDobl_Complex_VecMats_io;        use DoblDobl_Complex_VecMats_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs_io;        use QuadDobl_Complex_VecVecs_io;
with QuadDobl_Complex_VecMats;
with QuadDobl_Complex_VecMats_io;        use QuadDobl_Complex_VecMats_io;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_Vectors_io; use Standard_Complex_Series_Vectors_io;
with Standard_Complex_Series_Matrices;
with Standard_Complex_Vector_Series;
with Standard_Complex_Vector_Series_io;  use Standard_Complex_Vector_Series_io;
with Standard_Complex_Matrix_Series;
with Standard_Complex_Matrix_Series_io;  use Standard_Complex_Matrix_Series_io;
with Standard_Random_Series_Vectors;
with Standard_Random_Series_Matrices;
with Standard_Series_Matrix_Solvers;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Vectors_io; use DoblDobl_Complex_Series_Vectors_io;
with DoblDobl_Complex_Series_Matrices;
with DoblDobl_Complex_Vector_Series;
with DoblDobl_Complex_Vector_Series_io;  use DoblDobl_Complex_Vector_Series_io;
with DoblDobl_Complex_Matrix_Series;
with DoblDobl_Complex_Matrix_Series_io;  use DoblDobl_Complex_Matrix_Series_io;
with DoblDobl_Random_Series_Vectors;
with DoblDobl_Random_Series_Matrices;
with DoblDobl_Series_Matrix_Solvers;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Vectors_io; use QuadDobl_Complex_Series_Vectors_io;
with QuadDobl_Complex_Series_Matrices;
with QuadDobl_Complex_Vector_Series;
with QuadDobl_Complex_Vector_Series_io;  use QuadDobl_Complex_Vector_Series_io;
with QuadDobl_Complex_Matrix_Series;
with QuadDobl_Complex_Matrix_Series_io;  use QuadDobl_Complex_Matrix_Series_io;
with QuadDobl_Random_Series_Vectors;
with QuadDobl_Random_Series_Matrices;
with QuadDobl_Series_Matrix_Solvers;
with Series_Coefficient_Vectors;
with Evaluation_Differentiation_Errors;  use Evaluation_Differentiation_Errors;
with Multitasked_Series_Linearization;   use Multitasked_Series_Linearization;

procedure ts_mtserlin is

-- DESCRIPTION :
--   Tests the linearization of solving linear systems of truncated series
--   with multitasking.

  function Error ( xs : Standard_Complex_Vector_Series.Vector;
                   bscff : Standard_Complex_VecVecs.VecVec;
                   output : boolean := true ) return double_float is

  -- DESCRIPTION :
  --   Returns the sum of all differences between the coefficients in xs
  --   and those in bscff.  If output, then the coefficient vectors are
  --   written to screen for both xs and bscff.

    err : double_float;

  begin
    if output then
      put_line("The generated leading vector series of the solution :");
      put_line(xs.cff(0));
      put_line("The computed leading vector series of the solution :");
      put_line(bscff(0));
    end if;
    err := Difference(xs.cff(0),bscff(0));
    for k in 1..xs.deg loop
      if output then
        put("The generated term "); put(k,1);
        put_line(" of the vector series of the solution :");
        put_line(xs.cff(k));
        put("The computed term "); put(k,1);
        put_line(" of the vector series of the solution :");
        put_line(bscff(k));
      end if;
      err := err + Difference(xs.cff(k),bscff(k));
    end loop;
    return err;
  end Error;

  function Error ( xs : DoblDobl_Complex_Vector_Series.Vector;
                   bscff : DoblDobl_Complex_VecVecs.VecVec;
                   output : boolean := true ) return double_double is

  -- DESCRIPTION :
  --   Returns the sum of all differences between the coefficients in xs
  --   and those in bscff.  If output, then the coefficient vectors are
  --   written to screen for both xs and bscff.

    err : double_double;

  begin
    if output then
      put_line("The generated leading vector series of the solution :");
      put_line(xs.cff(0));
      put_line("The computed leading vector series of the solution :");
      put_line(bscff(0));
    end if;
    err := Difference(xs.cff(0),bscff(0));
    for k in 1..xs.deg loop
      if output then
        put("The generated term "); put(k,1);
        put_line(" of the vector series of the solution :");
        put_line(xs.cff(k));
        put("The computed term "); put(k,1);
        put_line(" of the vector series of the solution :");
        put_line(bscff(k));
      end if;
      err := err + Difference(xs.cff(k),bscff(k));
    end loop;
    return err;
  end Error;

  function Error ( xs : QuadDobl_Complex_Vector_Series.Vector;
                   bscff : QuadDobl_Complex_VecVecs.VecVec;
                   output : boolean := true ) return quad_double is

  -- DESCRIPTION :
  --   Returns the sum of all differences between the coefficients in xs
  --   and those in bscff.  If output, then the coefficient vectors are
  --   written to screen for both xs and bscff.

    err : quad_double;

  begin
    if output then
      put_line("The generated leading vector series of the solution :");
      put_line(xs.cff(0));
      put_line("The computed leading vector series of the solution :");
      put_line(bscff(0));
    end if;
    err := Difference(xs.cff(0),bscff(0));
    for k in 1..xs.deg loop
      if output then
        put("The generated term "); put(k,1);
        put_line(" of the vector series of the solution :");
        put_line(xs.cff(k));
        put("The computed term "); put(k,1);
        put_line(" of the vector series of the solution :");
        put_line(bscff(k));
      end if;
      err := err + Difference(xs.cff(k),bscff(k));
    end loop;
    return err;
  end Error;

  procedure Standard_Test ( nbt,n,d : in integer32 ) is

  -- DESCRIPTION :
  --   Generates an n-by-n matrix of series of degree d,
  --   with complex coefficients in standard double precision.
  --   Converts an n-by-n matrix of series of degree d with standard
  --   double precision complex coefficients into a matrix series.

    use Standard_Complex_Series_Matrices;
    use Standard_Series_Matrix_Solvers;

    sA : constant Standard_Complex_Series_Matrices.Matrix(1..n,1..n)
       := Standard_Random_Series_Matrices.Random_Series_Matrix(1,n,1,n,d);
    As : constant Standard_Complex_Matrix_Series.Matrix 
       := Standard_Complex_Matrix_Series.Create(sA); 
    vm : constant Standard_Complex_VecMats.VecMat(0..As.deg)
       := Series_Coefficient_Vectors.Standard_Series_Coefficients(As);
    sx : constant Standard_Complex_Series_Vectors.Vector(1..n)
       := Standard_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    xs : constant Standard_Complex_Vector_Series.Vector(d)
       := Standard_Complex_Vector_Series.Create(sx);
    sb : constant Standard_Complex_Series_Vectors.Vector(1..n) := sA*sx;
    bs : constant Standard_Complex_Vector_Series.Vector(d)
       := Standard_Complex_Vector_Series.Create(sb);
    sbcff : constant Standard_Complex_VecVecs.VecVec(1..n)
          := Series_Coefficient_Vectors.Standard_Series_Coefficients(sb);
    bscff : constant Standard_Complex_VecVecs.VecVec(0..bs.deg)
          := Series_Coefficient_Vectors.Standard_Series_Coefficients(bs);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    wrk : constant Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..n);
    info : integer32;
    ans : character;
    nbrotp,output : boolean;
    err : double_float;

  begin
    put("Output of numbers ? (y/n) "); Ask_Yes_or_No(ans);
    nbrotp := (ans = 'y');
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    if nbrotp then
      put_line("The coefficients of the matrix series :"); put(As);
      put_line("The coefficient matrices : "); put(vm);
      put_line("The exact solution x :"); put_line(sx);
      put_line("The coefficients of the vector series x :"); put(xs);
      put_line("The right hand side vector b :"); put_line(sb);
      put_line("The coefficients of b : "); put_line(sbcff);
      put_line("The coefficients of the vector series b :"); put(bs);
      put_line("The coefficients of the vector series b :"); put_line(bscff);
    end if;
    if nbt > 1 then
      Multitasked_Solve_by_lufac(nbt,vm,bscff,ipvt,info,output);
    elsif nbt = 1 then
      put("Run multitasked code ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Multitasked_Solve_by_lufac(nbt,vm,bscff,ipvt,info,output);
       else Solve_by_lufac(vm,bscff,ipvt,info,wrk);
      end if;
    end if;
    put("info : "); put(info,1); new_line;
    err := Error(xs,bscff,nbrotp);
    put("Sum of errors :"); put(err,3); new_line;
  end Standard_Test;

  procedure DoblDobl_Test ( nbt,n,d : in integer32 ) is

  -- DESCRIPTION :
  --   Generates an n-by-n matrix of series of degree d,
  --   with complex coefficients in standard double precision.
  --   Converts an n-by-n matrix of series of degree d with double
  --   double precision complex coefficients into a matrix series.

    use DoblDobl_Complex_Series_Matrices;
    use DoblDobl_Series_Matrix_Solvers;

    sA : constant DoblDobl_Complex_Series_Matrices.Matrix(1..n,1..n)
       := DoblDobl_Random_Series_Matrices.Random_Series_Matrix(1,n,1,n,d);
    As : constant DoblDobl_Complex_Matrix_Series.Matrix 
       := DoblDobl_Complex_Matrix_Series.Create(sA); 
    vm : constant DoblDobl_Complex_VecMats.VecMat(0..As.deg)
       := Series_Coefficient_Vectors.DoblDobl_Series_Coefficients(As);
    sx : constant DoblDobl_Complex_Series_Vectors.Vector(1..n)
       := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    xs : constant DoblDobl_Complex_Vector_Series.Vector(d)
       := DoblDobl_Complex_Vector_Series.Create(sx);
    sb : constant DoblDobl_Complex_Series_Vectors.Vector(1..n) := sA*sx;
    bs : constant DoblDobl_Complex_Vector_Series.Vector(d)
       := DoblDobl_Complex_Vector_Series.Create(sb);
    sbcff : constant DoblDobl_Complex_VecVecs.VecVec(1..n)
          := Series_Coefficient_Vectors.DoblDobl_Series_Coefficients(sb);
    bscff : constant DoblDobl_Complex_VecVecs.VecVec(0..bs.deg)
          := Series_Coefficient_Vectors.DoblDobl_Series_Coefficients(bs);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    wrk : constant DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector(1..n);
    info : integer32;
    ans : character;
    nbrotp,output : boolean;
    err : double_double;

  begin
    put("Output of numbers ? (y/n) "); Ask_Yes_or_No(ans);
    nbrotp := (ans = 'y');
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    if nbrotp then
      put_line("The coefficients of the matrix series :"); put(As);
      put_line("The coefficient matrices : "); put(vm);
      put_line("The exact solution x :"); put_line(sx);
      put_line("The coefficients of the vector series x :"); put(xs);
      put_line("The right hand side vector b :"); put_line(sb);
      put_line("The coefficients of b : "); put_line(sbcff);
      put_line("The coefficients of the vector series b :"); put(bs);
      put_line("The coefficients of the vector series b :"); put_line(bscff);
    end if;
    if nbt > 1 then
      Multitasked_Solve_by_lufac(nbt,vm,bscff,ipvt,info,output);
    elsif nbt = 1 then
      put("Run multitasked code ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Multitasked_Solve_by_lufac(nbt,vm,bscff,ipvt,info,output);
       else Solve_by_lufac(vm,bscff,ipvt,info,wrk);
      end if;
    end if;
    put("info : "); put(info,1); new_line;
    err := Error(xs,bscff,nbrotp);
    put("Sum of errors : "); put(err,3); new_line;
  end DoblDobl_Test;

  procedure QuadDobl_Test ( nbt,n,d : in integer32 ) is

  -- DESCRIPTION :
  --   Generates an n-by-n matrix of series of degree d,
  --   with complex coefficients in standard double precision.
  --   Converts an n-by-n matrix of series of degree d with quad
  --   double precision complex coefficients into a matrix series.

    use QuadDobl_Complex_Series_Matrices;
    use QuadDobl_Series_Matrix_Solvers;

    sA : constant QuadDobl_Complex_Series_Matrices.Matrix(1..n,1..n)
       := QuadDobl_Random_Series_Matrices.Random_Series_Matrix(1,n,1,n,d);
    As : constant QuadDobl_Complex_Matrix_Series.Matrix 
       := QuadDobl_Complex_Matrix_Series.Create(sA); 
    vm : constant QuadDobl_Complex_VecMats.VecMat(0..As.deg)
       := Series_Coefficient_Vectors.QuadDobl_Series_Coefficients(As);
    sx : constant QuadDobl_Complex_Series_Vectors.Vector(1..n)
       := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    xs : constant QuadDobl_Complex_Vector_Series.Vector(d)
       := QuadDobl_Complex_Vector_Series.Create(sx);
    sb : constant QuadDobl_Complex_Series_Vectors.Vector(1..n) := sA*sx;
    bs : constant QuadDobl_Complex_Vector_Series.Vector(d)
       := QuadDobl_Complex_Vector_Series.Create(sb);
    sbcff : constant QuadDobl_Complex_VecVecs.VecVec(1..n)
          := Series_Coefficient_Vectors.QuadDobl_Series_Coefficients(sb);
    bscff : constant QuadDobl_Complex_VecVecs.VecVec(0..bs.deg)
          := Series_Coefficient_Vectors.QuadDobl_Series_Coefficients(bs);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    wrk : constant QuadDobl_Complex_Vectors.Link_to_Vector
        := new QuadDobl_Complex_Vectors.Vector(1..n);
    info : integer32;
    ans : character;
    nbrotp,output : boolean;
    err : quad_double;

  begin
    put("Output of numbers ? (y/n) "); Ask_Yes_or_No(ans);
    nbrotp := (ans = 'y');
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    if nbrotp then
      put_line("The coefficients of the matrix series :"); put(As);
      put_line("The coefficient matrices : "); put(vm);
      put_line("The exact solution x :"); put_line(sx);
      put_line("The coefficients of the vector series x :"); put(xs);
      put_line("The right hand side vector b :"); put_line(sb);
      put_line("The coefficients of b : "); put_line(sbcff);
      put_line("The coefficients of the vector series b :"); put(bs);
      put_line("The coefficients of the vector series b :"); put_line(bscff);
    end if;
    if nbt > 1 then
      Multitasked_Solve_by_lufac(nbt,vm,bscff,ipvt,info,output);
    elsif nbt = 1 then
      put("Run multitasked code ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Multitasked_Solve_by_lufac(nbt,vm,bscff,ipvt,info,output);
       else Solve_by_lufac(vm,bscff,ipvt,info,wrk);
      end if;
    end if;
    put("info : "); put(info,1); new_line;
    err := Error(xs,bscff,nbrotp);
    put("Sum of errors : "); put(err,3); new_line;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension of the linear system,
  --   the degrees of the series in the system, and the number of tasks.

    nbt,dim,deg : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Testing the linearization of systems of power series ...");
    put("  Give the number of equations and variables : "); get(dim);
    put("  Give the degree of the series : "); get(deg);
    put("  Give the number of tasks : "); get(nbt);
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Standard_Test(nbt,dim,deg);
      when '1' => DoblDobl_Test(nbt,dim,deg);
      when '2' => QuadDobl_Test(nbt,dim,deg);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtserlin;
