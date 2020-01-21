with text_io;                            use text_io;
with duration_io;
with Ada.Calendar;
with Time_Stamps;
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

  procedure Show_Speedup ( serial,multi : in duration ) is

  -- DESCRIPTION :
  --   On input are the elapsed serial and multitasked times.
  --   If serial /= 0.0, then the speedup is shown.

    speedup : Duration;
 
  begin
    if serial + 1.0 /= 1.0 then
      speedup := serial/multi;
      put("The speedup : ");
      duration_io.put(speedup,1,3); new_line;
    end if;
  end Show_Speedup;

  procedure Standard_Run
              ( nbt,dim : in integer32;
                vm : in Standard_Complex_VecMats.VecMat;
                vb : in Standard_Complex_VecVecs.VecVec;
                xs : in Standard_Complex_Vector_Series.Vector;
                mltelp,serelp : in out Duration;
                output,nbrotp : in boolean ) is

  -- DESCRIPTION :
  --   Does a run with nbt tasks in double precision.
  --   Solves the linearized matrix series system defined by vm and vb.
  --   Prints the elapsed time and if defined, the speedup.
  --   Prints the difference between the computed and the generated solution.

  -- ON ENTRY :
  --   nbt      the number of tasks is 1 or larger;
  --   dim      dimension of the matrices in vm;
  --   vm       matrices in the series equation;
  --   vb       right hand side vector of the equation;
  --   xs       the generated solution to the matrix series equation;
  --   serelp   the previous elapsed wall clock time of a serial run;
  --   mltelp   the previous elapsed wall clock time of a multitasked run;
  --   output   if true, the multitasked run is verbose, else silent;
  --   nbrotp   if true, all numbers are shown, else not.

  -- ON RETURN :
  --   serelp   updated elapsed wall clock time of a serial run,
  --            if nbt = 1 and the user did not want multitasking;
  --   mltelp   updated elapsed wall clock time of a multitasked run,
  --            if nbt > 1.

    ans : character;
    err : double_float;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;
    wrk : constant Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..dim);
    wks : constant Standard_Complex_VecVecs.VecVec(1..nbt)
        := Allocate_Work_Space(nbt,dim);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;

    use Ada.Calendar;

  begin
    if nbt > 1 then
      multstart := Ada.Calendar.Clock;
      Multitasked_Solve_by_lufac(nbt,vm,vb,ipvt,info,wks,output);
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
      put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      Show_Speedup(serelp,mltelp);
    elsif nbt = 1 then
      put("Run multitasked code ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        multstart := Ada.Calendar.Clock;
        Multitasked_Solve_by_lufac(nbt,vm,vb,ipvt,info,wks,output);
        multstop := Ada.Calendar.Clock;
        mltelp := multstop - multstart;
        put_line("-> Elapsed time with one task :");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
        Show_Speedup(serelp,mltelp);
      else
        seristart := Ada.Calendar.Clock;
        Standard_Series_Matrix_Solvers.Solve_by_lufac(vm,vb,ipvt,info,wrk);
        seristop := Ada.Calendar.Clock;
        serelp := seristop - seristart;
        put_line("-> Elapsed time without multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      end if;
    end if;
    put("info : "); put(info,1); new_line;
    err := Error(xs,vb,nbrotp);
    put("Sum of errors :"); put(err,3); new_line;
  end Standard_Run;

  procedure DoblDobl_Run
              ( nbt,dim : in integer32;
                vm : in DoblDobl_Complex_VecMats.VecMat;
                vb : in DoblDobl_Complex_VecVecs.VecVec;
                xs : in DoblDobl_Complex_Vector_Series.Vector;
                mltelp,serelp : in out Duration;
                output,nbrotp : in boolean ) is

  -- DESCRIPTION :
  --   Does a run with nbt tasks in double double precision.
  --   Solves the linearized matrix series system defined by vm and vb.
  --   Prints the elapsed time and if defined, the speedup.
  --   Prints the difference between the computed and the generated solution.

  -- ON ENTRY :
  --   nbt      the number of tasks is 1 or larger;
  --   dim      dimension of the matrices in vm;
  --   vm       matrices in the series equation;
  --   vb       right hand side vector of the equation;
  --   xs       the generated solution to the matrix series equation;
  --   serelp   the previous elapsed wall clock time of a serial run;
  --   mltelp   the previous elapsed wall clock time of a multitasked run;
  --   output   if true, the multitasked run is verbose, else silent;
  --   nbrotp   if true, all numbers are shown, else not.

  -- ON RETURN :
  --   serelp   updated elapsed wall clock time of a serial run,
  --            if nbt = 1 and the user did not want multitasking;
  --   mltelp   updated elapsed wall clock time of a multitasked run,
  --            if nbt > 1.

    ans : character;
    err : double_double;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;
    wrk : constant DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector(1..dim);
    wks : constant DoblDobl_Complex_VecVecs.VecVec(1..nbt)
        := Allocate_Work_Space(nbt,dim);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;

    use Ada.Calendar;

  begin
    if nbt > 1 then
      multstart := Ada.Calendar.Clock;
      Multitasked_Solve_by_lufac(nbt,vm,vb,ipvt,info,wks,output);
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
      put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      Show_Speedup(serelp,mltelp);
    elsif nbt = 1 then
      put("Run multitasked code ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        multstart := Ada.Calendar.Clock;
        Multitasked_Solve_by_lufac(nbt,vm,vb,ipvt,info,wks,output);
        multstop := Ada.Calendar.Clock;
        mltelp := multstop - multstart;
        put_line("-> Elapsed time with one task :");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
        Show_Speedup(serelp,mltelp);
      else
        seristart := Ada.Calendar.Clock;
        DoblDobl_Series_Matrix_Solvers.Solve_by_lufac(vm,vb,ipvt,info,wrk);
        seristop := Ada.Calendar.Clock;
        serelp := seristop - seristart;
        put_line("-> Elapsed time without multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      end if;
    end if;
    put("info : "); put(info,1); new_line;
    err := Error(xs,vb,nbrotp);
    put("Sum of errors : "); put(err,3); new_line;
  end DoblDobl_Run;

  procedure QuadDobl_Run
              ( nbt,dim : in integer32;
                vm : in QuadDobl_Complex_VecMats.VecMat;
                vb : in QuadDobl_Complex_VecVecs.VecVec;
                xs : in QuadDobl_Complex_Vector_Series.Vector;
                mltelp,serelp : in out Duration;
                output,nbrotp : in boolean ) is

  -- DESCRIPTION :
  --   Does a run with nbt tasks in double double precision.
  --   Solves the linearized matrix series system defined by vm and vb.
  --   Prints the elapsed time and if defined, the speedup.
  --   Prints the difference between the computed and the generated solution.

  -- ON ENTRY :
  --   nbt      the number of tasks is 1 or larger;
  --   dim      dimension of the matrices in vm;
  --   vm       matrices in the series equation;
  --   vb       right hand side vector of the equation;
  --   xs       the generated solution to the matrix series equation;
  --   serelp   the previous elapsed wall clock time of a serial run;
  --   mltelp   the previous elapsed wall clock time of a multitasked run;
  --   output   if true, the multitasked run is verbose, else silent;
  --   nbrotp   if true, all numbers are shown, else not.

  -- ON RETURN :
  --   serelp   updated elapsed wall clock time of a serial run,
  --            if nbt = 1 and the user did not want multitasking;
  --   mltelp   updated elapsed wall clock time of a multitasked run,
  --            if nbt > 1.

    ans : character;
    err : quad_double;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;
    wrk : constant QuadDobl_Complex_Vectors.Link_to_Vector
        := new QuadDobl_Complex_Vectors.Vector(1..dim);
    wks : constant QuadDobl_Complex_VecVecs.VecVec(1..nbt)
        := Allocate_Work_Space(nbt,dim);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;

    use Ada.Calendar;

  begin
    if nbt > 1 then
      multstart := Ada.Calendar.Clock;
      Multitasked_Solve_by_lufac(nbt,vm,vb,ipvt,info,wks,output);
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
      put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      Show_Speedup(serelp,mltelp);
    elsif nbt = 1 then
      put("Run multitasked code ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        multstart := Ada.Calendar.Clock;
        Multitasked_Solve_by_lufac(nbt,vm,vb,ipvt,info,wks,output);
        multstop := Ada.Calendar.Clock;
        mltelp := multstop - multstart;
        put_line("-> Elapsed time with one task :");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
        Show_Speedup(serelp,mltelp);
      else
        seristart := Ada.Calendar.Clock;
        QuadDobl_Series_Matrix_Solvers.Solve_by_lufac(vm,vb,ipvt,info,wrk);
        seristop := Ada.Calendar.Clock;
        serelp := seristop - seristart;
        put_line("-> Elapsed time without multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      end if;
    end if;
    put("info : "); put(info,1); new_line;
    err := Error(xs,vb,nbrotp);
    put("Sum of errors : "); put(err,3); new_line;
  end QuadDobl_Run;

  procedure Standard_Test ( n,d : in integer32 ) is

  -- DESCRIPTION :
  --   Generates an n-by-n matrix of series of degree d,
  --   with complex coefficients and a solution of the same dimension
  --   and degree, in double precision,
  --   Prompts then the user for the number of tasks and runs the test.

    use Standard_Complex_Series_Matrices; -- for the sA*sx operation

    nbt : integer32 := 0;
    sA : constant Standard_Complex_Series_Matrices.Matrix(1..n,1..n)
       := Standard_Random_Series_Matrices.Random_Series_Matrix(1,n,1,n,d);
    As : constant Standard_Complex_Matrix_Series.Matrix 
       := Standard_Complex_Matrix_Series.Create(sA); 
    vmbackup : constant Standard_Complex_VecMats.VecMat(0..As.deg)
             := Series_Coefficient_Vectors.Standard_Series_Coefficients(As);
    vm : Standard_Complex_VecMats.VecMat(vmbackup'range);
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
    bswrk : Standard_Complex_VecVecs.VecVec(bscff'range);
    ans : character;
    nbrotp,output : boolean;
    mult_elapsed,seri_elapsed : Duration := 0.0;

  begin
    put("Output of numbers ? (y/n) "); Ask_Yes_or_No(ans);
    nbrotp := (ans = 'y');
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    if nbrotp then
      put_line("The coefficients of the matrix series :"); put(As);
      put_line("The coefficient matrices : "); put(vmbackup);
      put_line("The exact solution x :"); put_line(sx);
      put_line("The coefficients of the vector series x :"); put(xs);
      put_line("The right hand side vector b :"); put_line(sb);
      put_line("The coefficients of b : "); put_line(sbcff);
      put_line("The coefficients of the vector series b :"); put(bs);
      put_line("The coefficients of the vector series b :"); put_line(bscff);
    end if;
    loop
      new_line;
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      Standard_Complex_VecMats.Copy(vmbackup,vm);
      Standard_Complex_VecVecs.Copy(bscff,bswrk);
      Standard_Run(nbt,n,vm,bswrk,xs,mult_elapsed,seri_elapsed,output,nbrotp);
    end loop;
  end Standard_Test;

  procedure DoblDobl_Test ( n,d : in integer32 ) is

  -- DESCRIPTION :
  --   Generates an n-by-n matrix of series of degree d,
  --   with complex coefficients and a solution of the same dimension
  --   and degree, in double double precision,
  --   Prompts then the user for the number of tasks and runs the test.

    use DoblDobl_Complex_Series_Matrices; -- for the sA*sx operation

    nbt : integer32 := 0;
    sA : constant DoblDobl_Complex_Series_Matrices.Matrix(1..n,1..n)
       := DoblDobl_Random_Series_Matrices.Random_Series_Matrix(1,n,1,n,d);
    As : constant DoblDobl_Complex_Matrix_Series.Matrix 
       := DoblDobl_Complex_Matrix_Series.Create(sA); 
    vmbackup : constant DoblDobl_Complex_VecMats.VecMat(0..As.deg)
             := Series_Coefficient_Vectors.DoblDobl_Series_Coefficients(As);
    vm : DoblDobl_Complex_VecMats.VecMat(vmbackup'range);
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
    bswrk : DoblDobl_Complex_VecVecs.VecVec(bscff'range);
    ans : character;
    nbrotp,output : boolean;
    mult_elapsed,seri_elapsed : Duration := 0.0;

  begin
    put("Output of numbers ? (y/n) "); Ask_Yes_or_No(ans);
    nbrotp := (ans = 'y');
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    if nbrotp then
      put_line("The coefficients of the matrix series :"); put(As);
      put_line("The coefficient matrices : "); put(vmbackup);
      put_line("The exact solution x :"); put_line(sx);
      put_line("The coefficients of the vector series x :"); put(xs);
      put_line("The right hand side vector b :"); put_line(sb);
      put_line("The coefficients of b : "); put_line(sbcff);
      put_line("The coefficients of the vector series b :"); put(bs);
      put_line("The coefficients of the vector series b :"); put_line(bscff);
    end if;
    loop
      new_line;
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      DoblDobl_Complex_VecMats.Copy(vmbackup,vm);
      DoblDobl_Complex_VecVecs.Copy(bscff,bswrk);
      DoblDobl_Run(nbt,n,vm,bswrk,xs,mult_elapsed,seri_elapsed,output,nbrotp);
    end loop;
  end DoblDobl_Test;

  procedure QuadDobl_Test ( n,d : in integer32 ) is

  -- DESCRIPTION :
  --   Generates an n-by-n matrix of series of degree d,
  --   with complex coefficients and a solution of the same dimension
  --   and degree, in quad double precision,
  --   Prompts then the user for the number of tasks.

    use QuadDobl_Complex_Series_Matrices; -- for the sA*sx operation

    nbt : integer32 := 0;
    sA : constant QuadDobl_Complex_Series_Matrices.Matrix(1..n,1..n)
       := QuadDobl_Random_Series_Matrices.Random_Series_Matrix(1,n,1,n,d);
    As : constant QuadDobl_Complex_Matrix_Series.Matrix 
       := QuadDobl_Complex_Matrix_Series.Create(sA); 
    vmbackup : constant QuadDobl_Complex_VecMats.VecMat(0..As.deg)
             := Series_Coefficient_Vectors.QuadDobl_Series_Coefficients(As);
    vm : QuadDobl_Complex_VecMats.VecMat(vmbackup'range);
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
    bswrk : QuadDobl_Complex_VecVecs.VecVec(bscff'range);
    ans : character;
    nbrotp,output : boolean;
    seri_elapsed,mult_elapsed : Duration := 0.0;

  begin
    put("Output of numbers ? (y/n) "); Ask_Yes_or_No(ans);
    nbrotp := (ans = 'y');
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    if nbrotp then
      put_line("The coefficients of the matrix series :"); put(As);
      put_line("The coefficient matrices : "); put(vmbackup);
      put_line("The exact solution x :"); put_line(sx);
      put_line("The coefficients of the vector series x :"); put(xs);
      put_line("The right hand side vector b :"); put_line(sb);
      put_line("The coefficients of b : "); put_line(sbcff);
      put_line("The coefficients of the vector series b :"); put(bs);
      put_line("The coefficients of the vector series b :"); put_line(bscff);
    end if;
    loop
      new_line;
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt = 0);
      QuadDobl_Complex_VecMats.Copy(vmbackup,vm);
      QuadDobl_Complex_VecVecs.Copy(bscff,bswrk);
      QuadDobl_Run(nbt,n,vm,bswrk,xs,mult_elapsed,seri_elapsed,output,nbrotp);
    end loop;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension of the linear system,
  --   the degrees of the series in the system, and the number of tasks.

    dim,deg : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Testing the linearization of systems of power series ...");
    put("  Give the number of equations and variables : "); get(dim);
    put("  Give the degree of the series : "); get(deg);
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Standard_Test(dim,deg);
      when '1' => DoblDobl_Test(dim,deg);
      when '2' => QuadDobl_Test(dim,deg);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtserlin;
