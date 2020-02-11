with text_io;                            use text_io;
with duration_io;
with Ada.Calendar;
with Communications_with_User;           use Communications_with_User;
with Time_Stamps;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecMats;
with QuadDobl_Complex_Vectors_cv;
with Standard_Complex_Series_Vectors;
with Standard_Random_Series_Vectors;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Random_Series_Vectors;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Random_Series_Vectors;
with Series_Coefficient_Vectors;         use Series_Coefficient_Vectors;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Evaluation_Differentiation_Errors;  use Evaluation_Differentiation_Errors;
with Multitasked_AlgoDiff_Convolutions;  use Multitasked_AlgoDiff_Convolutions;

procedure ts_mtadcnv is

-- DESCRIPTION :
--   Development of algorithmic differentiation to evaluate and differentiate
--   polynomial systems at truncated power series with multitasking.

  procedure Standard_Random_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to evaluate and differentiate with multitasking,
  --   in double precision.  Compares with the one task run.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

    use Standard_Speelpenning_Convolutions;

    c : constant Circuits
      := Standard_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    x : constant Standard_Complex_Series_Vectors.Vector(1..dim)
      := Standard_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant Standard_Complex_VecVecs.VecVec(1..dim)
         := Standard_Series_Coefficients(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Maxima(c,dim);
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant Standard_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant Standard_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant Standard_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant Standard_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant Standard_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    err : double_float;
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup : Duration;
    ans : character;

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("Running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    EvalDiff(c,xcff,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      multstart := Ada.Calendar.Clock;
      Standard_Multitasked_EvalDiff(nbt,c,xcff,mxe,pwt,vy2,vm2,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      err := Difference(vy1,vy2);
      put("  the error of evaluation : "); put(err,3); new_line;
      err := Difference(vm1,vm2);
      put("  the error of differentiation : "); put(err,3); new_line;
      put_line("-> Elapsed time on multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        put("The speedup : ");
        duration_io.put(speedup,1,3); new_line;
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end Standard_Random_Test;

  procedure DoblDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to evaluate and differentiate with multitasking,
  --   in double double precision.  Compares with the one task run.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

    use DoblDobl_Speelpenning_Convolutions;

    c : constant Circuits
      := DoblDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    x : constant DoblDobl_Complex_Series_Vectors.Vector(1..dim)
      := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
         := DoblDobl_Series_Coefficients(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Maxima(c,dim);
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant DoblDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant DoblDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant DoblDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant DoblDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant DoblDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    err : double_double;
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup : Duration;
    ans : character;

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("Running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    EvalDiff(c,xcff,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    new_line;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      multstart := Ada.Calendar.Clock;
      DoblDobl_Multitasked_EvalDiff(nbt,c,xcff,mxe,pwt,vy2,vm2,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      err := Difference(vy1,vy2);
      put("  the error of evaluation : "); put(err,3); new_line;
      err := Difference(vm1,vm2);
      put("  the error of differentiation : "); put(err,3); new_line;
      put_line("-> Elapsed time on multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        put("The speedup : ");
        duration_io.put(speedup,1,3); new_line;
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end DoblDobl_Random_Test;

  procedure QuadDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32; nbt : in out integer32;
                output : in boolean := true ) is

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to evaluate and differentiate with multitasking,
  --   in double double precision.  Compares with the one task run.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

    use QuadDobl_Speelpenning_Convolutions;

    c : constant Circuits
      := QuadDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    x : constant QuadDobl_Complex_Series_Vectors.Vector(1..dim)
      := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
         := QuadDobl_Series_Coefficients(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Maxima(c,dim);
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant QuadDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant QuadDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant QuadDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant QuadDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant QuadDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    err : quad_double;
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup : Duration;
    ans : character;

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("Running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    EvalDiff(c,xcff,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    loop
      new_line;
      put("Running with "); put(nbt,1); put_line(" tasks ...");
      multstart := Ada.Calendar.Clock;
      QuadDobl_Multitasked_EvalDiff(nbt,c,xcff,mxe,pwt,vy2,vm2,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      err := Difference(vy1,vy2);
      put("  the error of evaluation : "); put(err,3); new_line;
      err := Difference(vm1,vm2);
      put("  the error of differentiation : "); put(err,3); new_line;
      new_line;
      put_line("-> Elapsed time on multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        put("The speedup : ");
        duration_io.put(speedup,1,3); new_line;
      end if;
      put("Continue with different number of tasks ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give the number of tasks : "); get(nbt);
    end loop;
  end QuadDobl_Random_Test;

  procedure Standard_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in Standard_Speelpenning_Convolutions.Circuits;
                x : in Standard_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Runs a benchmark in double precision.

  -- ON ENTRY :
  --   file     opened for output, to write timings and speedups;
  --   deg      degree of the power series;
  --   nbruns   number of multitasked runs;
  --   inc      increment on the number of tasks;
  --   c        some circuits;
  --   x        where to evaluate c;
  --   mxe      maximum exponent indices;
  --   verbose  if true, then timings are written to screen,
  --            otherwise, no output is written to screen.

    use Standard_Speelpenning_Convolutions;

    dim : constant integer32 := x'last;
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant Standard_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant Standard_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant Standard_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant Standard_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant Standard_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    put_line(file,"double double precision");
    if verbose
     then put_line("Running on one task ...");
    end if;
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,x);
    EvalDiff(c,x,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      multstart := Ada.Calendar.Clock;
      Standard_Multitasked_EvalDiff(nbt,c,x,mxe,pwt,vy2,vm2,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        new_line;
        put_line("-> Elapsed time on multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        if verbose then
          put("The speedup : ");
          duration_io.put(speedup,1,3); new_line;
        end if;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mult_elapsed,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3); new_line(file);
      nbt := nbt + inc;
    end loop;
  end Standard_Benchmark;

  procedure DoblDobl_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Runs a benchmark in double double precision.

  -- ON ENTRY :
  --   file     opened for output, to write timings and speedups;
  --   deg      degree of the power series;
  --   nbruns   number of multitasked runs;
  --   inc      increment on the number of tasks;
  --   c        some circuits;
  --   x        where to evaluate c;
  --   mxe      maximum exponent indices;
  --   verbose  if true, then timings are written to screen,
  --            otherwise, no output is written to screen.

    use DoblDobl_Speelpenning_Convolutions;

    dim : constant integer32 := x'last;
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant DoblDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant DoblDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant DoblDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant DoblDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant DoblDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    put_line(file,"double double precision");
    if verbose
     then put_line("Running on one task ...");
    end if;
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,x);
    EvalDiff(c,x,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      multstart := Ada.Calendar.Clock;
      DoblDobl_Multitasked_EvalDiff(nbt,c,x,mxe,pwt,vy2,vm2,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        new_line;
        put_line("-> Elapsed time on multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        if verbose then
          put("The speedup : ");
          duration_io.put(speedup,1,3); new_line;
        end if;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mult_elapsed,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3); new_line(file);
      nbt := nbt + inc;
    end loop;
  end DoblDobl_Benchmark;

  procedure QuadDobl_Benchmark
              ( file : in file_type; deg,nbruns,inc : in integer32;
                c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Runs a benchmark in quad double precision.

  -- ON ENTRY :
  --   file     opened for output, to write timings and speedups;
  --   deg      degree of the power series;
  --   nbruns   number of multitasked runs;
  --   inc      increment on the number of tasks;
  --   c        some circuits;
  --   x        where to evaluate c;
  --   mxe      maximum exponent indices;
  --   verbose  if true, then timings are written to screen,
  --            otherwise, no output is written to screen.

    use QuadDobl_Speelpenning_Convolutions;

    dim : constant integer32 := x'last;
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant QuadDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant QuadDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant QuadDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant QuadDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant QuadDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    put_line(file,"quad double precision");
    if verbose
     then put_line("Running on one task ...");
    end if;
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,x);
    EvalDiff(c,x,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,seri_elapsed,1,3); new_line(file);
    if verbose then
      put_line("-> Elapsed time without multitasking : ");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    end if;   
    for k in 1..nbruns loop
      if verbose then
        new_line;
        put("Running with "); put(nbt,1); put_line(" tasks ...");
      end if;
      multstart := Ada.Calendar.Clock;
      QuadDobl_Multitasked_EvalDiff(nbt,c,x,mxe,pwt,vy2,vm2,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        new_line;
        put_line("-> Elapsed time on multitasking : ");
        Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        if verbose then
          put("The speedup : ");
          duration_io.put(speedup,1,3); new_line;
        end if;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mult_elapsed,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3); new_line(file);
      nbt := nbt + inc;
    end loop;
  end QuadDobl_Benchmark;

  procedure Benchmark ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates random circuits in quad double precision,
  --   and runs benchmark tests in all three precisions.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    qd_c : constant QuadDobl_Speelpenning_Convolutions.Circuits
         := QuadDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    dd_c : constant DoblDobl_Speelpenning_Convolutions.Circuits
         := to_double_double(qd_c);
    d_c : constant Standard_Speelpenning_Convolutions.Circuits
        := to_double(qd_c);
    qd_x : constant QuadDobl_Complex_Series_Vectors.Vector(1..dim)
         := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    qd_xcff : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
            := QuadDobl_Series_Coefficients(qd_x);
    dd_xcff : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
            := QuadDobl_Complex_Vectors_cv.to_double_double(qd_xcff);
    d_xcff : constant Standard_Complex_VecVecs.VecVec(1..dim)
           := QuadDobl_Complex_Vectors_cv.to_double(qd_xcff);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := QuadDobl_Speelpenning_Convolutions.Exponent_Maxima(qd_c,dim);
    nbruns,inc : integer32 := 0;
    file : file_type;

  begin
    put("Give the number of multitasked runs : "); get(nbruns);
    put("Give the increment on the tasks : "); get(inc);
    skip_line;
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    Standard_Benchmark(file,deg,nbruns,inc,d_c,d_xcff,mxe,false);
    DoblDobl_Benchmark(file,deg,nbruns,inc,dd_c,dd_xcff,mxe,false);
    QuadDobl_Benchmark(file,deg,nbruns,inc,qd_c,qd_xcff,mxe,false);
  end Benchmark;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the number of tasks,
  --   for the dimensions of the problem,
  --   generates then the data for the problem, and
  --   launches the tasks.

    dim,deg,nbr,pwr,nbt : integer32 := 0;
    answer,precision : character;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree of the series : "); get(deg);
    put("Give the number of monomials : "); get(nbr);
    put("Give the highest power of each variable : "); get(pwr);
    put("Benchmarking for all precisions ? (y/n) "); Ask_Yes_or_No(answer);
    if answer  = 'y' then
      Benchmark(dim,deg,nbr,pwr);
    else
      put("Give the number of tasks : "); get(nbt);
      put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(answer);
      new_line;
      put_line("MENU for the working precision :");
      put_line("  0. standard double precision");
      put_line("  1. double double precision");
      put_line("  2. quad double precision");
      put("Type 0, 1, or 2 to select the precision : ");
      Ask_Alternative(precision,"012");
      new_line;
      case precision is
        when '0' => Standard_Random_Test(dim,deg,nbr,pwr,nbt,answer = 'y');
        when '1' => DoblDobl_Random_Test(dim,deg,nbr,pwr,nbt,answer = 'y');
        when '2' => QuadDobl_Random_Test(dim,deg,nbr,pwr,nbt,answer = 'y');
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_mtadcnv;
