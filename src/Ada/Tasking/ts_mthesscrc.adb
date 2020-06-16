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
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_VecMats;
with QuadDobl_Complex_Vectors_cv;
with Standard_Random_Vectors;
with DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;
with Standard_Vector_Splitters;
with DoblDobl_Vector_Splitters;
with QuadDobl_Vector_Splitters;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_System_and_Solutions_io;
with Evaluation_Differentiation_Errors;  use Evaluation_Differentiation_Errors;
with Standard_Complex_Circuits;
with Standard_Coefficient_Circuits;
with Standard_Circuit_Splitters;
with Standard_Circuit_Makers;
with DoblDobl_Complex_Circuits;
with DoblDobl_Circuit_Makers;
with QuadDobl_Complex_Circuits;
with QuadDobl_Circuit_Makers;
with Multitasked_Hessian_Circuits;

procedure ts_mthesscrc is

-- DESCRIPTION :
--   Tests the Hessian criterion, for systems of complex circuits,
--   in double, double double, and quad double arithmetic,
--   with multitasking for shared memory parallel computers.

  procedure Write_Singular_Values 
              ( values : in Standard_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Writes the singular values in values, line by line.

    lnk : Standard_Complex_Vectors.Link_to_Vector;
    val : double_float;

  begin
    for k in values'range loop
      lnk := values(k);
      for i in lnk'first..lnk'last-1 loop
        val := Standard_Complex_Numbers.REAL_PART(lnk(i));
        put(val,3);
      end loop;
      new_line;
    end loop;
  end Write_Singular_Values;

  procedure Write_Singular_Values 
              ( values : in DoblDobl_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Writes the singular values in values, line by line.

    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    val : double_double;

  begin
    for k in values'range loop
      lnk := values(k);
      for i in lnk'first..lnk'last-1 loop
        val := DoblDobl_Complex_Numbers.REAL_PART(lnk(i));
        put(" "); put(val,3);
      end loop;
      new_line;
    end loop;
  end Write_Singular_Values;

  procedure Write_Singular_Values 
              ( values : in QuadDobl_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Writes the singular values in values, line by line.

    lnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    val : quad_double;

  begin
    for k in values'range loop
      lnk := values(k);
      for i in lnk'first..lnk'last-1 loop
        val := QuadDobl_Complex_Numbers.REAL_PART(lnk(i));
        put(" "); put(val,3);
      end loop;
      new_line;
    end loop;
  end Write_Singular_Values;

  procedure Standard_Test
              ( s : in Standard_Complex_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                static,output : in boolean ) is

  -- DESCRIPTION :
  --   Runs the test in double precision.

  -- ON ENTRY :
  --   s        a system of convolution circuits;
  --   x        coefficients of a start solution;
  --   static   to apply static load balancing;
  --   output   to see all singular values and
  --            for output during the multitasked runs.

    dim : constant integer32 := s.dim;
    U,V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    e : Standard_Complex_Vectors.Vector(1..dim);
    svl : constant Standard_Complex_VecVecs.VecVec(0..s.neq)
        := Standard_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant Standard_Complex_VecMats.VecMat(1..s.neq)
       := Standard_Complex_Circuits.Allocate(s.neq,s.dim);
    values : Standard_Complex_VecVecs.VecVec(0..s.neq)
           := Standard_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    nbt : integer32 := 0;
    err,eta : double_float;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;

    use Ada.Calendar;

  begin
    put_line("Computing first without multitasking ...");
    seristart := Ada.Calendar.Clock;
    Standard_Complex_Circuits.Singular_Values(s,x,vh,U,V,e,svl);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking :");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    if output
     then put_line("All singular values :"); Write_Singular_Values(svl);
    end if;
    loop
      new_line;
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt <= 0);
      multstart := Ada.Calendar.Clock;
      Multitasked_Hessian_Circuits.Multitasked_Singular_Values
        (nbt,s,x,values,static,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if output then
        put_line("All singular values computed by multitasking :");
        Write_Singular_Values(values);
      end if;
      err := Difference(svl,values);
      eta := Multitasked_Hessian_Circuits.Standard_Distance(svl);
      put("The difference error : "); put(err,3); new_line;
      put("Estimate for the distance : "); put(eta,3); new_line;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("The speedup : "); duration_io.put(speedup,1,3);
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
    end loop;
  end Standard_Test;

  procedure Standard_Coefficient_Test
              ( s : in Standard_Coefficient_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                output : in boolean ) is

  -- DESCRIPTION :
  --   Runs the test in double precision.

  -- ON ENTRY :
  --   s        a system of convolution circuits;
  --   x        coefficients of a start solution;
  --   output   to see all singular values and
  --            for output during the multitasked runs.

    dim : constant integer32 := s.dim;
    xr,xi : Standard_Floating_Vectors.Link_to_Vector;
    U,V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    e : Standard_Complex_Vectors.Vector(1..dim);
    svl : constant Standard_Complex_VecVecs.VecVec(0..s.neq)
        := Standard_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant Standard_Complex_VecMats.VecMat(1..s.neq)
       := Standard_Complex_Circuits.Allocate(s.neq,s.dim);
    values : Standard_Complex_VecVecs.VecVec(0..s.neq)
           := Standard_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    nbt : integer32 := 0;
    err,eta : double_float;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;

    use Ada.Calendar;

  begin
    Standard_Vector_Splitters.Split_Complex(x,xr,xi);
    put_line("Computing first without multitasking ...");
    seristart := Ada.Calendar.Clock;
    Standard_Coefficient_Circuits.Singular_Values(s,xr,xi,vh,U,V,e,svl);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking :");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    if output
     then put_line("All singular values :"); Write_Singular_Values(svl);
    end if;
    loop
      new_line;
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt <= 0);
      multstart := Ada.Calendar.Clock;
      Multitasked_Hessian_Circuits.Multitasked_Singular_Values
        (nbt,s,xr,xi,values,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if output then
        put_line("All singular values computed by multitasking :");
        Write_Singular_Values(values);
      end if;
      err := Difference(svl,values);
      eta := Multitasked_Hessian_Circuits.Standard_Distance(svl);
      put("The difference error : "); put(err,3); new_line;
      put("Estimate for the distance : "); put(eta,3); new_line;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("The speedup : "); duration_io.put(speedup,1,3);
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
    end loop;
  end Standard_Coefficient_Test;

  procedure DoblDobl_Test
              ( s : in DoblDobl_Complex_Circuits.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                static,output : in boolean ) is

  -- DESCRIPTION :
  --   Runs the test in double double precision.

  -- ON ENTRY :
  --   s        a system of convolution circuits;
  --   x        coefficients of a start solution;
  --   static   to apply static load balancing;
  --   output   to see all singular values and
  --            for output during the multitasked runs.

    dim : constant integer32 := s.dim;
    U,V : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : DoblDobl_Complex_Vectors.Vector(1..dim);
    svl : constant DoblDobl_Complex_VecVecs.VecVec(0..s.neq)
        := DoblDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant DoblDobl_Complex_VecMats.VecMat(1..s.neq)
       := DoblDobl_Complex_Circuits.Allocate(s.neq,s.dim);
    values : DoblDobl_Complex_VecVecs.VecVec(0..s.neq)
           := DoblDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    nbt : integer32 := 0;
    err,eta : double_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;

    use Ada.Calendar;

  begin
    put_line("Computing first without multitasking ...");
    seristart := Ada.Calendar.Clock;
    DoblDobl_Complex_Circuits.Singular_Values(s,x,vh,U,V,e,svl);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking :");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    if output
     then put_line("All singular values :"); Write_Singular_Values(svl);
    end if;
    loop
      new_line;
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt <= 0);
      multstart := Ada.Calendar.Clock;
      Multitasked_Hessian_Circuits.Multitasked_Singular_Values
        (nbt,s,x,values,static,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if output then
        put_line("All singular values computed by multitasking :");
        Write_Singular_Values(values);
      end if;
      err := Difference(svl,values);
      eta := Multitasked_Hessian_Circuits.DoblDobl_Distance(svl);
      put("The difference error : "); put(err,3); new_line;
      put("Estimate for the distance : "); put(eta,3); new_line;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("The speedup : "); duration_io.put(speedup,1,3);
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
    end loop;
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( s : in QuadDobl_Complex_Circuits.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                static,output : in boolean ) is

  -- DESCRIPTION :
  --   Runs the test in quad double precision.

  -- ON ENTRY :
  --   s        a system of convolution circuits;
  --   x        coefficients of a start solution;
  --   static   to apply static load balancing;
  --   output   to see all singular values and
  --            for output during the multitasked runs.

    dim : constant integer32 := s.dim;
    U,V : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : QuadDobl_Complex_Vectors.Vector(1..dim);
    svl : constant QuadDobl_Complex_VecVecs.VecVec(0..s.neq)
        := QuadDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant QuadDobl_Complex_VecMats.VecMat(1..s.neq)
       := QuadDobl_Complex_Circuits.Allocate(s.neq,s.dim);
    values : QuadDobl_Complex_VecVecs.VecVec(0..s.neq)
           := QuadDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    nbt : integer32 := 0;
    err,eta : quad_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;

    use Ada.Calendar;

  begin
    put_line("Computing first without multitasking ...");
    seristart := Ada.Calendar.Clock;
    QuadDobl_Complex_Circuits.Singular_Values(s,x,vh,U,V,e,svl);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking :");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    if output
     then put_line("All singular values :"); Write_Singular_Values(svl);
    end if;
    loop
      new_line;
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt <= 0);
      multstart := Ada.Calendar.Clock;
      Multitasked_Hessian_Circuits.Multitasked_Singular_Values
        (nbt,s,x,values,static,output);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if output then
        put_line("All singular values computed by multitasking :");
        Write_Singular_Values(values);
      end if;
      err := Difference(svl,values);
      eta := Multitasked_Hessian_Circuits.QuadDobl_Distance(svl);
      put("The difference error : "); put(err,3); new_line;
      put("Estimate for the distance : "); put(eta,3); new_line;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
        put("The speedup : "); duration_io.put(speedup,1,3);
        put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
      end if;
    end loop;
  end QuadDobl_Test;

  procedure Standard_Random_Test
              ( dim,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in double precision
  --   and then launches the test.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    s : Standard_Complex_Circuits.Link_to_System;
    cs : Standard_Coefficient_Circuits.Link_to_System;
    v : constant Standard_Complex_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    x : constant Standard_Complex_Vectors.Link_to_Vector
      := new Standard_Complex_Vectors.Vector'(v);
    ans : character;
    static,output : boolean;

  begin
    new_line;
    put_line("Generating a random complex circuit system ...");
    s := Standard_Circuit_Makers.Random_Complex_System(dim,nbr,dim,pwr);
    new_line;
    put("Extra output wanted ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    new_line;
    put("Test coefficient circuits ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      cs := Standard_Circuit_Splitters.Split(s);
      Standard_Coefficient_Test(cs,x,output);
    else
      new_line;
      put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
      static := (ans = 'y');
      Standard_Test(s,x,static,output);
    end if;
  end Standard_Random_Test;

  procedure DoblDobl_Random_Test
              ( dim,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in double double precision,
  --   and then launches the test.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    s : DoblDobl_Complex_Circuits.Link_to_System;
    v : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    x : constant DoblDobl_Complex_Vectors.Link_to_Vector
      := new DoblDobl_Complex_Vectors.Vector'(v);
    ans : character;
    static,output : boolean;

  begin
    new_line;
    put_line("Generating a random complex circuit system ...");
    s := DoblDobl_Circuit_Makers.Random_Complex_System(dim,nbr,dim,pwr);
    new_line;
    put("Extra output wanted ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
    static := (ans = 'y');
    DoblDobl_Test(s,x,static,output);
  end DoblDobl_Random_Test;

  procedure QuadDobl_Random_Test
              ( dim,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in quad double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    s : QuadDobl_Complex_Circuits.Link_to_System;
    v : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    x : constant QuadDobl_Complex_Vectors.Link_to_Vector
      := new QuadDobl_Complex_Vectors.Vector'(v);
    ans : character;
    static,output : boolean;

  begin
    new_line;
    put_line("Generating a random complex circuit system ...");
    s := QuadDobl_Circuit_Makers.Random_Complex_System(dim,nbr,dim,pwr);
    new_line;
    put("Extra output wanted ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
    static := (ans = 'y');
    QuadDobl_Test(s,x,static,output);
  end QuadDobl_Random_Test;

  procedure Standard_User_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with solutions
  --   and launches the test in double precision.

    use Standard_Complex_Circuits;

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    ans : character;
    static,output : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    Standard_System_and_Solutions_io.get(lp,sols);
    ls := Standard_Complex_Solutions.Head_Of(sols);
    declare
      s : constant Link_to_System
        := Standard_Circuit_Makers.Make_Complex_System(lp);
      x : constant Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector'(ls.v);
    begin
      new_line;
      put("Extra output wanted ? (y/n) "); Ask_Yes_or_No(ans);
      output := (ans = 'y');
      put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
      static := (ans = 'y');
      Standard_Test(s,x,static,output);
    end;
  end Standard_User_Test;

  procedure DoblDobl_User_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with solutions
  --   and launches the test in double double precision.

    use DoblDobl_Complex_Circuits;

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    ans : character;
    static,output : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    DoblDobl_System_and_Solutions_io.get(lp,sols);
    ls := DoblDobl_Complex_Solutions.Head_Of(sols);
    declare
      s : constant Link_to_System
        := DoblDobl_Circuit_Makers.Make_Complex_System(lp);
      x : constant DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector'(ls.v);
    begin
      new_line;
      put("Extra output wanted ? (y/n) "); Ask_Yes_or_No(ans);
      output := (ans = 'y');
      put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
      static := (ans = 'y');
      DoblDobl_Test(s,x,static,output);
    end;
  end DoblDobl_User_Test;

  procedure QuadDobl_User_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with solutions
  --   and launches the test in quad double precision.

    use QuadDobl_Complex_Circuits;

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    ans : character;
    static,output : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    QuadDobl_System_and_Solutions_io.get(lp,sols);
    ls := QuadDobl_Complex_Solutions.Head_Of(sols);
    declare
      s : constant Link_to_System
        := QuadDobl_Circuit_Makers.Make_Complex_System(lp);
      x : constant QuadDobl_Complex_Vectors.Link_to_Vector
        := new QuadDobl_Complex_Vectors.Vector'(ls.v);
    begin
      new_line;
      put("Extra output wanted ? (y/n) "); Ask_Yes_or_No(ans);
      output := (ans = 'y');
      put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
      static := (ans = 'y');
      QuadDobl_Test(s,x,static,output);
    end;
  end QuadDobl_User_Test;

  procedure Standard_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in Standard_Complex_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Runs a benchmark test in double precision.

  -- REQUIRED : if nbruns = 0, then nbtseq /= null.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbruns   the number of multitasked runs,
  --            if zero, then nbtseq will be used;
  --   inc      increment on the number of tasks, if nbruns /= 0;
  --   nbtseq   sequence of number of tasks for multitasked runs;
  --   s        system in one parameter;
  --   x        some point to evaluate at;
  --   verbose  if extra output is needed.

    dim : constant integer32 := s.dim;
    U,V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    e : Standard_Complex_Vectors.Vector(1..dim);
    svl : constant Standard_Complex_VecVecs.VecVec(0..s.neq)
        := Standard_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant Standard_Complex_VecMats.VecMat(1..s.neq)
       := Standard_Complex_Circuits.Allocate(s.neq,s.dim);
    values : Standard_Complex_VecVecs.VecVec(0..s.neq);
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    if verbose
     then put_line("Computing first without multitasking ...");
    end if;
    put_line(file,"double precision"); flush(file);
    seristart := Ada.Calendar.Clock;
    Standard_Complex_Circuits.Singular_Values(s,x,vh,U,V,e,svl);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    if verbose then
      put_line("-> Elapsed time without multitasking :");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      put_line("running in double precision ...");
    end if;
    put(file,"  1 : "); duration_io.put(file,seri_elapsed,1,3);
    put_line(file," : 1.000"); flush(file);
    if nbruns /= 0 then
      for k in 1..nbruns loop
        values := Standard_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
        multstart := Ada.Calendar.Clock;
        Multitasked_Hessian_Circuits.Multitasked_Singular_Values
          (nbt,s,x,values,false,false);
        multstop := Ada.Calendar.Clock;
        mult_elapsed := multstop - multstart;
        if verbose then
          put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
          Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
        end if;
        put(file,nbt,3); put(file," : ");
        duration_io.put(file,mult_elapsed,1,3); put(file," : ");
        if seri_elapsed + 1.0 /= 1.0 then
          speedup := seri_elapsed/mult_elapsed;
          efficiency := speedup/duration(nbt);
          efficiency := duration(100)*efficiency;
          if verbose then
            put("The speedup : "); duration_io.put(speedup,1,3);
            put("  efficiency : "); duration_io.put(efficiency,1,3); new_line;
          end if;
          duration_io.put(file,speedup,1,3); put(file," : ");
          duration_io.put(file,efficiency,2,2); new_line(file); flush(file);
        end if;
        nbt := nbt + inc;
        Standard_Complex_VecVecs.Clear(values);
      end loop;
    else
      for k in nbtseq'range loop
        nbt := nbtseq(k);
        values := Standard_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
        multstart := Ada.Calendar.Clock;
        Multitasked_Hessian_Circuits.Multitasked_Singular_Values
          (nbt,s,x,values,false,false);
        multstop := Ada.Calendar.Clock;
        mult_elapsed := multstop - multstart;
        if verbose then
          put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
          Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
        end if;
        put(file,nbt,3); put(file," : ");
        duration_io.put(file,mult_elapsed,1,3); put(file," : ");
        if seri_elapsed + 1.0 /= 1.0 then
          speedup := seri_elapsed/mult_elapsed;
          efficiency := speedup/duration(nbt);
          efficiency := duration(100)*efficiency;
          if verbose then
            put("The speedup : "); duration_io.put(speedup,1,3);
            put("  efficiency : "); duration_io.put(efficiency,1,3); new_line;
          end if;
          duration_io.put(file,speedup,1,3); put(file," : ");
          duration_io.put(file,efficiency,2,2); new_line(file); flush(file);
        end if;
        Standard_Complex_VecVecs.Clear(values);
      end loop;
    end if;
  end Standard_Benchmark;

  procedure Standard_Coefficient_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in Standard_Coefficient_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Runs a benchmark test in double precision,
  --   on coefficient circuits.

  -- REQUIRED : if nbruns = 0, then nbtseq /= null.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbruns   the number of multitasked runs,
  --            if zero, then nbtseq will be used;
  --   inc      increment on the number of tasks, if nbruns /= 0;
  --   nbtseq   sequence of number of tasks for multitasked runs;
  --   s        system in one parameter;
  --   x        some point to evaluate at;
  --   verbose  if extra output is needed.

    dim : constant integer32 := s.dim;
    xr,xi : Standard_Floating_Vectors.Link_to_Vector;
    U,V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    e : Standard_Complex_Vectors.Vector(1..dim);
    svl : constant Standard_Complex_VecVecs.VecVec(0..s.neq)
        := Standard_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant Standard_Complex_VecMats.VecMat(1..s.neq)
       := Standard_Complex_Circuits.Allocate(s.neq,s.dim);
    values : Standard_Complex_VecVecs.VecVec(0..s.neq);
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    if verbose
     then put_line("Computing first without multitasking ...");
    end if;
    put_line(file,"double precision on coefficient circuits"); flush(file);
    Standard_Vector_Splitters.Split_Complex(x,xr,xi);
    seristart := Ada.Calendar.Clock;
    Standard_Coefficient_Circuits.Singular_Values(s,xr,xi,vh,U,V,e,svl);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    if verbose then
      put_line("-> Elapsed time without multitasking :");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      put_line("running in double precision ...");
    end if;
    put(file,"  1 : "); duration_io.put(file,seri_elapsed,1,3);
    put_line(file," : 1.000"); flush(file);
    if nbruns /= 0 then
      for k in 1..nbruns loop
        values := Standard_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
        multstart := Ada.Calendar.Clock;
        Multitasked_Hessian_Circuits.Multitasked_Singular_Values
          (nbt,s,xr,xi,values,false);
        multstop := Ada.Calendar.Clock;
        mult_elapsed := multstop - multstart;
        if verbose then
          put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
          Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
        end if;
        put(file,nbt,3); put(file," : ");
        duration_io.put(file,mult_elapsed,1,3); put(file," : ");
        if seri_elapsed + 1.0 /= 1.0 then
          speedup := seri_elapsed/mult_elapsed;
          efficiency := speedup/duration(nbt);
          efficiency := duration(100)*efficiency;
          if verbose then
            put("The speedup : "); duration_io.put(speedup,1,3);
            put("  efficiency : "); duration_io.put(efficiency,1,3); new_line;
          end if;
          duration_io.put(file,speedup,1,3); put(file," : ");
          duration_io.put(file,efficiency,2,2); new_line(file); flush(file);
        end if;
        nbt := nbt + inc;
        Standard_Complex_VecVecs.Clear(values);
      end loop;
    else
      for k in nbtseq'range loop
        nbt := nbtseq(k);
        values := Standard_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
        multstart := Ada.Calendar.Clock;
        Multitasked_Hessian_Circuits.Multitasked_Singular_Values
          (nbt,s,xr,xi,values,false);
        multstop := Ada.Calendar.Clock;
        mult_elapsed := multstop - multstart;
        if verbose then
          put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
          Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
        end if;
        put(file,nbt,3); put(file," : ");
        duration_io.put(file,mult_elapsed,1,3); put(file," : ");
        if seri_elapsed + 1.0 /= 1.0 then
          speedup := seri_elapsed/mult_elapsed;
          efficiency := speedup/duration(nbt);
          efficiency := duration(100)*efficiency;
          if verbose then
            put("The speedup : "); duration_io.put(speedup,1,3);
            put("  efficiency : "); duration_io.put(efficiency,1,3); new_line;
          end if;
          duration_io.put(file,speedup,1,3); put(file," : ");
          duration_io.put(file,efficiency,2,2); new_line(file); flush(file);
        end if;
        Standard_Complex_VecVecs.Clear(values);
      end loop;
    end if;
  end Standard_Coefficient_Benchmark;

  procedure DoblDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in DoblDobl_Complex_Circuits.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Runs a benchmark test in double double precision.

  -- REQUIRED : if nbruns = 0, then nbtseq /= null.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbruns   the number of multitasked runs,
  --            if zero, then nbtseq will be used;
  --   inc      increment on the number of tasks, if nbruns /= 0;
  --   nbtseq   sequence of numbers for multitasked runs;
  --   s        system in one parameter;
  --   x        some point to evaluate at;
  --   verbose  if extra output is needed.

    dim : constant integer32 := s.dim;
    U,V : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : DoblDobl_Complex_Vectors.Vector(1..dim);
    svl : constant DoblDobl_Complex_VecVecs.VecVec(0..s.neq)
        := DoblDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant DoblDobl_Complex_VecMats.VecMat(1..s.neq)
       := DoblDobl_Complex_Circuits.Allocate(s.neq,s.dim);
    values : DoblDobl_Complex_VecVecs.VecVec(0..s.neq);
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    if verbose
     then put_line("Computing first without multitasking ...");
    end if;
    put_line(file,"double double precision"); flush(file);
    seristart := Ada.Calendar.Clock;
    DoblDobl_Complex_Circuits.Singular_Values(s,x,vh,U,V,e,svl);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    if verbose then
      put_line("-> Elapsed time without multitasking :");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      put_line("running in double double precision ...");
    end if;
    put(file,"  1 : "); duration_io.put(file,seri_elapsed,1,3);
    put_line(file," : 1.000"); flush(file);
    if nbruns /= 0 then
      for k in 1..nbruns loop
        values := DoblDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
        multstart := Ada.Calendar.Clock;
        Multitasked_Hessian_Circuits.Multitasked_Singular_Values
          (nbt,s,x,values,false,false);
        multstop := Ada.Calendar.Clock;
        mult_elapsed := multstop - multstart;
        if verbose then
          put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
          Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
        end if;
        put(file,nbt,3); put(file," : ");
        duration_io.put(file,mult_elapsed,1,3); put(file," : ");
        if seri_elapsed + 1.0 /= 1.0 then
          speedup := seri_elapsed/mult_elapsed;
          efficiency := speedup/duration(nbt);
          efficiency := duration(100)*efficiency;
          if verbose then
            put("The speedup : "); duration_io.put(speedup,1,3); 
            put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
          end if;
          duration_io.put(file,speedup,1,3); put(file," : ");
          duration_io.put(file,efficiency,2,2); new_line(file); flush(file);
        end if;
        nbt := nbt + inc;
        DoblDobl_Complex_VecVecs.Clear(values);
      end loop;
    else
      for k in nbtseq'range loop
        nbt := nbtseq(k);
        values := DoblDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
        multstart := Ada.Calendar.Clock;
        Multitasked_Hessian_Circuits.Multitasked_Singular_Values
          (nbt,s,x,values,false,false);
        multstop := Ada.Calendar.Clock;
        mult_elapsed := multstop - multstart;
        if verbose then
          put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
          Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
        end if;
        put(file,nbt,3); put(file," : ");
        duration_io.put(file,mult_elapsed,1,3); put(file," : ");
        if seri_elapsed + 1.0 /= 1.0 then
          speedup := seri_elapsed/mult_elapsed;
          efficiency := speedup/duration(nbt);
          efficiency := duration(100)*efficiency;
          if verbose then
            put("The speedup : "); duration_io.put(speedup,1,3); 
            put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
          end if;
          duration_io.put(file,speedup,1,3); put(file," : ");
          duration_io.put(file,efficiency,2,2); new_line(file); flush(file);
        end if;
        DoblDobl_Complex_VecVecs.Clear(values);
      end loop;
    end if;
  end DoblDobl_Benchmark;

  procedure QuadDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in QuadDobl_Complex_Circuits.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Runs a benchmark test in quad double precision.

  -- REQUIRED : if nbruns = 0, then nbtseq /= null.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   nbruns   the number of multitasked runs,
  --            if zero, then nbtseq will be used;
  --   inc      increment to the number of tasks, if nbruns /= 0;
  --   nbtseq   sequence of number of tasks for multitasked runs;
  --   s        system in one parameter;
  --   x        some point to evaluate at;
  --   verbose  if extra output needs to be written to screen.

    dim : constant integer32 := s.dim;
    U,V : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : QuadDobl_Complex_Vectors.Vector(1..dim);
    svl : constant QuadDobl_Complex_VecVecs.VecVec(0..s.neq)
        := QuadDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant QuadDobl_Complex_VecMats.VecMat(1..s.neq)
       := QuadDobl_Complex_Circuits.Allocate(s.neq,s.dim);
    values : QuadDobl_Complex_VecVecs.VecVec(0..s.neq);
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    if verbose
     then put_line("Computing first without multitasking ...");
    end if;
    put_line(file,"quad double precision"); flush(file);
    seristart := Ada.Calendar.Clock;
    QuadDobl_Complex_Circuits.Singular_Values(s,x,vh,U,V,e,svl);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    if verbose then
      put_line("-> Elapsed time without multitasking :");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      put_line("running in quad double precision ...");
    end if;
    put(file,"  1 : "); duration_io.put(file,seri_elapsed,1,3);
    put_line(file," : 1.000"); flush(file);
    if nbruns /= 0 then
      for k in 1..nbruns loop
        values := QuadDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
        multstart := Ada.Calendar.Clock;
        Multitasked_Hessian_Circuits.Multitasked_Singular_Values
          (nbt,s,x,values,false,false);
        multstop := Ada.Calendar.Clock;
        mult_elapsed := multstop - multstart;
        if verbose then
          put("-> Elapsed time with "); put(file,nbt,1); put_line(" tasks :");
          Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
        end if;
        put(file,nbt,3); put(file," : ");
        duration_io.put(file,mult_elapsed,1,3); put(file," : ");
        if seri_elapsed + 1.0 /= 1.0 then
          speedup := seri_elapsed/mult_elapsed;
          efficiency := speedup/duration(nbt);
          efficiency := duration(100)*efficiency;
          if verbose then
            put("The speedup : "); duration_io.put(speedup,1,3);
            put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
          end if;
          duration_io.put(file,speedup,1,3); put(file," : ");
          duration_io.put(file,efficiency,2,2); new_line(file); flush(file);
        end if;
        nbt := nbt + inc;
        QuadDobl_Complex_VecVecs.Clear(values);
      end loop;
    else
      for k in nbtseq'range loop
        nbt := nbtseq(k);
        values := QuadDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
        multstart := Ada.Calendar.Clock;
        Multitasked_Hessian_Circuits.Multitasked_Singular_Values
          (nbt,s,x,values,false,false);
        multstop := Ada.Calendar.Clock;
        mult_elapsed := multstop - multstart;
        if verbose then
          put("-> Elapsed time with "); put(file,nbt,1); put_line(" tasks :");
          Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
        end if;
        put(file,nbt,3); put(file," : ");
        duration_io.put(file,mult_elapsed,1,3); put(file," : ");
        if seri_elapsed + 1.0 /= 1.0 then
          speedup := seri_elapsed/mult_elapsed;
          efficiency := speedup/duration(nbt);
          efficiency := duration(100)*efficiency;
          if verbose then
            put("The speedup : "); duration_io.put(speedup,1,3);
            put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
          end if;
          duration_io.put(file,speedup,1,3); put(file," : ");
          duration_io.put(file,efficiency,2,2); new_line(file); flush(file);
        end if;
        QuadDobl_Complex_VecVecs.Clear(values);
      end loop;
    end if;
  end QuadDobl_Benchmark;

  function Prompt_for_Sequence
             ( max : in integer32 )
             return Standard_Integer_Vectors.Link_to_Vector is

  -- DESCRIPTION :
  --   Prompts the user for a sequence of numbers
  --   and the sequence is returned in a vector.
  --   The length of the sequence may not exceed max.

    res : Standard_Integer_Vectors.Link_to_Vector;
    seq : Standard_Integer_Vectors.Vector(1..max);
    cnt,nbr : integer32 := 0;
    ans : character;

  begin
    loop
      put_line("Reading a sequence of numbers ...");
      for k in 1..max loop
        put("-> give a number (0 to exit) : "); get(nbr);
        exit when (nbr = 0);
        cnt := cnt + 1;
        seq(cnt) := nbr;
      end loop;
      put("The sequence : "); put(seq(1..cnt)); new_line;
      put("Okay to exit ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans = 'y');
      put_line("Entering the sequence again ...");
      cnt := 0;
    end loop;
    res := new Standard_Integer_Vectors.Vector'(seq(1..cnt));
    return res;
  end Prompt_For_Sequence;

  procedure Benchmark ( dim,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random system of circuits in quad double precision,
  --   and runs benchmark tests in all three precisions.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    qds : QuadDobl_Complex_Circuits.Link_to_System;
    dds : DoblDobl_Complex_Circuits.Link_to_System;
    d_s : Standard_Complex_Circuits.Link_to_System;
    dcs : Standard_Coefficient_Circuits.Link_to_System;
    qdx : constant QuadDobl_Complex_Vectors.Vector(1..dim)
        := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    vqdx : constant QuadDobl_Complex_Vectors.Link_to_Vector
         := new QuadDobl_Complex_Vectors.Vector'(qdx);
    ddx : constant DoblDobl_Complex_Vectors.Vector(1..dim)
        := QuadDobl_Complex_Vectors_cv.QuadDobl_Complex_to_DoblDobl(qdx);
    vddx : constant DoblDobl_Complex_Vectors.Link_to_Vector
         := new DoblDobl_Complex_Vectors.Vector'(ddx);
    dx : constant Standard_Complex_Vectors.Vector(1..dim)
        := QuadDobl_Complex_Vectors_cv.QuadDobl_Complex_to_Standard(qdx);
    vdx : constant Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector'(dx);
    file : file_type;
    nbruns,inc : integer32 := 0;
    nbtseq : Standard_Integer_Vectors.Link_to_Vector;
 
  begin
    new_line;
    put("Give the number of multitasked runs (0 for sequence) : ");
    get(nbruns);
    if nbruns /= 0
     then put("Give the increment on the tasks : "); get(inc); skip_line;
     else nbtseq := Prompt_for_Sequence(44);
    end if;
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    qds := QuadDobl_Circuit_Makers.Random_Complex_System(dim,nbr,dim,pwr);
    dds := DoblDobl_Circuit_Makers.to_double_double(qds);
    d_s := Standard_Circuit_Makers.to_double(qds);
    dcs := Standard_Circuit_Splitters.Split(d_s);
    Standard_Benchmark(file,nbruns,inc,nbtseq,d_s,vdx);
    Standard_Coefficient_Benchmark(file,nbruns,inc,nbtseq,dcs,vdx);
    DoblDobl_Benchmark(file,nbruns,inc,nbtseq,dds,vddx);
    QuadDobl_Benchmark(file,nbruns,inc,nbtseq,qds,vqdx);
  end Benchmark;

  procedure Benchmark
              ( p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   For the given polynomial system p with some solutions in sols,
  --   prompts the user for the parameters and then runs benchmarks.

    file : file_type;
    qds : QuadDobl_Complex_Circuits.Link_to_System;
    dds : DoblDobl_Complex_Circuits.Link_to_System;
    d_s : Standard_Complex_Circuits.Link_to_System;
    dcs : Standard_Coefficient_Circuits.Link_to_System;
    ls : constant QuadDobl_Complex_Solutions.Link_to_Solution
       := QuadDobl_Complex_Solutions.Head_Of(sols);
    qdx : constant QuadDobl_Complex_Vectors.Link_to_Vector
        := new QuadDobl_Complex_Vectors.Vector'(ls.v);
    dim : constant integer32 := ls.n;
    ddx : constant DoblDobl_Complex_Vectors.Vector(1..dim)
        := QuadDobl_Complex_Vectors_cv.QuadDobl_Complex_to_DoblDobl(ls.v);
    vddx : constant DoblDobl_Complex_Vectors.Link_to_Vector
         := new DoblDobl_Complex_Vectors.Vector'(ddx);
    dx : constant Standard_Complex_Vectors.Vector(1..dim)
        := QuadDobl_Complex_Vectors_cv.QuadDobl_Complex_to_Standard(ls.v);
    vdx : constant Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector'(dx);
    nbruns,inc : integer32 := 0;
    nbtseq : Standard_Integer_Vectors.Link_to_Vector;
 
  begin
    new_line;
    put_line("making circuits for the given system ...");
    qds := QuadDobl_Circuit_Makers.Make_Complex_System(p);
    new_line;
    put("Give the number of multitasked runs (0 for sequence) : ");
    get(nbruns);
    if nbruns /= 0
     then put("Give the increment on the tasks : "); get(inc); skip_line;
     else nbtseq := Prompt_for_Sequence(44);
    end if;
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    dds := DoblDobl_Circuit_Makers.to_double_double(qds);
    d_s := Standard_Circuit_Makers.to_double(qds);
    dcs := Standard_Circuit_Splitters.Split(d_s);
    Standard_Benchmark(file,nbruns,inc,nbtseq,d_s,vdx);
    Standard_Coefficient_Benchmark(file,nbruns,inc,nbtseq,dcs,vdx);
    DoblDobl_Benchmark(file,nbruns,inc,nbtseq,dds,vddx);
    QuadDobl_Benchmark(file,nbruns,inc,nbtseq,qds,qdx);
  end Benchmark;

  function Prompt_for_Precision return character is

  -- DESCRIPTION :
  --   Prompts the user for the work precision
  --   and returns '0', '1', or '2',
  --   for double, double double, or quad double precision.

    precision : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(precision,"012");
    return precision;
  end Prompt_for_Precision;

  procedure Main is

  -- DESCRIPTION :
  --   Launches the test after prompting for the parameters
  --   to generate a random problem.

    dim,nbr,pwr : integer32 := 0;
    ans,precision : character;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Testing the Hessian criterion ...");
    new_line;
    put("Generate a random system ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'n' then
      new_line;
      put("Benchmark for all precisions ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("Reading a polynomial system with solutions ...");
        QuadDobl_System_and_Solutions_io.get(lp,sols);
        Benchmark(lp,sols);
      else
        precision := Prompt_for_Precision;
        case precision is
          when '0' => Standard_User_Test;
          when '1' => DoblDobl_User_Test;
          when '2' => QuadDobl_User_Test;
          when others => null;
        end case;
      end if;
    else
      put("Give the dimension : "); get(dim);
      put("Give the number of terms in each circuit : "); get(nbr);
      put("Give the largest power of the variables : "); get(pwr);
      new_line;
      put("Benchmark for all precisions ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y' then
        Benchmark(dim,nbr,pwr);
      else
        precision := Prompt_for_Precision;
        case precision is
          when '0' => Standard_Random_Test(dim,nbr,pwr);
          when '1' => DoblDobl_Random_Test(dim,nbr,pwr);
          when '2' => QuadDobl_Random_Test(dim,nbr,pwr);
          when others => null;
        end case;
      end if;
    end if;
  end Main;

begin
  Main;
end ts_mthesscrc;
