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
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors_cv;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with System_Convolution_Circuits;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Evaluation_Differentiation_Errors;  use Evaluation_Differentiation_Errors;
with Hessian_Convolution_Circuits;       use Hessian_Convolution_Circuits;
with Jacobian_Convolution_Circuits;      use Jacobian_Convolution_Circuits;
with Multitasked_Hessian_Convolutions;   use Multitasked_Hessian_Convolutions;

procedure ts_mthessian is

-- DESCRIPTION :
--   Tests the development of the Hessian criterion,
--   in double, double double, and quad double arithmetic,
--   with multitasking for shared memory parallel computers.

  procedure Write_Singular_Values 
              ( values : in Standard_Complex_VecVecs.VecVec;
                jmvals : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Writes the singular values in values, line by line.
  --   The vector jmvals contains the singular values of the Jacobian.

    lnk : Standard_Complex_Vectors.Link_to_Vector;
    val : double_float;

  begin
    for k in values'range loop
      lnk := values(k);
      for i in lnk'range loop
        val := Standard_Complex_Numbers.REAL_PART(lnk(i));
        put(val,3);
      end loop;
      new_line;
    end loop;
    for i in jmvals'range loop
      val := Standard_Complex_Numbers.REAL_PART(jmvals(i));
      put(val,3);
    end loop;
    new_line;
  end Write_Singular_Values;

  procedure Write_Singular_Values 
              ( values : in DoblDobl_Complex_VecVecs.VecVec;
                jmvals : in DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Writes the singular values in values, line by line.
  --   The vector jmvals contains the singular values of the Jacobian.

    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    val : double_double;

  begin
    for k in values'range loop
      lnk := values(k);
      for i in lnk'range loop
        val := DoblDobl_Complex_Numbers.REAL_PART(lnk(i));
        put(" "); put(val,3);
      end loop;
      new_line;
    end loop;
    for i in jmvals'range loop
      val := DoblDobl_Complex_Numbers.REAL_PART(jmvals(i));
      put(" "); put(val,3);
    end loop;
    new_line;
  end Write_Singular_Values;

  procedure Write_Singular_Values 
              ( values : in QuadDobl_Complex_VecVecs.VecVec;
                jmvals : in QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Writes the singular values in values, line by line.
  --   The vector jmvals contains the singular values of the Jacobian.

    lnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    val : quad_double;

  begin
    for k in values'range loop
      lnk := values(k);
      for i in lnk'range loop
        val := QuadDobl_Complex_Numbers.REAL_PART(lnk(i));
        put(" "); put(val,3);
      end loop;
      new_line;
    end loop;
    for i in jmvals'range loop
      val := QuadDobl_Complex_Numbers.REAL_PART(jmvals(i));
      put(" "); put(val,3);
    end loop;
    new_line;
  end Write_Singular_Values;

  procedure Standard_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use Standard_Complex_VecVecs;
    use Standard_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    vx : Standard_Complex_Vectors.Vector(1..dim);
    A,U,V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    e : Standard_Complex_Vectors.Vector(1..dim);
    svl : constant Standard_Complex_VecVecs.VecVec(1..dim) := Allocate(dim,dim);
    values : Standard_Complex_VecVecs.VecVec(1..dim) := Allocate(dim,dim);
    jm1vls,jm2vls : Standard_Complex_Vectors.Vector(1..dim);
    ans : character;
    otp : boolean;
    nbt : integer32 := 0;
    err,eta : double_float;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup : Duration;

    use Ada.Calendar;

  begin
    new_line;
    put("Extra output wanted ? (y/n) ");
    Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    Standard_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    for i in 1..dim loop
      lnk := x(i);
      vx(i) := lnk(0);
    end loop;
    put_line("Computing first without multitasking ...");
    seristart := Ada.Calendar.Clock;
    Singular_Values(s.crc,vx,A,U,V,e,jm1vls);
    for i in 1..dim loop
      lnk := svl(i);
      Singular_Values(s.crc(i),vx,A,U,V,e,lnk.all);
    end loop;
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking :");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    if otp
     then put_line("All singular values :"); Write_Singular_Values(svl,jm1vls);
    end if;
    loop
      new_line;
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt <= 0);
      multstart := Ada.Calendar.Clock;
      Multitasked_Singular_Values(nbt,s,vx,jm2vls,values,otp);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if otp then
        put_line("All singular values computed by multitasking :");
        Write_Singular_Values(values,jm2vls);
      end if;
      err := Difference(svl,values);
      eta := Standard_Distance(jm2vls,values);
      put("The difference error : "); put(err,3); new_line;
      put("Estimate for the distance : "); put(eta,3); new_line;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        put("The speedup : ");
        duration_io.put(speedup,1,3); new_line;
      end if;
    end loop;
  end Standard_Random_Test;

  procedure DoblDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in double double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    vx : DoblDobl_Complex_Vectors.Vector(1..dim);
    A,U,V : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : DoblDobl_Complex_Vectors.Vector(1..dim);
    svl : constant DoblDobl_Complex_VecVecs.VecVec(1..dim) := Allocate(dim,dim);
    values : DoblDobl_Complex_VecVecs.VecVec(1..dim) := Allocate(dim,dim);
    jm1vls,jm2vls : DoblDobl_Complex_Vectors.Vector(1..dim);
    ans : character;
    otp : boolean;
    nbt : integer32 := 0;
    err,eta : double_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup : Duration;

    use Ada.Calendar;

  begin
    new_line;
    put("Extra output wanted ? (y/n) ");
    Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    DoblDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    for i in 1..dim loop
      lnk := x(i);
      vx(i) := lnk(0);
    end loop;
    put_line("Computing first without multitasking ...");
    seristart := Ada.Calendar.Clock;
    Singular_Values(s.crc,vx,A,U,V,e,jm1vls);
    for i in 1..dim loop
      lnk := svl(i);
      Singular_Values(s.crc(i),vx,A,U,V,e,lnk.all);
    end loop;
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking :");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    if otp
     then put_line("All singular values :"); Write_Singular_Values(svl,jm1vls);
    end if;
    loop
      new_line;
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt <= 0);
      multstart := Ada.Calendar.Clock;
      Multitasked_Singular_Values(nbt,s,vx,jm2vls,values,otp);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if otp then
        put_line("All singular values computed by multitasking :");
        Write_Singular_Values(values,jm2vls);
      end if;
      err := Difference(svl,values);
      eta := DoblDobl_Distance(jm2vls,values);
      put("The difference error : "); put(err,3); new_line;
      put("Estimate for the distance : "); put(eta,3); new_line;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        put("The speedup : ");
        duration_io.put(speedup,1,3); new_line;
      end if;
    end loop;
  end DoblDobl_Random_Test;

  procedure QuadDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in quad double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use QuadDobl_Complex_VecVecs;
    use QuadDobl_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    lnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    vx : QuadDobl_Complex_Vectors.Vector(1..dim);
    A,U,V : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : QuadDobl_Complex_Vectors.Vector(1..dim);
    svl : constant QuadDobl_Complex_VecVecs.VecVec(1..dim) := Allocate(dim,dim);
    values : QuadDobl_Complex_VecVecs.VecVec(1..dim) := Allocate(dim,dim);
    jm1vls,jm2vls : QuadDobl_Complex_Vectors.Vector(1..dim);
    ans : character;
    otp : boolean;
    nbt : integer32 := 0;
    err,eta : quad_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup : Duration;

    use Ada.Calendar;

  begin
    new_line;
    put("Extra output wanted ? (y/n) ");
    Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    QuadDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    for i in 1..dim loop
      lnk := x(i);
      vx(i) := lnk(0);
    end loop;
    put_line("Computing first without multitasking ...");
    seristart := Ada.Calendar.Clock;
    Singular_Values(s.crc,vx,A,U,V,e,jm1vls);
    for i in 1..dim loop
      lnk := svl(i);
      Singular_Values(s.crc(i),vx,A,U,V,e,lnk.all);
    end loop;
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put_line("-> Elapsed time without multitasking :");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    if otp
     then put_line("All singular values :"); Write_Singular_Values(svl,jm1vls);
    end if;
    loop
      new_line;
      put("Give the number of tasks (0 to exit) : "); get(nbt);
      exit when (nbt <= 0);
      multstart := Ada.Calendar.Clock;
      Multitasked_Singular_Values(nbt,s,vx,jm2vls,values,otp);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
      Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
      if otp then
        put_line("All singular values computed by multitasking :");
        Write_Singular_Values(values,jm2vls);
      end if;
      err := Difference(svl,values);
      eta := QuadDobl_Distance(jm2vls,values);
      put("The difference error : "); put(err,3); new_line;
      put("Estimate for the distance : "); put(eta,3); new_line;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        put("The speedup : ");
        duration_io.put(speedup,1,3); new_line;
      end if;
    end loop;
  end QuadDobl_Random_Test;

  function to_double_double
             ( v : QuadDobl_Complex_VecVecs.Link_to_VecVec )
             return DoblDobl_Complex_VecVecs.Link_to_VecVec is

  -- DESCRIPTION :
  --   Returns the vector v converted to double double precision.

  -- REQUIRED : v /= null;

    res : DoblDobl_Complex_VecVecs.Link_to_VecVec;
    ddv : DoblDobl_Complex_VecVecs.VecVec(v'range);

    use QuadDobl_Complex_Vectors_cv;

  begin
    for i in ddv'range loop
      declare
        vec : constant DoblDobl_Complex_Vectors.Vector
            := QuadDobl_Complex_to_DoblDobl(v(i).all);
      begin
        ddv(i) := new DoblDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    res := new DoblDobl_Complex_VecVecs.VecVec'(ddv);
    return res;
  end to_double_double;

  function to_double
             ( v : QuadDobl_Complex_VecVecs.Link_to_VecVec )
             return Standard_Complex_VecVecs.Link_to_VecVec is

  -- DESCRIPTION :
  --   Returns the vector v converted to double precision.

  -- REQUIRED : v /= null;

    res : Standard_Complex_VecVecs.Link_to_VecVec;
    ddv : Standard_Complex_VecVecs.VecVec(v'range);

    use QuadDobl_Complex_Vectors_cv;

  begin
    for i in ddv'range loop
      declare
        vec : constant Standard_Complex_Vectors.Vector
            := QuadDobl_Complex_to_Standard(v(i).all);
      begin
        ddv(i) := new Standard_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    res := new Standard_Complex_VecVecs.VecVec'(ddv);
    return res;
  end to_double;

  procedure Benchmark ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in quad double precision,
  --   and runs benchmark tests in all three precisions.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    qds : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    dds : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    d_s : Standard_Speelpenning_Convolutions.Link_to_System;
    qdx : QuadDobl_Complex_VecVecs.Link_to_VecVec;
    ddx : DoblDobl_Complex_VecVecs.Link_to_VecVec;
    d_x : Standard_Complex_VecVecs.Link_to_VecVec;
    qd_values : QuadDobl_Complex_VecVecs.VecVec(1..dim) := Allocate(dim,dim);
    dd_values : DoblDobl_Complex_VecVecs.VecVec(1..dim) := Allocate(dim,dim);
    d_values : Standard_Complex_VecVecs.VecVec(1..dim) := Allocate(dim,dim);
    qd_jmvls : QuadDobl_Complex_Vectors.Vector(1..dim);
    dd_jmvls : DoblDobl_Complex_Vectors.Vector(1..dim);
    d_jmvls : Standard_Complex_Vectors.Vector(1..dim);
    qdlnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    ddlnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    dlnk : Standard_Complex_Vectors.Link_to_Vector;
    qdvx : QuadDobl_Complex_Vectors.Vector(1..dim);
    ddvx : DoblDobl_Complex_Vectors.Vector(1..dim);
    dvx : Standard_Complex_Vectors.Vector(1..dim);
    nbt : integer32 := 0;
    start,stop : Ada.Calendar.Time;

  begin
    new_line;
    put("Give the number of tasks : "); get(nbt);
    QuadDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,qds,qdx);
    dds := System_Convolution_Circuits.to_double_double(qds);
    ddx := to_double_double(qdx);
    d_s := System_Convolution_Circuits.to_double(qds);
    d_x := to_double(qdx);
    for i in 1..dim loop
      qdlnk := qdx(i); ddlnk := ddx(i); dlnk := d_x(i);
      qdvx(i) := qdlnk(0); ddvx(i) := ddlnk(0); dvx(i) := dlnk(0);
    end loop;
    put_line("running in double precision ...");
    start := Ada.Calendar.Clock;
    Multitasked_Singular_Values(nbt,d_s,dvx,d_jmvls,d_values,false);
    stop := Ada.Calendar.Clock;
    put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
    Time_Stamps.Write_Elapsed_Time(standard_output,start,stop);
    put_line("running in double double precision ...");
    start := Ada.Calendar.Clock;
    Multitasked_Singular_Values(nbt,dds,ddvx,dd_jmvls,dd_values,false);
    stop := Ada.Calendar.Clock;
    put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
    Time_Stamps.Write_Elapsed_Time(standard_output,start,stop);
    put_line("running in quad double precision ...");
    start := Ada.Calendar.Clock;
    Multitasked_Singular_Values(nbt,qds,qdvx,qd_jmvls,qd_values,false);
    stop := Ada.Calendar.Clock;
    put("-> Elapsed time with "); put(nbt,1); put_line(" tasks :");
    Time_Stamps.Write_Elapsed_Time(standard_output,start,stop);
  end Benchmark;

  procedure Main is

  -- DESCRIPTION :
  --   Launches the test after prompting for the parameters
  --   to generate a random problem.

    dim,deg,nbr,pwr : integer32 := 0;
    ans,precision : character;

  begin
    new_line;
    put_line("Testing the Hessian criterion ...");
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree of the power series : "); get(deg);
    put("Give the number of terms in each circuit : "); get(nbr);
    put("Give the largest power of the variables : "); get(pwr);
    new_line;
    put("Benchmark for all precisions ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Benchmark(dim,deg,nbr,pwr);
    else
      new_line;
      put_line("MENU for the working precision :");
      put_line("  0. double precision");
      put_line("  1. double double precision");
      put_line("  2. quad double precision");
      put("Type 0, 1, or 2 to select the precision : ");
      Ask_Alternative(precision,"012");
      case precision is
        when '0' => Standard_Random_Test(dim,deg,nbr,pwr);
        when '1' => DoblDobl_Random_Test(dim,deg,nbr,pwr);
        when '2' => QuadDobl_Random_Test(dim,deg,nbr,pwr);
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_mthessian;
