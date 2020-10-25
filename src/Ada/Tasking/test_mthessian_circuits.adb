with duration_io;
with Ada.Calendar;
with Time_Stamps;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Triple_Double_Numbers_io;           use Triple_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Penta_Double_Numbers_io;            use Penta_Double_Numbers_io;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Octo_Double_Numbers_io;             use Octo_Double_Numbers_io;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Deca_Double_Numbers_io;             use Deca_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with TripDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with PentDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Random_Vectors;
with DoblDobl_Random_Vectors;
with TripDobl_Random_Vectors;
with QuadDobl_Random_Vectors;
with PentDobl_Random_Vectors;
with OctoDobl_Random_Vectors;
with DecaDobl_Random_Vectors;
with Standard_Vector_Splitters;
with DoblDobl_Vector_Splitters;
with TripDobl_Vector_Splitters;
with QuadDobl_Vector_Splitters;
with PentDobl_Vector_Splitters;
with OctoDobl_Vector_Splitters;
with DecaDobl_Vector_Splitters;
with DecaDobl_Complex_Vectors_cv;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with TripDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with PentDobl_Complex_Poly_Systems;
with OctoDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_System_and_Solutions_io;
with TripDobl_Complex_Solutions;
with TripDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_System_and_Solutions_io;
with PentDobl_Complex_Solutions;
with PentDobl_System_and_Solutions_io;
with OctoDobl_Complex_Solutions;
with OctoDobl_System_and_Solutions_io;
with DecaDobl_System_and_Solutions_io;
with Evaluation_Differentiation_Errors;  use Evaluation_Differentiation_Errors;
with Standard_Circuit_Splitters;
with Standard_Circuit_Makers;
with DoblDobl_Circuit_Makers;
with TripDobl_Circuit_Makers;
with QuadDobl_Circuit_Makers;
with PentDobl_Circuit_Makers;
with OctoDobl_Circuit_Makers;
with DecaDobl_Circuit_Makers;
with Multitasked_Hessian_Circuits;

package body Test_mtHessian_Circuits is

  procedure Write_Singular_Values 
              ( values : in Standard_Complex_VecVecs.VecVec ) is

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
              ( values : in TripDobl_Complex_VecVecs.VecVec ) is

    lnk : TripDobl_Complex_Vectors.Link_to_Vector;
    val : triple_double;

  begin
    for k in values'range loop
      lnk := values(k);
      for i in lnk'first..lnk'last-1 loop
        val := TripDobl_Complex_Numbers.REAL_PART(lnk(i));
        put(" "); put(val,3);
      end loop;
      new_line;
    end loop;
  end Write_Singular_Values;

  procedure Write_Singular_Values 
              ( values : in QuadDobl_Complex_VecVecs.VecVec ) is

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

  procedure Write_Singular_Values 
              ( values : in PentDobl_Complex_VecVecs.VecVec ) is

    lnk : PentDobl_Complex_Vectors.Link_to_Vector;
    val : penta_double;

  begin
    for k in values'range loop
      lnk := values(k);
      for i in lnk'first..lnk'last-1 loop
        val := PentDobl_Complex_Numbers.REAL_PART(lnk(i));
        put(" "); put(val,3);
      end loop;
      new_line;
    end loop;
  end Write_Singular_Values;

  procedure Write_Singular_Values 
              ( values : in OctoDobl_Complex_VecVecs.VecVec ) is

    lnk : OctoDobl_Complex_Vectors.Link_to_Vector;
    val : octo_double;

  begin
    for k in values'range loop
      lnk := values(k);
      for i in lnk'first..lnk'last-1 loop
        val := OctoDobl_Complex_Numbers.REAL_PART(lnk(i));
        put(" "); put(val,3);
      end loop;
      new_line;
    end loop;
  end Write_Singular_Values;

  procedure Write_Singular_Values 
              ( values : in DecaDobl_Complex_VecVecs.VecVec ) is

    lnk : DecaDobl_Complex_Vectors.Link_to_Vector;
    val : deca_double;

  begin
    for k in values'range loop
      lnk := values(k);
      for i in lnk'first..lnk'last-1 loop
        val := DecaDobl_Complex_Numbers.REAL_PART(lnk(i));
        put(" "); put(val,3);
      end loop;
      new_line;
    end loop;
  end Write_Singular_Values;

  procedure Standard_Test
              ( s : in Standard_Complex_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                static,output : in boolean ) is

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

  procedure TripDobl_Test
              ( s : in TripDobl_Complex_Circuits.Link_to_System;
                x : in TripDobl_Complex_Vectors.Link_to_Vector;
                static,output : in boolean ) is

    dim : constant integer32 := s.dim;
    U,V : TripDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : TripDobl_Complex_Vectors.Vector(1..dim);
    svl : constant TripDobl_Complex_VecVecs.VecVec(0..s.neq)
        := TripDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant TripDobl_Complex_VecMats.VecMat(1..s.neq)
       := TripDobl_Complex_Circuits.Allocate(s.neq,s.dim);
    values : TripDobl_Complex_VecVecs.VecVec(0..s.neq)
           := TripDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    nbt : integer32 := 0;
    err,eta : triple_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;

    use Ada.Calendar;

  begin
    put_line("Computing first without multitasking ...");
    seristart := Ada.Calendar.Clock;
    TripDobl_Complex_Circuits.Singular_Values(s,x,vh,U,V,e,svl);
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
      eta := Multitasked_Hessian_Circuits.TripDobl_Distance(svl);
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
  end TripDobl_Test;

  procedure QuadDobl_Test
              ( s : in QuadDobl_Complex_Circuits.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                static,output : in boolean ) is

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

  procedure PentDobl_Test
              ( s : in PentDobl_Complex_Circuits.Link_to_System;
                x : in PentDobl_Complex_Vectors.Link_to_Vector;
                static,output : in boolean ) is

    dim : constant integer32 := s.dim;
    U,V : PentDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : PentDobl_Complex_Vectors.Vector(1..dim);
    svl : constant PentDobl_Complex_VecVecs.VecVec(0..s.neq)
        := PentDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant PentDobl_Complex_VecMats.VecMat(1..s.neq)
       := PentDobl_Complex_Circuits.Allocate(s.neq,s.dim);
    values : PentDobl_Complex_VecVecs.VecVec(0..s.neq)
           := PentDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    nbt : integer32 := 0;
    err,eta : penta_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;

    use Ada.Calendar;

  begin
    put_line("Computing first without multitasking ...");
    seristart := Ada.Calendar.Clock;
    PentDobl_Complex_Circuits.Singular_Values(s,x,vh,U,V,e,svl);
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
      eta := Multitasked_Hessian_Circuits.PentDobl_Distance(svl);
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
  end PentDobl_Test;

  procedure OctoDobl_Test
              ( s : in OctoDobl_Complex_Circuits.Link_to_System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                static,output : in boolean ) is

    dim : constant integer32 := s.dim;
    U,V : OctoDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : OctoDobl_Complex_Vectors.Vector(1..dim);
    svl : constant OctoDobl_Complex_VecVecs.VecVec(0..s.neq)
        := OctoDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant OctoDobl_Complex_VecMats.VecMat(1..s.neq)
       := OctoDobl_Complex_Circuits.Allocate(s.neq,s.dim);
    values : OctoDobl_Complex_VecVecs.VecVec(0..s.neq)
           := OctoDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    nbt : integer32 := 0;
    err,eta : octo_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;

    use Ada.Calendar;

  begin
    put_line("Computing first without multitasking ...");
    seristart := Ada.Calendar.Clock;
    OctoDobl_Complex_Circuits.Singular_Values(s,x,vh,U,V,e,svl);
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
      eta := Multitasked_Hessian_Circuits.OctoDobl_Distance(svl);
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
  end OctoDobl_Test;

  procedure DecaDobl_Test
              ( s : in DecaDobl_Complex_Circuits.Link_to_System;
                x : in DecaDobl_Complex_Vectors.Link_to_Vector;
                static,output : in boolean ) is

    dim : constant integer32 := s.dim;
    U,V : DecaDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : DecaDobl_Complex_Vectors.Vector(1..dim);
    svl : constant DecaDobl_Complex_VecVecs.VecVec(0..s.neq)
        := DecaDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant DecaDobl_Complex_VecMats.VecMat(1..s.neq)
       := DecaDobl_Complex_Circuits.Allocate(s.neq,s.dim);
    values : DecaDobl_Complex_VecVecs.VecVec(0..s.neq)
           := DecaDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    nbt : integer32 := 0;
    err,eta : deca_double;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;

    use Ada.Calendar;

  begin
    put_line("Computing first without multitasking ...");
    seristart := Ada.Calendar.Clock;
    DecaDobl_Complex_Circuits.Singular_Values(s,x,vh,U,V,e,svl);
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
      eta := Multitasked_Hessian_Circuits.DecaDobl_Distance(svl);
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
  end DecaDobl_Test;

  procedure Standard_Random_Test
              ( dim,nbr,pwr : in integer32 ) is

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

  procedure TripDobl_Random_Test
              ( dim,nbr,pwr : in integer32 ) is

    s : TripDobl_Complex_Circuits.Link_to_System;
    v : constant TripDobl_Complex_Vectors.Vector(1..dim)
      := TripDobl_Random_Vectors.Random_Vector(1,dim);
    x : constant TripDobl_Complex_Vectors.Link_to_Vector
      := new TripDobl_Complex_Vectors.Vector'(v);
    ans : character;
    static,output : boolean;

  begin
    new_line;
    put_line("Generating a random complex circuit system ...");
    s := TripDobl_Circuit_Makers.Random_Complex_System(dim,nbr,dim,pwr);
    new_line;
    put("Extra output wanted ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
    static := (ans = 'y');
    TripDobl_Test(s,x,static,output);
  end TripDobl_Random_Test;

  procedure QuadDobl_Random_Test
              ( dim,nbr,pwr : in integer32 ) is

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

  procedure PentDobl_Random_Test
              ( dim,nbr,pwr : in integer32 ) is

    s : PentDobl_Complex_Circuits.Link_to_System;
    v : constant PentDobl_Complex_Vectors.Vector(1..dim)
      := PentDobl_Random_Vectors.Random_Vector(1,dim);
    x : constant PentDobl_Complex_Vectors.Link_to_Vector
      := new PentDobl_Complex_Vectors.Vector'(v);
    ans : character;
    static,output : boolean;

  begin
    new_line;
    put_line("Generating a random complex circuit system ...");
    s := PentDobl_Circuit_Makers.Random_Complex_System(dim,nbr,dim,pwr);
    new_line;
    put("Extra output wanted ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
    static := (ans = 'y');
    PentDobl_Test(s,x,static,output);
  end PentDobl_Random_Test;

  procedure OctoDobl_Random_Test
              ( dim,nbr,pwr : in integer32 ) is

    s : OctoDobl_Complex_Circuits.Link_to_System;
    v : constant OctoDobl_Complex_Vectors.Vector(1..dim)
      := OctoDobl_Random_Vectors.Random_Vector(1,dim);
    x : constant OctoDobl_Complex_Vectors.Link_to_Vector
      := new OctoDobl_Complex_Vectors.Vector'(v);
    ans : character;
    static,output : boolean;

  begin
    new_line;
    put_line("Generating a random complex circuit system ...");
    s := OctoDobl_Circuit_Makers.Random_Complex_System(dim,nbr,dim,pwr);
    new_line;
    put("Extra output wanted ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
    static := (ans = 'y');
    OctoDobl_Test(s,x,static,output);
  end OctoDobl_Random_Test;

  procedure DecaDobl_Random_Test
              ( dim,nbr,pwr : in integer32 ) is

    s : DecaDobl_Complex_Circuits.Link_to_System;
    v : constant DecaDobl_Complex_Vectors.Vector(1..dim)
      := DecaDobl_Random_Vectors.Random_Vector(1,dim);
    x : constant DecaDobl_Complex_Vectors.Link_to_Vector
      := new DecaDobl_Complex_Vectors.Vector'(v);
    ans : character;
    static,output : boolean;

  begin
    new_line;
    put_line("Generating a random complex circuit system ...");
    s := DecaDobl_Circuit_Makers.Random_Complex_System(dim,nbr,dim,pwr);
    new_line;
    put("Extra output wanted ? (y/n) "); Ask_Yes_or_No(ans);
    output := (ans = 'y');
    put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
    static := (ans = 'y');
    DecaDobl_Test(s,x,static,output);
  end DecaDobl_Random_Test;

  procedure Standard_User_Test is

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

  procedure TripDobl_User_Test is

    use TripDobl_Complex_Circuits;

    lp : TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : TripDobl_Complex_Solutions.Solution_List;
    ls : TripDobl_Complex_Solutions.Link_to_Solution;
    ans : character;
    static,output : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    TripDobl_System_and_Solutions_io.get(lp,sols);
    ls := TripDobl_Complex_Solutions.Head_Of(sols);
    declare
      s : constant Link_to_System
        := TripDobl_Circuit_Makers.Make_Complex_System(lp);
      x : constant TripDobl_Complex_Vectors.Link_to_Vector
        := new TripDobl_Complex_Vectors.Vector'(ls.v);
    begin
      new_line;
      put("Extra output wanted ? (y/n) "); Ask_Yes_or_No(ans);
      output := (ans = 'y');
      put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
      static := (ans = 'y');
      TripDobl_Test(s,x,static,output);
    end;
  end TripDobl_User_Test;

  procedure QuadDobl_User_Test is

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

  procedure PentDobl_User_Test is

    use PentDobl_Complex_Circuits;

    lp : PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : PentDobl_Complex_Solutions.Solution_List;
    ls : PentDobl_Complex_Solutions.Link_to_Solution;
    ans : character;
    static,output : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    PentDobl_System_and_Solutions_io.get(lp,sols);
    ls := PentDobl_Complex_Solutions.Head_Of(sols);
    declare
      s : constant Link_to_System
        := PentDobl_Circuit_Makers.Make_Complex_System(lp);
      x : constant PentDobl_Complex_Vectors.Link_to_Vector
        := new PentDobl_Complex_Vectors.Vector'(ls.v);
    begin
      new_line;
      put("Extra output wanted ? (y/n) "); Ask_Yes_or_No(ans);
      output := (ans = 'y');
      put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
      static := (ans = 'y');
      PentDobl_Test(s,x,static,output);
    end;
  end PentDobl_User_Test;

  procedure OctoDobl_User_Test is

    use OctoDobl_Complex_Circuits;

    lp : OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : OctoDobl_Complex_Solutions.Solution_List;
    ls : OctoDobl_Complex_Solutions.Link_to_Solution;
    ans : character;
    static,output : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    OctoDobl_System_and_Solutions_io.get(lp,sols);
    ls := OctoDobl_Complex_Solutions.Head_Of(sols);
    declare
      s : constant Link_to_System
        := OctoDobl_Circuit_Makers.Make_Complex_System(lp);
      x : constant OctoDobl_Complex_Vectors.Link_to_Vector
        := new OctoDobl_Complex_Vectors.Vector'(ls.v);
    begin
      new_line;
      put("Extra output wanted ? (y/n) "); Ask_Yes_or_No(ans);
      output := (ans = 'y');
      put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
      static := (ans = 'y');
      OctoDobl_Test(s,x,static,output);
    end;
  end OctoDobl_User_Test;

  procedure DecaDobl_User_Test is

    use DecaDobl_Complex_Circuits;

    lp : DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DecaDobl_Complex_Solutions.Solution_List;
    ls : DecaDobl_Complex_Solutions.Link_to_Solution;
    ans : character;
    static,output : boolean;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    DecaDobl_System_and_Solutions_io.get(lp,sols);
    ls := DecaDobl_Complex_Solutions.Head_Of(sols);
    declare
      s : constant Link_to_System
        := DecaDobl_Circuit_Makers.Make_Complex_System(lp);
      x : constant DecaDobl_Complex_Vectors.Link_to_Vector
        := new DecaDobl_Complex_Vectors.Vector'(ls.v);
    begin
      new_line;
      put("Extra output wanted ? (y/n) "); Ask_Yes_or_No(ans);
      output := (ans = 'y');
      put("Static load balancing ? (y/n) "); Ask_Yes_or_No(ans);
      static := (ans = 'y');
      DecaDobl_Test(s,x,static,output);
    end;
  end DecaDobl_User_Test;

  procedure Standard_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in Standard_Complex_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false ) is

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

  procedure TripDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in TripDobl_Complex_Circuits.Link_to_System;
                x : in TripDobl_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false ) is

    dim : constant integer32 := s.dim;
    U,V : TripDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : TripDobl_Complex_Vectors.Vector(1..dim);
    svl : constant TripDobl_Complex_VecVecs.VecVec(0..s.neq)
        := TripDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant TripDobl_Complex_VecMats.VecMat(1..s.neq)
       := TripDobl_Complex_Circuits.Allocate(s.neq,s.dim);
    values : TripDobl_Complex_VecVecs.VecVec(0..s.neq);
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    if verbose
     then put_line("Computing first without multitasking ...");
    end if;
    put_line(file,"triple double precision"); flush(file);
    seristart := Ada.Calendar.Clock;
    TripDobl_Complex_Circuits.Singular_Values(s,x,vh,U,V,e,svl);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    if verbose then
      put_line("-> Elapsed time without multitasking :");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      put_line("running in triple double precision ...");
    end if;
    put(file,"  1 : "); duration_io.put(file,seri_elapsed,1,3);
    put_line(file," : 1.000"); flush(file);
    if nbruns /= 0 then
      for k in 1..nbruns loop
        values := TripDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
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
        TripDobl_Complex_VecVecs.Clear(values);
      end loop;
    else
      for k in nbtseq'range loop
        nbt := nbtseq(k);
        values := TripDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
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
        TripDobl_Complex_VecVecs.Clear(values);
      end loop;
    end if;
  end TripDobl_Benchmark;

  procedure QuadDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in QuadDobl_Complex_Circuits.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false ) is

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

  procedure PentDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in PentDobl_Complex_Circuits.Link_to_System;
                x : in PentDobl_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false ) is

    dim : constant integer32 := s.dim;
    U,V : PentDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : PentDobl_Complex_Vectors.Vector(1..dim);
    svl : constant PentDobl_Complex_VecVecs.VecVec(0..s.neq)
        := PentDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant PentDobl_Complex_VecMats.VecMat(1..s.neq)
       := PentDobl_Complex_Circuits.Allocate(s.neq,s.dim);
    values : PentDobl_Complex_VecVecs.VecVec(0..s.neq);
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    if verbose
     then put_line("Computing first without multitasking ...");
    end if;
    put_line(file,"penta double precision"); flush(file);
    seristart := Ada.Calendar.Clock;
    PentDobl_Complex_Circuits.Singular_Values(s,x,vh,U,V,e,svl);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    if verbose then
      put_line("-> Elapsed time without multitasking :");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      put_line("running in penta double precision ...");
    end if;
    put(file,"  1 : "); duration_io.put(file,seri_elapsed,1,3);
    put_line(file," : 1.000"); flush(file);
    if nbruns /= 0 then
      for k in 1..nbruns loop
        values := PentDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
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
        PentDobl_Complex_VecVecs.Clear(values);
      end loop;
    else
      for k in nbtseq'range loop
        nbt := nbtseq(k);
        values := PentDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
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
        PentDobl_Complex_VecVecs.Clear(values);
      end loop;
    end if;
  end PentDobl_Benchmark;

  procedure OctoDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in OctoDobl_Complex_Circuits.Link_to_System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false ) is

    dim : constant integer32 := s.dim;
    U,V : OctoDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : OctoDobl_Complex_Vectors.Vector(1..dim);
    svl : constant OctoDobl_Complex_VecVecs.VecVec(0..s.neq)
        := OctoDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant OctoDobl_Complex_VecMats.VecMat(1..s.neq)
       := OctoDobl_Complex_Circuits.Allocate(s.neq,s.dim);
    values : OctoDobl_Complex_VecVecs.VecVec(0..s.neq);
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    if verbose
     then put_line("Computing first without multitasking ...");
    end if;
    put_line(file,"octo double precision"); flush(file);
    seristart := Ada.Calendar.Clock;
    OctoDobl_Complex_Circuits.Singular_Values(s,x,vh,U,V,e,svl);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    if verbose then
      put_line("-> Elapsed time without multitasking :");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      put_line("running in octo double precision ...");
    end if;
    put(file,"  1 : "); duration_io.put(file,seri_elapsed,1,3);
    put_line(file," : 1.000"); flush(file);
    if nbruns /= 0 then
      for k in 1..nbruns loop
        values := OctoDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
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
        OctoDobl_Complex_VecVecs.Clear(values);
      end loop;
    else
      for k in nbtseq'range loop
        nbt := nbtseq(k);
        values := OctoDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
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
        OctoDobl_Complex_VecVecs.Clear(values);
      end loop;
    end if;
  end OctoDobl_Benchmark;

  procedure DecaDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                nbtseq : in Standard_Integer_Vectors.Link_to_Vector;
                s : in DecaDobl_Complex_Circuits.Link_to_System;
                x : in DecaDobl_Complex_Vectors.Link_to_Vector;
                verbose : in boolean := false ) is

    dim : constant integer32 := s.dim;
    U,V : DecaDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : DecaDobl_Complex_Vectors.Vector(1..dim);
    svl : constant DecaDobl_Complex_VecVecs.VecVec(0..s.neq)
        := DecaDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
    vh : constant DecaDobl_Complex_VecMats.VecMat(1..s.neq)
       := DecaDobl_Complex_Circuits.Allocate(s.neq,s.dim);
    values : DecaDobl_Complex_VecVecs.VecVec(0..s.neq);
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup,efficiency : Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    if verbose
     then put_line("Computing first without multitasking ...");
    end if;
    put_line(file,"deca double precision"); flush(file);
    seristart := Ada.Calendar.Clock;
    DecaDobl_Complex_Circuits.Singular_Values(s,x,vh,U,V,e,svl);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    if verbose then
      put_line("-> Elapsed time without multitasking :");
      Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
      put_line("running in deca double precision ...");
    end if;
    put(file,"  1 : "); duration_io.put(file,seri_elapsed,1,3);
    put_line(file," : 1.000"); flush(file);
    if nbruns /= 0 then
      for k in 1..nbruns loop
        values := DecaDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
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
        DecaDobl_Complex_VecVecs.Clear(values);
      end loop;
    else
      for k in nbtseq'range loop
        nbt := nbtseq(k);
        values := DecaDobl_Vector_Splitters.Allocate(s.neq,dim+1,0,1);
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
        DecaDobl_Complex_VecVecs.Clear(values);
      end loop;
    end if;
  end DecaDobl_Benchmark;

  function Prompt_for_Sequence
             ( max : in integer32 )
             return Standard_Integer_Vectors.Link_to_Vector is

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

    das : DecaDobl_Complex_Circuits.Link_to_System;
    ods : OctoDobl_Complex_Circuits.Link_to_System;
    pds : PentDobl_Complex_Circuits.Link_to_System;
    qds : QuadDobl_Complex_Circuits.Link_to_System;
    tds : TripDobl_Complex_Circuits.Link_to_System;
    dds : DoblDobl_Complex_Circuits.Link_to_System;
    d_s : Standard_Complex_Circuits.Link_to_System;
   -- dcs : Standard_Coefficient_Circuits.Link_to_System;
    dax : constant DecaDobl_Complex_Vectors.Vector(1..dim)
        := DecaDobl_Random_Vectors.Random_Vector(1,dim);
    vdax : constant DecaDobl_Complex_Vectors.Link_to_Vector
         := new DecaDobl_Complex_Vectors.Vector'(dax);
    odx : constant OctoDobl_Complex_Vectors.Vector(1..dim)
        := DecaDobl_Complex_Vectors_cv.DecaDobl_Complex_to_OctoDobl(dax);
    vodx : constant OctoDobl_Complex_Vectors.Link_to_Vector
         := new OctoDobl_Complex_Vectors.Vector'(odx);
    pdx : constant PentDobl_Complex_Vectors.Vector(1..dim)
        := DecaDobl_Complex_Vectors_cv.DecaDobl_Complex_to_PentDobl(dax);
    vpdx : constant PentDobl_Complex_Vectors.Link_to_Vector
         := new PentDobl_Complex_Vectors.Vector'(pdx);
    qdx : constant QuadDobl_Complex_Vectors.Vector(1..dim)
        := DecaDobl_Complex_Vectors_cv.DecaDobl_Complex_to_QuadDobl(dax);
    vqdx : constant QuadDobl_Complex_Vectors.Link_to_Vector
         := new QuadDobl_Complex_Vectors.Vector'(qdx);
    tdx : constant TripDobl_Complex_Vectors.Vector(1..dim)
        := DecaDobl_Complex_Vectors_cv.DecaDobl_Complex_to_TripDobl(dax);
    vtdx : constant TripDobl_Complex_Vectors.Link_to_Vector
         := new TripDobl_Complex_Vectors.Vector'(tdx);
    ddx : constant DoblDobl_Complex_Vectors.Vector(1..dim)
        := DecaDobl_Complex_Vectors_cv.DecaDobl_Complex_to_DoblDobl(dax);
    vddx : constant DoblDobl_Complex_Vectors.Link_to_Vector
         := new DoblDobl_Complex_Vectors.Vector'(ddx);
    dx : constant Standard_Complex_Vectors.Vector(1..dim)
        := DecaDobl_Complex_Vectors_cv.DecaDobl_Complex_to_Standard(dax);
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
    das := DecaDobl_Circuit_Makers.Random_Complex_System(dim,nbr,dim,pwr);
    ods := DecaDobl_Circuit_Makers.to_octo_double(das);
    pds := DecaDobl_Circuit_Makers.to_penta_double(das);
    qds := DecaDobl_Circuit_Makers.to_quad_double(das);
    tds := DecaDobl_Circuit_Makers.to_triple_double(das);
    dds := DecaDobl_Circuit_Makers.to_double_double(das);
    d_s := DecaDobl_Circuit_Makers.to_double(das);
   -- dcs := Standard_Circuit_Splitters.Split(d_s);
    Standard_Benchmark(file,nbruns,inc,nbtseq,d_s,vdx);
   -- Standard_Coefficient_Benchmark(file,nbruns,inc,nbtseq,dcs,vdx);
    DoblDobl_Benchmark(file,nbruns,inc,nbtseq,dds,vddx);
    TripDobl_Benchmark(file,nbruns,inc,nbtseq,tds,vtdx);
    QuadDobl_Benchmark(file,nbruns,inc,nbtseq,qds,vqdx);
    PentDobl_Benchmark(file,nbruns,inc,nbtseq,pds,vpdx);
    OctoDobl_Benchmark(file,nbruns,inc,nbtseq,ods,vodx);
    DecaDobl_Benchmark(file,nbruns,inc,nbtseq,das,vdax);
  end Benchmark;

  procedure Benchmark
              ( p : in DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in DecaDobl_Complex_Solutions.Solution_List ) is

    file : file_type;
    das : DecaDobl_Complex_Circuits.Link_to_System;
    ods : OctoDobl_Complex_Circuits.Link_to_System;
    pds : PentDobl_Complex_Circuits.Link_to_System;
    qds : QuadDobl_Complex_Circuits.Link_to_System;
    tds : TripDobl_Complex_Circuits.Link_to_System;
    dds : DoblDobl_Complex_Circuits.Link_to_System;
    d_s : Standard_Complex_Circuits.Link_to_System;
   -- dcs : Standard_Coefficient_Circuits.Link_to_System;
    ls : constant DecaDobl_Complex_Solutions.Link_to_Solution
       := DecaDobl_Complex_Solutions.Head_Of(sols);
    dax : constant DecaDobl_Complex_Vectors.Link_to_Vector
        := new DecaDobl_Complex_Vectors.Vector'(ls.v);
    dim : constant integer32 := ls.n;
    odx : constant OctoDobl_Complex_Vectors.Vector(1..dim)
        := DecaDobl_Complex_Vectors_cv.DecaDobl_Complex_to_OctoDobl(ls.v);
    vodx : constant OctoDobl_Complex_Vectors.Link_to_Vector
         := new OctoDobl_Complex_Vectors.Vector'(odx);
    pdx : constant PentDobl_Complex_Vectors.Vector(1..dim)
        := DecaDobl_Complex_Vectors_cv.DecaDobl_Complex_to_PentDobl(ls.v);
    vpdx : constant PentDobl_Complex_Vectors.Link_to_Vector
         := new PentDobl_Complex_Vectors.Vector'(pdx);
    qdx : constant QuadDobl_Complex_Vectors.Vector(1..dim)
        := DecaDobl_Complex_Vectors_cv.DecaDobl_Complex_to_QuadDobl(ls.v);
    vqdx : constant QuadDobl_Complex_Vectors.Link_to_Vector
         := new QuadDobl_Complex_Vectors.Vector'(qdx);
    tdx : constant TripDobl_Complex_Vectors.Vector(1..dim)
        := DecaDobl_Complex_Vectors_cv.DecaDobl_Complex_to_TripDobl(ls.v);
    vtdx : constant TripDobl_Complex_Vectors.Link_to_Vector
         := new TripDobl_Complex_Vectors.Vector'(tdx);
    ddx : constant DoblDobl_Complex_Vectors.Vector(1..dim)
        := DecaDobl_Complex_Vectors_cv.DecaDobl_Complex_to_DoblDobl(ls.v);
    vddx : constant DoblDobl_Complex_Vectors.Link_to_Vector
         := new DoblDobl_Complex_Vectors.Vector'(ddx);
    dx : constant Standard_Complex_Vectors.Vector(1..dim)
        := DecaDobl_Complex_Vectors_cv.DecaDobl_Complex_to_Standard(ls.v);
    vdx : constant Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector'(dx);
    nbruns,inc : integer32 := 0;
    nbtseq : Standard_Integer_Vectors.Link_to_Vector;
 
  begin
    new_line;
    put_line("making circuits for the given system ...");
    das := DecaDobl_Circuit_Makers.Make_Complex_System(p);
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
    ods := DecaDobl_Circuit_Makers.to_octo_double(das);
    pds := DecaDobl_Circuit_Makers.to_penta_double(das);
    qds := DecaDobl_Circuit_Makers.to_quad_double(das);
    tds := DecaDobl_Circuit_Makers.to_triple_double(das);
    dds := DecaDobl_Circuit_Makers.to_double_double(das);
    d_s := Standard_Circuit_Makers.to_double(qds);
   -- dcs := Standard_Circuit_Splitters.Split(d_s);
    Standard_Benchmark(file,nbruns,inc,nbtseq,d_s,vdx);
   -- Standard_Coefficient_Benchmark(file,nbruns,inc,nbtseq,dcs,vdx);
    DoblDobl_Benchmark(file,nbruns,inc,nbtseq,dds,vddx);
    TripDobl_Benchmark(file,nbruns,inc,nbtseq,tds,vtdx);
    QuadDobl_Benchmark(file,nbruns,inc,nbtseq,qds,vqdx);
    PentDobl_Benchmark(file,nbruns,inc,nbtseq,pds,vpdx);
    OctoDobl_Benchmark(file,nbruns,inc,nbtseq,ods,vodx);
    DecaDobl_Benchmark(file,nbruns,inc,nbtseq,das,dax);
  end Benchmark;

  function Prompt_for_Precision return character is

    precision : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  1. double precision");
    put_line("  2. double double precision");
    put_line("  3. triple double precision");
    put_line("  4. quad double precision");
    put_line("  5. penta double precision");
    put_line("  6. octo double precision");
    put_line("  7. deca double precision");
    put("Type 1, 2, 3, 4, 5, 6, or 7 to select the precision : ");
    Ask_Alternative(precision,"1234567");
    return precision;
  end Prompt_for_Precision;

  procedure Main is

    dim,nbr,pwr : integer32 := 0;
    ans,precision : character;
    lp : DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DecaDobl_Complex_Solutions.Solution_List;

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
        DecaDobl_System_and_Solutions_io.get(lp,sols);
        Benchmark(lp,sols);
      else
        precision := Prompt_for_Precision;
        case precision is
          when '1' => Standard_User_Test;
          when '2' => DoblDobl_User_Test;
          when '3' => TripDobl_User_Test;
          when '4' => QuadDobl_User_Test;
          when '5' => PentDobl_User_Test;
          when '6' => OctoDobl_User_Test;
          when '7' => DecaDobl_User_Test;
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
          when '1' => Standard_Random_Test(dim,nbr,pwr);
          when '2' => DoblDobl_Random_Test(dim,nbr,pwr);
          when '3' => TripDobl_Random_Test(dim,nbr,pwr);
          when '4' => QuadDobl_Random_Test(dim,nbr,pwr);
          when '5' => PentDobl_Random_Test(dim,nbr,pwr);
          when '6' => OctoDobl_Random_Test(dim,nbr,pwr);
          when '7' => DecaDobl_Random_Test(dim,nbr,pwr);
          when others => null;
        end case;
      end if;
    end if;
  end Main;

end Test_mtHessian_Circuits;
