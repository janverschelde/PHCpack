with duration_io;
with Ada.Calendar;
with Time_Stamps;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
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
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors_cv;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_System_and_Solutions_io;
with System_Convolution_Circuits;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Evaluation_Differentiation_Errors;  use Evaluation_Differentiation_Errors;
with Hessian_Convolution_Circuits;       use Hessian_Convolution_Circuits;
with Jacobian_Convolution_Circuits;      use Jacobian_Convolution_Circuits;
with Multitasked_Hessian_Convolutions;   use Multitasked_Hessian_Convolutions;
with Standard_Newton_Convolutions;
with DoblDobl_Newton_Convolutions;
with QuadDobl_Newton_Convolutions;

package body Test_mtHessian_Convolutions is

  procedure Write_Singular_Values 
              ( values : in Standard_Complex_VecVecs.VecVec;
                jmvals : in Standard_Complex_Vectors.Vector ) is

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

  procedure Standard_Test
              ( dim : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_VecVecs.Link_to_VecVec ) is

    use Standard_Complex_VecVecs;
    use Standard_Speelpenning_Convolutions;

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
  end Standard_Test;

  procedure DoblDobl_Test
              ( dim : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_VecVecs.Link_to_VecVec ) is

    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Speelpenning_Convolutions;

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
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( dim : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_VecVecs.Link_to_VecVec ) is

    use QuadDobl_Complex_VecVecs;
    use QuadDobl_Speelpenning_Convolutions;

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
  end QuadDobl_Test;

  procedure Standard_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

    s : Standard_Speelpenning_Convolutions.Link_to_System;
    x : Standard_Complex_VecVecs.Link_to_VecVec;

  begin
    new_line;
    put_line("Generating a random Newton homotopy ...");
    Standard_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    Standard_Test(dim,s,x);
  end Standard_Random_Test;

  procedure DoblDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

    s : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    x : DoblDobl_Complex_VecVecs.Link_to_VecVec;

  begin
    new_line;
    put_line("Generating a random Newton homotopy ...");
    DoblDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    DoblDobl_Test(dim,s,x);
  end DoblDobl_Random_Test;

  procedure QuadDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

    s : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    x : QuadDobl_Complex_VecVecs.Link_to_VecVec;

  begin
    new_line;
    put_line("Generating a random Newton homotopy ...");
    QuadDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    QuadDobl_Test(dim,s,x);
  end QuadDobl_Random_Test;

  procedure Standard_User_Test is

    use System_Convolution_Circuits;
    use Standard_Speelpenning_Convolutions;

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    dim,deg : integer32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    Standard_System_and_Solutions_io.get(lp,sols);
    ls := Standard_Complex_Solutions.Head_Of(sols);
    dim := ls.n;
    new_line;
    put("Give the degree of the series : "); get(deg);
    declare
      c : constant Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(deg));
      s : constant Link_to_System := Create(c,dim,deg);
      scf : constant Standard_Complex_VecVecs.VecVec(1..dim)
          := Standard_Newton_Convolutions.Series_Coefficients(ls.v,deg);
      x : constant Standard_Complex_VecVecs.Link_to_VecVec
        := new Standard_Complex_VecVecs.VecVec'(scf);
    begin
      Standard_Test(dim,s,x);
    end;
  end Standard_User_Test;

  procedure DoblDobl_User_Test is

    use System_Convolution_Circuits;
    use DoblDobl_Speelpenning_Convolutions;

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    dim,deg : integer32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    DoblDobl_System_and_Solutions_io.get(lp,sols);
    ls := DoblDobl_Complex_Solutions.Head_Of(sols);
    dim := ls.n;
    new_line;
    put("Give the degree of the series : "); get(deg);
    declare
      c : constant Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(deg));
      s : constant Link_to_System := Create(c,dim,deg);
      scf : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
          := DoblDobl_Newton_Convolutions.Series_Coefficients(ls.v,deg);
      x : constant DoblDobl_Complex_VecVecs.Link_to_VecVec
        := new DoblDobl_Complex_VecVecs.VecVec'(scf);
    begin
      DoblDobl_Test(dim,s,x);
    end;
  end DoblDobl_User_Test;

  procedure QuadDobl_User_Test is

    use System_Convolution_Circuits;
    use QuadDobl_Speelpenning_Convolutions;

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    dim,deg : integer32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    QuadDobl_System_and_Solutions_io.get(lp,sols);
    ls := QuadDobl_Complex_Solutions.Head_Of(sols);
    dim := ls.n;
    new_line;
    put("Give the degree of the series : "); get(deg);
    declare
      c : constant Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(deg));
      s : constant Link_to_System := Create(c,dim,deg);
      scf : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
          := QuadDobl_Newton_Convolutions.Series_Coefficients(ls.v,deg);
      x : constant QuadDobl_Complex_VecVecs.Link_to_VecVec
        := new QuadDobl_Complex_VecVecs.VecVec'(scf);
    begin
      QuadDobl_Test(dim,s,x);
    end;
  end QuadDobl_User_Test;

  procedure Standard_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_Vectors.Vector;
                verbose : in boolean := false ) is

    dim : constant integer32 := s.dim;
    A,U,V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    e : Standard_Complex_Vectors.Vector(1..dim);
    svl : constant Standard_Complex_VecVecs.VecVec(1..dim) := Allocate(dim,dim);
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    jm1vls : Standard_Complex_Vectors.Vector(1..dim);
    jm2vls : Standard_Complex_Vectors.Vector(1..dim);
    values : Standard_Complex_VecVecs.VecVec(1..dim);
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup :  Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    if verbose
     then put_line(file,"Computing first without multitasking ...");
     else put_line(file,"double precision");
    end if;
    seristart := Ada.Calendar.Clock;
    Singular_Values(s.crc,x,A,U,V,e,jm1vls);
    for i in 1..dim loop
      lnk := svl(i);
      Singular_Values(s.crc(i),x,A,U,V,e,lnk.all);
    end loop;
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    if verbose then
      put_line(file,"-> Elapsed time without multitasking :");
      Time_Stamps.Write_Elapsed_Time(file,seristart,seristop);
      put_line(file,"running in double precision ...");
    else
      put(file,"  1 : ");
      duration_io.put(file,seri_elapsed,1,3); put_line(file," : 1.000");
    end if;
    for k in 1..nbruns loop
      values := Allocate(dim,dim);
      multstart := Ada.Calendar.Clock;
      Multitasked_Singular_Values(nbt,s,x,jm2vls,values,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        put(file,"-> Elapsed time with "); put(file,nbt,1);
        put_line(file," tasks :");
        Time_Stamps.Write_Elapsed_Time(file,multstart,multstop);
      else
        put(file,nbt,3); put(file," : ");
        duration_io.put(file,mult_elapsed,1,3); put(file," : ");
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        if verbose
         then put(file,"The speedup : ");
        end if;
        duration_io.put(file,speedup,1,3); new_line(file);
      end if;
      nbt := nbt + inc;
      Standard_Complex_VecVecs.Clear(values);
    end loop;
  end Standard_Benchmark;

  procedure DoblDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Vector;
                verbose : in boolean := false ) is

    dim : constant integer32 := s.dim;
    A,U,V : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : DoblDobl_Complex_Vectors.Vector(1..dim);
    svl : constant DoblDobl_Complex_VecVecs.VecVec(1..dim) := Allocate(dim,dim);
    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    jm1vls : DoblDobl_Complex_Vectors.Vector(1..dim);
    jm2vls : DoblDobl_Complex_Vectors.Vector(1..dim);
    values : DoblDobl_Complex_VecVecs.VecVec(1..dim);
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup :  Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    if verbose
     then put_line(file,"Computing first without multitasking ...");
     else put_line(file,"double double precision");
    end if;
    seristart := Ada.Calendar.Clock;
    Singular_Values(s.crc,x,A,U,V,e,jm1vls);
    for i in 1..dim loop
      lnk := svl(i);
      Singular_Values(s.crc(i),x,A,U,V,e,lnk.all);
    end loop;
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    if verbose then
      put_line(file,"-> Elapsed time without multitasking :");
      Time_Stamps.Write_Elapsed_Time(file,seristart,seristop);
      put_line(file,"running in double double precision ...");
    else
      put(file,"  1 : ");
      duration_io.put(file,seri_elapsed,1,3); put_line(file," : 1.000");
    end if;
    for k in 1..nbruns loop
      values := Allocate(dim,dim);
      multstart := Ada.Calendar.Clock;
      Multitasked_Singular_Values(nbt,s,x,jm2vls,values,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        put(file,"-> Elapsed time with "); put(file,nbt,1);
        put_line(file," tasks :");
        Time_Stamps.Write_Elapsed_Time(file,multstart,multstop);
      else
        put(file,nbt,3); put(file," : ");
        duration_io.put(file,mult_elapsed,1,3); put(file," : ");
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        if verbose
         then put(file,"The speedup : ");
        end if;
        duration_io.put(file,speedup,1,3); new_line(file);
      end if;
      nbt := nbt + inc;
      DoblDobl_Complex_VecVecs.Clear(values);
    end loop;
  end DoblDobl_Benchmark;

  procedure QuadDobl_Benchmark
              ( file : in file_type; nbruns,inc : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Vector;
                verbose : in boolean := false ) is

    dim : constant integer32 := s.dim;
    A,U,V : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : QuadDobl_Complex_Vectors.Vector(1..dim);
    svl : constant QuadDobl_Complex_VecVecs.VecVec(1..dim) := Allocate(dim,dim);
    lnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    jm1vls : QuadDobl_Complex_Vectors.Vector(1..dim);
    jm2vls : QuadDobl_Complex_Vectors.Vector(1..dim);
    values : QuadDobl_Complex_VecVecs.VecVec(1..dim);
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup :  Duration;
    nbt : integer32 := 2;

    use Ada.Calendar;

  begin
    if verbose
     then put_line(file,"Computing first without multitasking ...");
     else put_line(file,"quad double precision");
    end if;
    seristart := Ada.Calendar.Clock;
    Singular_Values(s.crc,x,A,U,V,e,jm1vls);
    for i in 1..dim loop
      lnk := svl(i);
      Singular_Values(s.crc(i),x,A,U,V,e,lnk.all);
    end loop;
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    if verbose then
      put_line(file,"-> Elapsed time without multitasking :");
      Time_Stamps.Write_Elapsed_Time(file,seristart,seristop);
      put_line(file,"running in quad double precision ...");
    else
      put(file,"  1 : ");
      duration_io.put(file,seri_elapsed,1,3); put_line(file," : 1.000");
    end if;
    for k in 1..nbruns loop
      values := Allocate(dim,dim);
      multstart := Ada.Calendar.Clock;
      Multitasked_Singular_Values(nbt,s,x,jm2vls,values,false);
      multstop := Ada.Calendar.Clock;
      mult_elapsed := multstop - multstart;
      if verbose then
        put(file,"-> Elapsed time with "); put(file,nbt,1);
        put_line(file," tasks :");
        Time_Stamps.Write_Elapsed_Time(file,multstart,multstop);
      else
        put(file,nbt,3); put(file," : ");
        duration_io.put(file,mult_elapsed,1,3); put(file," : ");
      end if;
      if seri_elapsed + 1.0 /= 1.0 then
        speedup := seri_elapsed/mult_elapsed;
        if verbose
         then put(file,"The speedup : ");
        end if;
        duration_io.put(file,speedup,1,3); new_line(file);
      end if;
      nbt := nbt + inc;
      QuadDobl_Complex_VecVecs.Clear(values);
    end loop;
  end QuadDobl_Benchmark;

  procedure Benchmark ( dim,deg,nbr,pwr : in integer32 ) is

    qds : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    dds : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    d_s : Standard_Speelpenning_Convolutions.Link_to_System;
    qdx : QuadDobl_Complex_VecVecs.Link_to_VecVec;
    ddx : DoblDobl_Complex_VecVecs.Link_to_VecVec;
    d_x : Standard_Complex_VecVecs.Link_to_VecVec;
    qdlnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    ddlnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    dlnk : Standard_Complex_Vectors.Link_to_Vector;
    qdvx : QuadDobl_Complex_Vectors.Vector(1..dim);
    ddvx : DoblDobl_Complex_Vectors.Vector(1..dim);
    dvx : Standard_Complex_Vectors.Vector(1..dim);
    file : file_type;
    nbruns,inc : integer32 := 0;
 
  begin
    new_line;
    put("Give the number of multitasked runs : "); get(nbruns);
    put("Give the increment on the tasks : "); get(inc); skip_line;
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    QuadDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,qds,qdx);
    dds := System_Convolution_Circuits.to_double_double(qds);
    ddx := QuadDobl_Complex_Vectors_cv.to_double_double(qdx);
    d_s := System_Convolution_Circuits.to_double(qds);
    d_x := QuadDobl_Complex_Vectors_cv.to_double(qdx);
    for i in 1..dim loop
      qdlnk := qdx(i); ddlnk := ddx(i); dlnk := d_x(i);
      qdvx(i) := qdlnk(0); ddvx(i) := ddlnk(0); dvx(i) := dlnk(0);
    end loop;
    Standard_Benchmark(file,nbruns,inc,d_s,dvx);
    DoblDobl_Benchmark(file,nbruns,inc,dds,ddvx);
    QuadDobl_Benchmark(file,nbruns,inc,qds,qdvx);
  end Benchmark;

  function Prompt_for_Precision return character is

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

    dim,deg,nbr,pwr : integer32 := 0;
    ans,precision : character;

  begin
    new_line;
    put_line("Testing the Hessian criterion ...");
    new_line;
    put("Generate a random system ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'n' then
      precision := Prompt_for_Precision;
      case precision is
        when '0' => Standard_User_Test;
        when '1' => DoblDobl_User_Test;
        when '2' => QuadDobl_User_Test;
        when others => null;
      end case;
    else
      put("Give the dimension : "); get(dim);
      put("Give the degree of the power series : "); get(deg);
      put("Give the number of terms in each circuit : "); get(nbr);
      put("Give the largest power of the variables : "); get(pwr);
      new_line;
      put("Benchmark for all precisions ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y' then
        Benchmark(dim,deg,nbr,pwr);
      else
        precision := Prompt_for_Precision;
        case precision is
          when '0' => Standard_Random_Test(dim,deg,nbr,pwr);
          when '1' => DoblDobl_Random_Test(dim,deg,nbr,pwr);
          when '2' => QuadDobl_Random_Test(dim,deg,nbr,pwr);
          when others => null;
        end case;
      end if;
    end if;
  end Main;

end Test_mtHessian_Convolutions;
