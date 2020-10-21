with duration_io;
with Ada.Calendar;
with Communications_with_User;           use Communications_with_User;
with Time_Stamps;
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
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Random_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with DoblDobl_Random_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with QuadDobl_Random_Vectors;
with Evaluation_Differentiation_Errors;  use Evaluation_Differentiation_Errors;
with Standard_Rational_Approximations;
with DoblDobl_Rational_Approximations;
with QuadDobl_Rational_Approximations;
with Multitasked_Pade_Approximations;    use Multitasked_Pade_Approximations;

package body Test_mtPade_Approximations is

  procedure Standard_Allocate
              ( nbr,numdeg,dendeg : in integer32;
                cff : out Standard_Complex_VecVecs.VecVec;
                numcff1,numcff2 : out Standard_Complex_VecVecs.VecVec;
                dencff1,dencff2 : out Standard_Complex_VecVecs.VecVec ) is

    dim : constant integer32 := numdeg + dendeg;

  begin
    for k in 1..nbr loop
      declare
        kcff : constant Standard_Complex_Vectors.Vector(0..dim)
             := Standard_Random_Vectors.Random_Vector(0,dim);
        knum1 : constant Standard_Complex_Vectors.Vector(0..numdeg)
              := (0..numdeg => Standard_Complex_Numbers.Create(integer(0)));
        knum2 : constant Standard_Complex_Vectors.Vector(0..numdeg)
              := (0..numdeg => Standard_Complex_Numbers.Create(integer(0)));
        kden1 : constant Standard_Complex_Vectors.Vector(0..dendeg)
              := (0..numdeg => Standard_Complex_Numbers.Create(integer(0)));
        kden2 : constant Standard_Complex_Vectors.Vector(0..dendeg)
              := (0..numdeg => Standard_Complex_Numbers.Create(integer(0)));
      begin
        cff(k) := new Standard_Complex_Vectors.Vector'(kcff);
        numcff1(k) := new Standard_Complex_Vectors.Vector'(knum1);
        numcff2(k) := new Standard_Complex_Vectors.Vector'(knum2);
        dencff1(k) := new Standard_Complex_Vectors.Vector'(kden1);
        dencff2(k) := new Standard_Complex_Vectors.Vector'(kden2);
      end;
    end loop;
  end Standard_Allocate;

  procedure DoblDobl_Allocate
              ( nbr,numdeg,dendeg : in integer32;
                cff : out DoblDobl_Complex_VecVecs.VecVec;
                numcff1,numcff2 : out DoblDobl_Complex_VecVecs.VecVec;
                dencff1,dencff2 : out DoblDobl_Complex_VecVecs.VecVec ) is

    dim : constant integer32 := numdeg + dendeg;

  begin
    for k in 1..nbr loop
      declare
        kcff : constant DoblDobl_Complex_Vectors.Vector(0..dim)
             := DoblDobl_Random_Vectors.Random_Vector(0,dim);
        knum1 : constant DoblDobl_Complex_Vectors.Vector(0..numdeg)
              := (0..numdeg => DoblDobl_Complex_Numbers.Create(integer(0)));
        knum2 : constant DoblDobl_Complex_Vectors.Vector(0..numdeg)
              := (0..numdeg => DoblDobl_Complex_Numbers.Create(integer(0)));
        kden1 : constant DoblDobl_Complex_Vectors.Vector(0..dendeg)
              := (0..numdeg => DoblDobl_Complex_Numbers.Create(integer(0)));
        kden2 : constant DoblDobl_Complex_Vectors.Vector(0..dendeg)
              := (0..numdeg => DoblDobl_Complex_Numbers.Create(integer(0)));
      begin
        cff(k) := new DoblDobl_Complex_Vectors.Vector'(kcff);
        numcff1(k) := new DoblDobl_Complex_Vectors.Vector'(knum1);
        numcff2(k) := new DoblDobl_Complex_Vectors.Vector'(knum2);
        dencff1(k) := new DoblDobl_Complex_Vectors.Vector'(kden1);
        dencff2(k) := new DoblDobl_Complex_Vectors.Vector'(kden2);
      end;
    end loop;
  end DoblDobl_Allocate;

  procedure QuadDobl_Allocate
              ( nbr,numdeg,dendeg : in integer32;
                cff : out QuadDobl_Complex_VecVecs.VecVec;
                numcff1,numcff2 : out QuadDobl_Complex_VecVecs.VecVec;
                dencff1,dencff2 : out QuadDobl_Complex_VecVecs.VecVec ) is

    dim : constant integer32 := numdeg + dendeg;

  begin
    null;
    for k in 1..nbr loop
      declare
        kcff : constant QuadDobl_Complex_Vectors.Vector(0..dim)
             := QuadDobl_Random_Vectors.Random_Vector(0,dim);
        knum1 : constant QuadDobl_Complex_Vectors.Vector(0..numdeg)
              := (0..numdeg => QuadDobl_Complex_Numbers.Create(integer(0)));
        knum2 : constant QuadDobl_Complex_Vectors.Vector(0..numdeg)
              := (0..numdeg => QuadDobl_Complex_Numbers.Create(integer(0)));
        kden1 : constant QuadDobl_Complex_Vectors.Vector(0..dendeg)
              := (0..numdeg => QuadDobl_Complex_Numbers.Create(integer(0)));
        kden2 : constant QuadDobl_Complex_Vectors.Vector(0..dendeg)
              := (0..numdeg => QuadDobl_Complex_Numbers.Create(integer(0)));
      begin
        cff(k) := new QuadDobl_Complex_Vectors.Vector'(kcff);
        numcff1(k) := new QuadDobl_Complex_Vectors.Vector'(knum1);
        numcff2(k) := new QuadDobl_Complex_Vectors.Vector'(knum2);
        dencff1(k) := new QuadDobl_Complex_Vectors.Vector'(kden1);
        dencff2(k) := new QuadDobl_Complex_Vectors.Vector'(kden2);
      end;
    end loop;
  end QuadDobl_Allocate;

  procedure Standard_Test ( nbr,numdeg,dendeg : in integer32 ) is

    cff : Standard_Complex_VecVecs.VecVec(1..nbr);
    numcff1,numcff2 : Standard_Complex_VecVecs.VecVec(1..nbr);
    dencff1,dencff2 : Standard_Complex_VecVecs.VecVec(1..nbr);
    mat : Standard_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : Standard_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;
    nbtasks : integer32 := 0;
    err : double_float;
    ans : character;
    t : constant double_float := 0.1;
    eva1,eva2 : Standard_Complex_Vectors.Vector(1..nbr);
    otp : boolean;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    serelp,mltelp,speedup,efficiency : duration;

    use Ada.Calendar;
    use Standard_Complex_Numbers;
    use Standard_Rational_Approximations;

  begin
    new_line;
    put_line("Allocating and generating data ...");
    Standard_Allocate(nbr,numdeg,dendeg,cff,numcff1,numcff2,dencff1,dencff2);
    new_line;
    put_line("Computing without multitasking ...");
    seristart := Ada.Calendar.Clock;
    Pade_Vector(numdeg,dendeg,cff,numcff1,dencff1,mat,rhs,ipvt,info,false);
    Evaluate(numcff1,dencff1,t,eva1);
    seristop := Ada.Calendar.Clock;
    serelp := seristop - seristart;
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    new_line;
    put("Give the number of tasks : "); get(nbtasks);
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    multstart := Ada.Calendar.Clock;
    Standard_Multitasking(nbtasks,numdeg,dendeg,cff,numcff2,dencff2,t,eva2,otp);
    multstop := Ada.Calendar.Clock;
    mltelp := multstop - multstart;
    Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
    if serelp + 1.0 /= 1.0 then
      speedup := serelp/mltelp;
      put("The speedup : "); duration_io.put(speedup,1,3);
      efficiency := speedup/duration(nbtasks);
      efficiency := duration(100)*efficiency;
      put("  the efficiency : "); duration_io.put(efficiency,2,2); new_line;
    end if;
    err := Difference(numcff1,numcff2);
    put("  the error in numerator coefficients : "); put(err,3); new_line;
    err := Difference(dencff1,dencff2);
    put("  the error in denominator coefficients : "); put(err,3); new_line;
    err := 0.0;
    for k in 1..nbr loop
      err := err + AbsVal(eva1(k) - eva2(k));
    end loop;
    put("  the error in the evaluation : "); put(err,3); new_line;
  end Standard_Test;

  procedure DoblDobl_Test ( nbr,numdeg,dendeg : in integer32 ) is

    cff : DoblDobl_Complex_VecVecs.VecVec(1..nbr);
    numcff1,numcff2 : DoblDobl_Complex_VecVecs.VecVec(1..nbr);
    dencff1,dencff2 : DoblDobl_Complex_VecVecs.VecVec(1..nbr);
    mat : DoblDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : DoblDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;
    nbtasks : integer32 := 0;
    err : double_double;
    ans : character;
    t : constant double_double := create(0.1);
    eva1,eva2 : DoblDobl_Complex_Vectors.Vector(1..nbr);
    otp : boolean;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    serelp,mltelp,speedup,efficiency : duration;

    use Ada.Calendar;
    use DoblDobl_Complex_Numbers;
    use DoblDobl_Rational_Approximations;

  begin
    new_line;
    put_line("Allocating and generating data ...");
    DoblDobl_Allocate(nbr,numdeg,dendeg,cff,numcff1,numcff2,dencff1,dencff2);
    new_line;
    put_line("Computing without multitasking ...");
    seristart := Ada.Calendar.Clock;
    Pade_Vector(numdeg,dendeg,cff,numcff1,dencff1,mat,rhs,ipvt,info,false);
    Evaluate(numcff1,dencff1,t,eva1);
    seristop := Ada.Calendar.Clock;
    serelp := seristop - seristart;
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    new_line;
    put("Give the number of tasks : "); get(nbtasks);
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    multstart := Ada.Calendar.Clock;
    DoblDobl_Multitasking(nbtasks,numdeg,dendeg,cff,numcff2,dencff2,t,eva2,otp);
    multstop := Ada.Calendar.Clock;
    mltelp := multstop - multstart;
    Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
    if serelp + 1.0 /= 1.0 then
      speedup := serelp/mltelp;
      put("The speedup : "); duration_io.put(speedup,1,3);
      efficiency := speedup/duration(nbtasks);
      efficiency := duration(100)*efficiency;
      put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
    end if;
    err := Difference(numcff1,numcff2);
    put("  the error in numerator coefficients : "); put(err,3); new_line;
    err := Difference(dencff1,dencff2);
    put("  the error in denominator coefficients : "); put(err,3); new_line;
    err := create(0.0);
    for k in 1..nbr loop
      err := err + AbsVal(eva1(k) - eva2(k));
    end loop;
    put("  the error in the evaluation : "); put(err,3); new_line;
  end DoblDobl_Test;

  procedure QuadDobl_Test ( nbr,numdeg,dendeg : in integer32 ) is

    cff : QuadDobl_Complex_VecVecs.VecVec(1..nbr);
    numcff1,numcff2 : QuadDobl_Complex_VecVecs.VecVec(1..nbr);
    dencff1,dencff2 : QuadDobl_Complex_VecVecs.VecVec(1..nbr);
    mat : QuadDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : QuadDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;
    nbtasks : integer32 := 0;
    err : quad_double;
    ans : character;
    t : constant quad_double := create(0.1);
    eva1,eva2 : QuadDobl_Complex_Vectors.Vector(1..nbr);
    otp : boolean;
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    serelp,mltelp,speedup,efficiency : duration;

    use Ada.Calendar;
    use QuadDobl_Complex_Numbers;
    use QuadDobl_Rational_Approximations;

  begin
    new_line;
    put_line("Allocating and generating data ...");
    QuadDobl_Allocate(nbr,numdeg,dendeg,cff,numcff1,numcff2,dencff1,dencff2);
    new_line;
    put_line("Computing without multitasking ...");
    seristart := Ada.Calendar.Clock;
    Pade_Vector(numdeg,dendeg,cff,numcff1,dencff1,mat,rhs,ipvt,info,false);
    Evaluate(numcff1,dencff1,t,eva1);
    seristop := Ada.Calendar.Clock;
    serelp := seristop - seristart;
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    new_line;
    put("Give the number of tasks : "); get(nbtasks);
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    multstart := Ada.Calendar.Clock;
    QuadDobl_Multitasking(nbtasks,numdeg,dendeg,cff,numcff2,dencff2,t,eva2,otp);
    multstop := Ada.Calendar.Clock;
    mltelp := multstop - multstart;
    Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
    if serelp + 1.0 /= 1.0 then
      speedup := serelp/mltelp;
      put("The speedup : "); duration_io.put(speedup,1,3);
      efficiency := speedup/duration(nbtasks);
      efficiency := duration(100)*efficiency;
      put("  efficiency : "); duration_io.put(efficiency,2,2); new_line;
    end if;
    err := Difference(numcff1,numcff2);
    put("  the error in numerator coefficients : "); put(err,3); new_line;
    err := Difference(dencff1,dencff2);
    put("  the error in denominator coefficients : "); put(err,3); new_line;
    err := create(0.0);
    for k in 1..nbr loop
      err := err + AbsVal(eva1(k) - eva2(k));
    end loop;
    put("  the error in the evaluation : "); put(err,3); new_line;
  end QuadDobl_Test;

  procedure Standard_Benchmark
              ( file : in file_type;
                nbr,numdeg,dendeg,nbruns,inc : in integer32 ) is

    cff : Standard_Complex_VecVecs.VecVec(1..nbr);
    numcff1,numcff2 : Standard_Complex_VecVecs.VecVec(1..nbr);
    dencff1,dencff2 : Standard_Complex_VecVecs.VecVec(1..nbr);
    mat : Standard_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : Standard_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;
    t : constant double_float := 0.1;
    eva1,eva2 : Standard_Complex_Vectors.Vector(1..nbr);
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    serelp,mltelp,speedup,efficiency : duration;
    nbt : integer32 := 2;

    use Ada.Calendar;
    use Standard_Rational_Approximations;

  begin
    put_line(file,"double precision");
    Standard_Allocate(nbr,numdeg,dendeg,cff,numcff1,numcff2,dencff1,dencff2);
    seristart := Ada.Calendar.Clock;
    Pade_Vector(numdeg,dendeg,cff,numcff1,dencff1,mat,rhs,ipvt,info,false);
    Evaluate(numcff1,dencff1,t,eva1);
    seristop := Ada.Calendar.Clock;
    serelp := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,serelp,1,3); new_line(file); flush(file);
    for k in 1..nbruns loop
      multstart := Ada.Calendar.Clock;
      Standard_Multitasking(nbt,numdeg,dendeg,cff,numcff2,dencff2,t,eva2,false);
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
      if serelp + 1.0 /= 1.0 then
        speedup := serelp/mltelp;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
      else
         speedup := 0.0; efficiency := 0.0;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mltelp,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3);
      put(file," : "); duration_io.put(file,efficiency,2,2);
      new_line(file); flush(file);
      nbt := nbt + inc;
    end loop;
  end Standard_Benchmark;

  procedure DoblDobl_Benchmark
              ( file : in file_type;
                nbr,numdeg,dendeg,nbruns,inc : in integer32 ) is

    cff : DoblDobl_Complex_VecVecs.VecVec(1..nbr);
    numcff1,numcff2 : DoblDobl_Complex_VecVecs.VecVec(1..nbr);
    dencff1,dencff2 : DoblDobl_Complex_VecVecs.VecVec(1..nbr);
    mat : DoblDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : DoblDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;
    t : constant double_double := create(0.1);
    eva1,eva2 : DoblDobl_Complex_Vectors.Vector(1..nbr);
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    serelp,mltelp,speedup,efficiency : duration;
    nbt : integer32 := 2;

    use Ada.Calendar;
    use DoblDobl_Rational_Approximations;

  begin
    put_line(file,"double double precision");
    DoblDobl_Allocate(nbr,numdeg,dendeg,cff,numcff1,numcff2,dencff1,dencff2);
    seristart := Ada.Calendar.Clock;
    Pade_Vector(numdeg,dendeg,cff,numcff1,dencff1,mat,rhs,ipvt,info,false);
    Evaluate(numcff1,dencff1,t,eva1);
    seristop := Ada.Calendar.Clock;
    serelp := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,serelp,1,3); new_line(file); flush(file);
    for k in 1..nbruns loop
      multstart := Ada.Calendar.Clock;
      DoblDobl_Multitasking(nbt,numdeg,dendeg,cff,numcff2,dencff2,t,eva2,false);
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
      if serelp + 1.0 /= 1.0 then
        speedup := serelp/mltelp;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
      else
        speedup := 0.0; efficiency := 0.0;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mltelp,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3);
      put(file," : "); duration_io.put(file,efficiency,2,2);
      new_line(file); flush(file);
      nbt := nbt + inc;
    end loop;
  end DoblDobl_Benchmark;

  procedure QuadDobl_Benchmark
              ( file : in file_type;
                nbr,numdeg,dendeg,nbruns,inc : in integer32 ) is

    cff : QuadDobl_Complex_VecVecs.VecVec(1..nbr);
    numcff1,numcff2 : QuadDobl_Complex_VecVecs.VecVec(1..nbr);
    dencff1,dencff2 : QuadDobl_Complex_VecVecs.VecVec(1..nbr);
    mat : QuadDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : QuadDobl_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;
    t : constant quad_double := create(0.1);
    eva1,eva2 : QuadDobl_Complex_Vectors.Vector(1..nbr);
    seristart,seristop,multstart,multstop : Ada.Calendar.Time;
    serelp,mltelp,speedup,efficiency : duration;
    nbt : integer32 := 2;

    use Ada.Calendar;
    use QuadDobl_Rational_Approximations;

  begin
    put_line(file,"quad double precision");
    QuadDobl_Allocate(nbr,numdeg,dendeg,cff,numcff1,numcff2,dencff1,dencff2);
    seristart := Ada.Calendar.Clock;
    Pade_Vector(numdeg,dendeg,cff,numcff1,dencff1,mat,rhs,ipvt,info,false);
    Evaluate(numcff1,dencff1,t,eva1);
    seristop := Ada.Calendar.Clock;
    serelp := seristop - seristart;
    put(file,"  1 : ");
    duration_io.put(file,serelp,1,3); new_line(file); flush(file);
    for k in 1..nbruns loop
      multstart := Ada.Calendar.Clock;
      QuadDobl_Multitasking(nbt,numdeg,dendeg,cff,numcff2,dencff2,t,eva2,false);
      multstop := Ada.Calendar.Clock;
      mltelp := multstop - multstart;
      if serelp + 1.0 /= 1.0 then
        speedup := serelp/mltelp;
        efficiency := speedup/duration(nbt);
        efficiency := duration(100)*efficiency;
      else
        speedup := 0.0; efficiency := 0.0;
      end if;
      put(file,nbt,3);
      put(file," : "); duration_io.put(file,mltelp,1,3);
      put(file," : "); duration_io.put(file,speedup,1,3);
      put(file," : "); duration_io.put(file,efficiency,2,2);
      new_line(file); flush(file);
      nbt := nbt + inc;
    end loop;
  end QuadDobl_Benchmark;

  procedure Benchmark ( nbr,numdeg,dendeg : in integer32 ) is

    file : file_type;
    nbruns,inc : integer32 := 0;
 
  begin
    new_line;
    put("Give the number of multitasked runs : "); get(nbruns);
    put("Give the increment on the tasks : "); get(inc);
    skip_line;
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    put(file,"numerator degree : "); put(file,numdeg,1);
    put(file,"  denominator degree : "); put(file,dendeg,1);
    new_line(file);
    Standard_Benchmark(file,nbr,numdeg,dendeg,nbruns,inc);
    DoblDobl_Benchmark(file,nbr,numdeg,dendeg,nbruns,inc);
    QuadDobl_Benchmark(file,nbr,numdeg,dendeg,nbruns,inc);
  end Benchmark;

  procedure Main is

    numdeg,dendeg,nbr : integer32 := 0;
    ans,precision : character;

  begin
    new_line;
    put_line("Testing the construction of rational approximations ...");
    new_line;
    put("Give the degree of the numerator : "); get(numdeg);
    put("Give the degree of the denominator : "); get(dendeg);
    put("Give the number of components : "); get(nbr);
    new_line;
    put("Benchmark for all precisions ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Benchmark(nbr,numdeg,dendeg);
    else
      new_line;
      put_line("MENU for the working precision :");
      put_line("  0. double precision");
      put_line("  1. double double precision");
      put_line("  2. quad double precision");
      put("Type 0, 1, or 2 to select the precision : ");
      Ask_Alternative(precision,"012");
      case precision is
        when '0' => Standard_Test(nbr,numdeg,dendeg);
        when '1' => DoblDobl_Test(nbr,numdeg,dendeg);
        when '2' => QuadDobl_Test(nbr,numdeg,dendeg);
        when others => null;
      end case;
    end if;
  end Main;

end Test_mtPade_Approximations;
