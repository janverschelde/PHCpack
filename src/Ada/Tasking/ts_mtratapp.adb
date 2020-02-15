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
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Random_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with DoblDobl_Random_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with QuadDobl_Random_Vectors;
with Evaluation_Differentiation_Errors;  use Evaluation_Differentiation_Errors;
with Standard_Rational_Approximations;
with DoblDobl_Rational_Approximations;
with QuadDobl_Rational_Approximations;
with Multitasked_Pade_Approximations;    use Multitasked_Pade_Approximations;

procedure ts_mtratapp is

-- DESCRIPTION :
--   Development of multitasked algorithms for rational approximations.

  procedure Standard_Test ( nbr,numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Makes a vector of range 1..nbr of random coefficients to construct
  --   rational approximations of the given degrees, in double precision.

  -- ON ENTRY :
  --   nbr       the number of components in the vector of approximants;
  --   numdeg    degree of the numerators;
  --   dendeg    degree of the denominators.

    dim : constant integer32 := numdeg + dendeg;
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
    serelp,mltelp,speedup : duration;

    use Ada.Calendar;
    use Standard_Complex_Numbers;
    use Standard_Rational_Approximations;

  begin
    new_line;
    put_line("Allocating and generating data ...");
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
      put("The speedup : "); duration_io.put(speedup,1,3); new_line;
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

  -- DESCRIPTION :
  --   Makes a vector of range 1..nbr of random coefficients to construct
  --   rational approximations of the given degrees,
  --   in double double precision.

  -- ON ENTRY :
  --   nbr       the number of components in the vector of approximants;
  --   numdeg    degree of the numerators;
  --   dendeg    degree of the denominators.

    dim : constant integer32 := numdeg + dendeg;
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
    serelp,mltelp,speedup : duration;

    use Ada.Calendar;
    use DoblDobl_Complex_Numbers;
    use DoblDobl_Rational_Approximations;

  begin
    new_line;
    put_line("Allocating and generating data ...");
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
      put("The speedup : "); duration_io.put(speedup,1,3); new_line;
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

  -- DESCRIPTION :
  --   Makes a vector of range 1..nbr of random coefficients to construct
  --   rational approximations of the given degrees,
  --   in double double precision.

  -- ON ENTRY :
  --   nbr       the number of components in the vector of approximants;
  --   numdeg    degree of the numerators;
  --   dendeg    degree of the denominators.

    dim : constant integer32 := numdeg + dendeg;
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
    serelp,mltelp,speedup : duration;

    use Ada.Calendar;
    use QuadDobl_Complex_Numbers;
    use QuadDobl_Rational_Approximations;

  begin
    new_line;
    put_line("Allocating and generating data ...");
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
      put("The speedup : "); duration_io.put(speedup,1,3); new_line;
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

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension of the problem
  --   and then launches the test.

    numdeg,dendeg,nbr : integer32 := 0;
    precision : character;

  begin
    new_line;
    put_line("Testing the construction of rational approximations ...");
    new_line;
    put("Give the degree of the numerator : "); get(numdeg);
    put("Give the degree of the denominator : "); get(dendeg);
    put("Give the number of components : "); get(nbr);
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
  end Main;

begin
  Main;
end ts_mtratapp;
