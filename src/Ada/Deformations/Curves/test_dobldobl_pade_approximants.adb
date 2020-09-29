with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with DoblDobl_Mathematical_Functions;
with Standard_Natural_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Rational_Approximations;
with DoblDobl_Complex_Series_Vectors_io;
with DoblDobl_CSeries_Vector_Functions;
with DoblDobl_Homotopy;
with Homotopy_Series_Readers;
with DoblDobl_Pade_Approximants_io;
with Homotopy_Pade_Approximants;

package body Test_DoblDobl_Pade_Approximants is

  function dobldobl_log_series
             ( dim : integer32 ) return DoblDobl_Complex_Vectors.Vector is

    use DoblDobl_Complex_Numbers;

    res : DoblDobl_Complex_Vectors.Vector(0..dim);
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);
    val : double_double;

  begin
    res(0) := Create(zero);
    for k in 1..dim loop
      val := create(k);
      val := one/val;
      res(k) := Create(val);
      if k mod 2 = 0
       then Min(res(k));
      end if;
    end loop;
    return res;
  end dobldobl_log_series;

  function dobldobl_invfactorial
             ( n : integer32 )
             return DoblDobl_Complex_Numbers.Complex_Number is

    use DoblDobl_Complex_Numbers;

    res : Complex_Number;
    fac : integer32 := 1;
    dd_fac : double_double;
    invfac : double_double;

  begin
    for k in 2..n loop
      fac := fac*k;
    end loop;
    dd_fac := create(fac);
    invfac := 1.0/dd_fac;
    res := Create(invfac);
    return res;
  end dobldobl_invfactorial;

  function dobldobl_exp_series
             ( dim : integer32 ) return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(0..dim);

  begin
    for k in res'range loop
      res(k) := dobldobl_invfactorial(k);
    end loop;
    return res;
  end dobldobl_exp_series;

  function dobldobl_sin_series
             ( dim : integer32 ) return DoblDobl_Complex_Vectors.Vector is

    use DoblDobl_Complex_Numbers;

    res : DoblDobl_Complex_Vectors.Vector(0..dim);
    zero : constant double_double := create(0.0);
    plus : boolean := true;

  begin
    for k in res'range loop
      if k mod 2 = 0 then
        res(k) := Create(zero);
      else
        res(k) := dobldobl_invfactorial(k);
        if not plus
         then Min(res(k)); 
        end if;
        plus := not plus;
      end if;
    end loop;
    return res;
  end dobldobl_sin_series;

  function dobldobl_cos_series
             ( dim : integer32 ) return DoblDobl_Complex_Vectors.Vector is

    use DoblDobl_Complex_Numbers;

    res : DoblDobl_Complex_Vectors.Vector(0..dim);
    zero : constant double_double := create(0.0);
    plus : boolean := true;

  begin
    for k in res'range loop
      if k mod 2 = 1 then
        res(k) := Create(zero);
      else
        res(k) := dobldobl_invfactorial(k);
        if not plus
         then Min(res(k)); 
        end if;
        plus := not plus;
      end if;
    end loop;
    return res;
  end dobldobl_cos_series;

  procedure DoblDobl_Line_Test ( numdeg,dendeg : in integer32 ) is

    use DoblDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : DoblDobl_Complex_Vectors.Vector(0..dim);
    num : DoblDobl_Complex_Vectors.Vector(0..numdeg);
    den : DoblDobl_Complex_Vectors.Vector(0..dendeg);
    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);
    info : integer32;

  begin
    cff := (0..dim => Create(zero));
    cff(1) := Create(one);
    DoblDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
    end if;
  end DoblDobl_Line_Test;

  procedure DoblDobl_log_Test ( numdeg,dendeg : in integer32 ) is

    use DoblDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant DoblDobl_Complex_Vectors.Vector(0..dim)
        := DoblDobl_log_series(dim);
    num : DoblDobl_Complex_Vectors.Vector(0..numdeg);
    den : DoblDobl_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    nbr : constant double_double := create(0.1);
    pnt : constant Complex_Number := Create(nbr);
    eva : Complex_Number;
    arg : constant double_double := create(1.1);
    chkpnt : constant double_double := Double_Double_Numbers.log(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    DoblDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := DoblDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 1.1      : "); put(eva); new_line;
      put("The value of log(1.1) : "); put(chkpnt); new_line;
    end if;
  end DoblDobl_log_Test;

  procedure DoblDobl_sin_Test ( numdeg,dendeg : in integer32 ) is

    use DoblDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant DoblDobl_Complex_Vectors.Vector(0..dim)
        := dobldobl_sin_series(dim);
    num : DoblDobl_Complex_Vectors.Vector(0..numdeg);
    den : DoblDobl_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    arg : constant double_double := create(0.1);
    pnt : constant Complex_Number := Create(arg);
    eva : Complex_Number;
    chkpnt : constant double_double
           := DoblDobl_Mathematical_Functions.SIN(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    DoblDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := DoblDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      : "); put(eva); new_line;
      put("The value of sin(0.1) : "); put(chkpnt); new_line;
    end if;
  end DoblDobl_sin_Test;

  procedure DoblDobl_exp_Test ( numdeg,dendeg : in integer32 ) is

    use DoblDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant DoblDobl_Complex_Vectors.Vector(0..dim)
        := dobldobl_exp_series(dim);
    num : DoblDobl_Complex_Vectors.Vector(0..numdeg);
    den : DoblDobl_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    arg : constant double_double := create(0.1);
    pnt : constant Complex_Number := Create(arg);
    eva : Complex_Number;
    chkpnt : constant double_double := Double_Double_Numbers.exp(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    DoblDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := DoblDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      : "); put(eva); new_line;
      put("The value of exp(0.1) : "); put(chkpnt); new_line;
    end if;
  end DoblDobl_exp_Test;

  procedure DoblDobl_cos_Test ( numdeg,dendeg : in integer32 ) is

    use DoblDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant DoblDobl_Complex_Vectors.Vector(0..dim)
        := dobldobl_cos_series(dim);
    num : DoblDobl_Complex_Vectors.Vector(0..numdeg);
    den : DoblDobl_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    arg : constant double_double := create(0.1);
    pnt : constant Complex_Number := Create(arg);
    eva : Complex_Number;
    chkpnt : constant double_double
           := DoblDobl_Mathematical_Functions.COS(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    DoblDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := DoblDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      : "); put(eva); new_line;
      put("The value of cos(0.1) : "); put(chkpnt); new_line;
    end if;
  end DoblDobl_cos_Test;

  procedure DoblDobl_Pade_Approximation
              ( nbequ,nbsteps : in integer32;
                srv : in DoblDobl_Complex_Series_Vectors.Vector;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector ) is

    arg : double_double;
    pnt : DoblDobl_Complex_Numbers.Complex_Number;
    value : DoblDobl_Complex_Vectors.Vector(1..nbequ);
    approx : DoblDobl_Complex_Vectors.Vector(pv'range);
    valsrv : DoblDobl_Complex_Vectors.Vector(srv'range);

  begin
    for k in 1..nbsteps loop
      arg := Double_Double_Numbers.create(double_float(k)*0.04);
      pnt := DoblDobl_Complex_Numbers.Create(arg);
      valsrv := DoblDobl_CSeries_Vector_Functions.Eval(srv,pnt);
      put_line("The value of the series approximation :");
      put_line(valsrv);
      approx := DoblDobl_Pade_Approximants.Eval(pv,pnt);
      put_line("The value of the rational approximation :");
      put_line(approx);
      value := DoblDobl_Homotopy.Eval(valsrv,pnt);
      put("Series approximation evaluated at "); put(arg,3); put_line(" :");
      put_line(value);
      value := DoblDobl_Homotopy.Eval(approx,pnt);
      put("Pade approximant evaluated at "); put(arg,3); put_line(" :");
      put_line(value);
    end loop;
  end DoblDobl_Pade_Approximation;

  procedure DoblDobl_Test_Homotopy is

    use DoblDobl_Complex_Polynomials;

    start,target : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..1);
    startpol,targetpol : Poly;
    trm : Term;
    tpow : constant natural32 := 1;
    nbr : double_double := create(1.0);
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Complex_Numbers.Create(nbr);

  begin
    trm.cf := DoblDobl_Complex_Numbers.Create(nbr);
    trm.dg := new Standard_Natural_Vectors.Vector'(1..1 => 0);
    trm.dg(1) := 2;
    startpol := Create(trm);
    trm.dg(1) := 0;
    Sub(startpol,trm);
    start(1) := startpol;
    trm.dg(1) := 2;
    nbr := create(3.0);
    trm.cf := DoblDobl_Complex_Numbers.Create(nbr);
    targetpol := Create(trm);
    trm.dg(1) := 0;
    nbr := create(1.5);
    trm.cf := DoblDobl_Complex_Numbers.Create(nbr);
    Sub(targetpol,trm);
    target(1) := targetpol;
    DoblDobl_Homotopy.Create(target,start,tpow,gamma);
  end DoblDobl_Test_Homotopy;

  procedure DoblDobl_Test_Start_Solutions
              ( sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    sol : DoblDobl_Complex_Solutions.Solution(1);
    nbr : double_double;
    two,one : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    nbr := create(0.0);
    sol.t := DoblDobl_Complex_Numbers.Create(nbr);
    sol.m := 1;
    nbr := create(-1.0);
    sol.v(1) := DoblDobl_Complex_Numbers.Create(nbr);
    sol.err := create(0.0);
    sol.rco := create(1.0);
    sol.res := create(0.0); 
    two := new DoblDobl_Complex_Solutions.Solution'(sol);
    DoblDobl_Complex_Solutions.Construct(two,sols);
    nbr := create(1.0);
    sol.v(1) := DoblDobl_Complex_Numbers.Create(nbr);
    one := new DoblDobl_Complex_Solutions.Solution'(sol);
    DoblDobl_Complex_Solutions.Construct(one,sols);
  end DoblDobl_Test_Start_Solutions;

  procedure DoblDobl_Test_Case
              ( sols : out DoblDobl_Complex_Solutions.Solution_List ) is

  begin
    DoblDobl_Test_Homotopy;
    DoblDobl_Test_Start_Solutions(sols);
  end DoblDobl_Test_Case;

  procedure DoblDobl_Pade_Homotopy
              ( numdeg,dendeg,nbeq,nbsteps : in integer32;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

    lnk : constant DoblDobl_Complex_Solutions.Link_to_Solution
        := DoblDobl_Complex_Solutions.Head_Of(sols);
    sol : constant DoblDobl_Complex_Solutions.Solution := lnk.all;
    nbt : constant natural32 := natural32(numdeg+dendeg+1);
    nit : constant natural32 := 4*nbt;
    srv : DoblDobl_Complex_Series_Vectors.Vector(sol.v'range);
    eva : DoblDobl_Complex_Series_Vectors.Vector(1..nbeq);
    pv : DoblDobl_Pade_Approximants.Pade_Vector(srv'range);

  begin
    Homotopy_Pade_Approximants.DoblDobl_Pade_Approximant
      (sol.v,nbeq+1,nbeq,numdeg,dendeg,nit,srv,eva,pv);
    put_line("The solution series :");
    DoblDobl_Complex_Series_Vectors_io.put(srv);
    put_line("The evaluated solution series :");
    DoblDobl_Complex_Series_Vectors_io.put(eva);
    DoblDobl_Pade_Approximation(nbeq,nbsteps,srv,pv);
    put_line("The Pade approximant :");
    for i in pv'range loop
      put_line(DoblDobl_Pade_Approximants_io.Write(pv(i)));
    end loop;
    DoblDobl_Pade_Approximants.Clear(pv);
  end DoblDobl_Pade_Homotopy;

  procedure DoblDobl_Homotopy_Test ( numdeg,dendeg : in integer32 ) is

    nbeq : integer32;
    nbsteps : integer32 := 5;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    ans : character;

  begin
    new_line;
    put("Run an example test case ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    if ans = 'y' then
      DoblDobl_Test_Case(sols); nbeq := 1; nbsteps := 50;
    else
      put("Random gamma ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
       -- Homotopy_Series_Readers.DoblDobl_Reader(nbeq,sols,tpow=>1);
        Homotopy_Series_Readers.DoblDobl_Reader(nbeq,sols);
      else
        declare
          one : constant double_double := create(1.0);
          gamma : constant DoblDobl_Complex_Numbers.Complex_Number
                := DoblDobl_Complex_Numbers.Create(one);
        begin
         -- Homotopy_Series_Readers.DoblDobl_Reader(nbeq,sols,1,gamma);
          Homotopy_Series_Readers.DoblDobl_Reader(nbeq,sols,gamma);
        end;
      end if;
      new_line;
      put("Give the number of steps : "); get(nbsteps);
      new_line;
    end if;
    DoblDobl_Pade_Homotopy(numdeg,dendeg,nbeq,nbsteps,sols);
  end DoblDobl_Homotopy_Test;

  procedure Main is

    degnum : integer32 := 0; -- degree numerator
    degden : integer32 := 0; -- degree denominator
    dim : integer32 := 0;    -- degnum + degden
    ans : character;

  begin
    new_line;
    put_line("MENU to test rational approximations for a function ...");
    put_line("  0. simple test on line x = t");
    put_line("  1. series for log(1+x)");
    put_line("  2. series for exp(x)");
    put_line("  3. series for sin(x)");
    put_line("  4. series for cos(x)");
    put_line("  5. artificial parameter homotopy");
    put("Type 0, 1, 2, 3, 4, or 5 to choose the test : ");
    Ask_Alternative(ans,"012345");
    new_line;
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    dim := degnum + degden;
    put("The dimension : "); put(dim,1); new_line;
    case ans is 
      when '0' => DoblDobl_Line_Test(degnum,degden);
      when '1' => DoblDobl_log_Test(degnum,degden);
      when '2' => DoblDobl_exp_Test(degnum,degden);
      when '3' => DoblDobl_sin_Test(degnum,degden);
      when '4' => DoblDobl_cos_Test(degnum,degden);
      when '5' => DoblDobl_Homotopy_Test(degnum,degden);
      when others => null;
    end case;
  end Main;

end Test_DoblDobl_Pade_Approximants;
