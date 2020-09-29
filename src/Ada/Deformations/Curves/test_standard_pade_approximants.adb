with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Mathematical_Functions;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Rational_Approximations;
with Standard_Complex_Series_Vectors_io;
with Standard_CSeries_Vector_Functions;
with Standard_Homotopy;
with Homotopy_Series_Readers;
with Standard_Pade_Approximants_io;
with Homotopy_Pade_Approximants;

package body Test_Standard_Pade_Approximants is

  function standard_log_series
             ( dim : integer32 ) return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Numbers;

    res : Standard_Complex_Vectors.Vector(0..dim);

  begin
    res(0) := Create(0.0);
    for k in 1..dim loop
      res(k) := Create(1.0/double_float(k));
      if k mod 2 = 0
       then Min(res(k));
      end if;
    end loop;
    return res;
  end standard_log_series;

  function standard_invfactorial
             ( n : integer32 )
             return Standard_Complex_Numbers.Complex_Number is

    use Standard_Complex_Numbers;

    res : Complex_Number;
    fac : integer32 := 1;
    invfac : double_float;

  begin
    for k in 2..n loop
      fac := fac*k;
    end loop;
    invfac := 1.0/double_float(fac);
    res := Create(invfac);
    return res;
  end standard_invfactorial;

  function standard_exp_series
             ( dim : integer32 ) return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(0..dim);

  begin
    for k in res'range loop
      res(k) := standard_invfactorial(k);
    end loop;
    return res;
  end standard_exp_series;

  function standard_sin_series
             ( dim : integer32 ) return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Numbers;

    res : Standard_Complex_Vectors.Vector(0..dim);
    plus : boolean := true;

  begin
    for k in res'range loop
      if k mod 2 = 0 then
        res(k) := Create(0.0);
      else
        res(k) := standard_invfactorial(k);
        if not plus
         then Min(res(k)); 
        end if;
        plus := not plus;
      end if;
    end loop;
    return res;
  end standard_sin_series;

  function standard_cos_series
             ( dim : integer32 ) return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Numbers;

    res : Standard_Complex_Vectors.Vector(0..dim);
    plus : boolean := true;

  begin
    for k in res'range loop
      if k mod 2 = 1 then
        res(k) := Create(0.0);
      else
        res(k) := standard_invfactorial(k);
        if not plus
         then Min(res(k)); 
        end if;
        plus := not plus;
      end if;
    end loop;
    return res;
  end standard_cos_series;

  procedure Standard_Line_Test ( numdeg,dendeg : in integer32 ) is

    use Standard_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : Standard_Complex_Vectors.Vector(0..dim);
    num : Standard_Complex_Vectors.Vector(0..numdeg);
    den : Standard_Complex_Vectors.Vector(0..dendeg);
    info : integer32;

  begin
    cff := (0..dim => Create(0.0));
    cff(1) := Create(1.0);
    Standard_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
    end if;
  end Standard_Line_Test;

  procedure Standard_log_Test ( numdeg,dendeg : in integer32 ) is

    use Standard_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant Standard_Complex_Vectors.Vector(0..dim)
        := standard_log_series(dim);
    num : Standard_Complex_Vectors.Vector(0..numdeg);
    den : Standard_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    pnt : constant Complex_Number := Create(0.1);
    eva : Complex_Number;
    chkpnt : constant double_float
           := Standard_Mathematical_Functions.LN(1.1);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    Standard_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := Standard_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 1.1      :"); put(eva); new_line;
      put("The value of log(1.1) :"); put(chkpnt); new_line;
    end if;
  end Standard_log_Test;

  procedure Standard_sin_Test ( numdeg,dendeg : in integer32 ) is

    use Standard_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant Standard_Complex_Vectors.Vector(0..dim)
        := standard_sin_series(dim);
    num : Standard_Complex_Vectors.Vector(0..numdeg);
    den : Standard_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    pnt : constant Complex_Number := Create(0.1);
    eva : Complex_Number;
    chkpnt : constant double_float
           := Standard_Mathematical_Functions.SIN(0.1);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    Standard_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := Standard_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      :"); put(eva); new_line;
      put("The value of sin(0.1) :"); put(chkpnt); new_line;
    end if;
  end Standard_Sin_Test;

  procedure Standard_exp_Test ( numdeg,dendeg : in integer32 ) is

    use Standard_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant Standard_Complex_Vectors.Vector(0..dim)
        := standard_exp_series(dim);
    num : Standard_Complex_Vectors.Vector(0..numdeg);
    den : Standard_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    pnt : constant Complex_Number := Create(0.1);
    eva : Complex_Number;
    chkpnt : constant double_float
           := Standard_Mathematical_Functions.EXP(0.1);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    Standard_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := Standard_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      :"); put(eva); new_line;
      put("The value of exp(0.1) :"); put(chkpnt); new_line;
    end if;
  end Standard_exp_Test;

  procedure Standard_cos_Test ( numdeg,dendeg : in integer32 ) is

    use Standard_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant Standard_Complex_Vectors.Vector(0..dim)
        := standard_cos_series(dim);
    num : Standard_Complex_Vectors.Vector(0..numdeg);
    den : Standard_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    pnt : constant Complex_Number := Create(0.1);
    eva : Complex_Number;
    chkpnt : constant double_float
           := Standard_Mathematical_Functions.COS(0.1);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    Standard_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := Standard_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      :"); put(eva); new_line;
      put("The value of cos(0.1) :"); put(chkpnt); new_line;
    end if;
  end Standard_cos_Test;

  procedure Standard_Pade_Approximation
              ( nbequ,nbsteps : in integer32;
                srv : in Standard_Complex_Series_Vectors.Vector;
                pv : in Standard_Pade_Approximants.Pade_Vector ) is

    arg : double_float;
    pnt : Standard_Complex_Numbers.Complex_Number;
    value : Standard_Complex_Vectors.Vector(1..nbequ);
    approx : Standard_Complex_Vectors.Vector(pv'range);
    valsrv : Standard_Complex_Vectors.Vector(srv'range);

  begin
    for k in 1..nbsteps loop
      arg := double_float(k)*0.04;
      pnt := Standard_Complex_Numbers.Create(arg);
      valsrv := Standard_CSeries_Vector_Functions.Eval(srv,pnt);
      put_line("The value of the series approximation :");
      put_line(valsrv);
      approx := Standard_Pade_Approximants.Eval(pv,pnt);
      put_line("The value of the rational approximation :");
      put_line(approx);
      value := Standard_Homotopy.Eval(valsrv,pnt);
      put("Series approximation evaluated at "); put(arg,3); put_line(" :");
      put_line(value);
      value := Standard_Homotopy.Eval(approx,pnt);
      put("Pade approximation evaluated at "); put(arg,3); put_line(" :");
      put_line(value);
    end loop;
  end Standard_Pade_Approximation;

  procedure Standard_Test_Homotopy is

    use Standard_Complex_Polynomials;

    start,target : Standard_Complex_Poly_Systems.Poly_Sys(1..1);
    startpol,targetpol : Poly;
    trm : Term;
    tpow : constant natural32 := 1;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Complex_Numbers.Create(1.0);

  begin
    trm.cf := Standard_Complex_Numbers.Create(1.0);
    trm.dg := new Standard_Natural_Vectors.Vector'(1..1 => 0);
    trm.dg(1) := 2;
    startpol := Create(trm);
    trm.dg(1) := 0;
    Sub(startpol,trm);
    start(1) := startpol;
    trm.dg(1) := 2;
    trm.cf := Standard_Complex_Numbers.Create(3.0);
    targetpol := Create(trm);
    trm.dg(1) := 0;
    trm.cf := Standard_Complex_Numbers.Create(1.5);
    Sub(targetpol,trm);
    target(1) := targetpol;
    Standard_Homotopy.Create(target,start,tpow,gamma);
  end Standard_Test_Homotopy;

  procedure Standard_Test_Start_Solutions
              ( sols : out Standard_Complex_Solutions.Solution_List ) is

    sol : Standard_Complex_Solutions.Solution(1);
    two,one : Standard_Complex_Solutions.Link_to_Solution;

  begin
    sol.t := Standard_Complex_Numbers.Create(0.0);
    sol.m := 1;
    sol.v(1) := Standard_Complex_Numbers.Create(-1.0);
    sol.err := 0.0; sol.rco := 1.0; sol.res := 0.0; 
    two := new Standard_Complex_Solutions.Solution'(sol);
    Standard_Complex_Solutions.Construct(two,sols);
    sol.v(1) := Standard_Complex_Numbers.Create(1.0);
    one := new Standard_Complex_Solutions.Solution'(sol);
    Standard_Complex_Solutions.Construct(one,sols);
  end Standard_Test_Start_Solutions;

  procedure Standard_Test_Case
              ( sols : out Standard_Complex_Solutions.Solution_List ) is
  begin
    Standard_Test_Homotopy;
    Standard_Test_Start_Solutions(sols);
  end Standard_Test_Case;

  procedure Standard_Pade_Homotopy
              ( numdeg,dendeg,nbeq,nbsteps : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List ) is

    lnk : constant Standard_Complex_Solutions.Link_to_Solution
        := Standard_Complex_Solutions.Head_Of(sols);
    sol : constant Standard_Complex_Solutions.Solution := lnk.all;
    nbt : constant natural32 := natural32(numdeg+dendeg+1);
    nit : constant natural32 := 4*nbt;
    srv : Standard_Complex_Series_Vectors.Vector(sol.v'range);
    eva : Standard_Complex_Series_Vectors.Vector(1..nbeq);
    pv : Standard_Pade_Approximants.Pade_Vector(srv'range);

  begin
    Homotopy_Pade_Approximants.Standard_Pade_Approximant
      (sol.v,nbeq+1,nbeq,numdeg,dendeg,nit,srv,eva,pv);
    put_line("The solution series :");
    Standard_Complex_Series_Vectors_io.put(srv);
    put_line("The evaluated solution series :");
    Standard_Complex_Series_Vectors_io.put(eva);
    Standard_Pade_Approximation(nbeq,nbsteps,srv,pv);
    put_line("The Pade approximant :");
    for i in pv'range loop
      put_line(Standard_Pade_Approximants_io.Write(pv(i)));
    end loop;
    Standard_Pade_Approximants.Clear(pv);
  end Standard_Pade_Homotopy;

  procedure Standard_Homotopy_Test ( numdeg,dendeg : in integer32 ) is

    nbeq : integer32;
    nbsteps : integer32 := 5;
    sols : Standard_Complex_Solutions.Solution_List;
    ans : character;

  begin
    new_line;
    put("Run an example test case ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    if ans = 'y' then
      Standard_Test_Case(sols); nbeq := 1; nbsteps := 50;
    else
      put("Random gamma ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
       -- Homotopy_Series_Readers.Standard_Reader(nbeq,sols,tpow=>1);
        Homotopy_Series_Readers.Standard_Reader(nbeq,sols);
      else
        declare
          gamma : constant Standard_Complex_Numbers.Complex_Number
                := Standard_Complex_Numbers.Create(1.0);
        begin
         -- Homotopy_Series_Readers.Standard_Reader(nbeq,sols,1,gamma);
          Homotopy_Series_Readers.Standard_Reader(nbeq,sols,gamma);
        end;
      end if;
      new_line;
      put("Give the number of steps : "); get(nbsteps);
      new_line;
    end if;
    Standard_Pade_Homotopy(numdeg,dendeg,nbeq,nbsteps,sols);
  end Standard_Homotopy_Test;

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
      when '0' => Standard_Line_Test(degnum,degden);
      when '1' => Standard_log_Test(degnum,degden);
      when '2' => Standard_exp_Test(degnum,degden);
      when '3' => Standard_sin_Test(degnum,degden);
      when '4' => Standard_cos_Test(degnum,degden);
      when '5' => Standard_Homotopy_Test(degnum,degden);
      when others => null;
    end case;
  end Main;

end Test_Standard_Pade_Approximants;
