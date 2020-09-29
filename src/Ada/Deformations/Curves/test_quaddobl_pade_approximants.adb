with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with QuadDobl_Mathematical_Functions;
with Standard_Natural_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Rational_Approximations;
with QuadDobl_Complex_Series_Vectors_io;
with QuadDobl_CSeries_Vector_Functions;
with QuadDobl_Homotopy;
with Homotopy_Series_Readers;
with QuadDobl_Pade_Approximants_io;
with Homotopy_Pade_Approximants;

package body Test_QuadDobl_Pade_Approximants is

  function quaddobl_log_series
             ( dim : integer32 ) return QuadDobl_Complex_Vectors.Vector is

    use QuadDobl_Complex_Numbers;

    res : QuadDobl_Complex_Vectors.Vector(0..dim);
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);
    val : quad_double;

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
  end quaddobl_log_series;

  function quaddobl_invfactorial
             ( n : integer32 )
             return QuadDobl_Complex_Numbers.Complex_Number is

    use QuadDobl_Complex_Numbers;

    res : Complex_Number;
    fac : integer32 := 1;
    dd_fac : quad_double;
    invfac : quad_double;

  begin
    for k in 2..n loop
      fac := fac*k;
    end loop;
    dd_fac := create(fac);
    invfac := 1.0/dd_fac;
    res := Create(invfac);
    return res;
  end quaddobl_invfactorial;

  function quaddobl_exp_series
             ( dim : integer32 ) return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(0..dim);

  begin
    for k in res'range loop
      res(k) := quaddobl_invfactorial(k);
    end loop;
    return res;
  end quaddobl_exp_series;

  function quaddobl_sin_series
             ( dim : integer32 ) return QuadDobl_Complex_Vectors.Vector is

    use QuadDobl_Complex_Numbers;

    res : QuadDobl_Complex_Vectors.Vector(0..dim);
    zero : constant quad_double := create(0.0);
    plus : boolean := true;

  begin
    for k in res'range loop
      if k mod 2 = 0 then
        res(k) := Create(zero);
      else
        res(k) := quaddobl_invfactorial(k);
        if not plus
         then Min(res(k)); 
        end if;
        plus := not plus;
      end if;
    end loop;
    return res;
  end quaddobl_sin_series;

  function quaddobl_cos_series
             ( dim : integer32 ) return QuadDobl_Complex_Vectors.Vector is

    use QuadDobl_Complex_Numbers;

    res : QuadDobl_Complex_Vectors.Vector(0..dim);
    zero : constant quad_double := create(0.0);
    plus : boolean := true;

  begin
    for k in res'range loop
      if k mod 2 = 1 then
        res(k) := Create(zero);
      else
        res(k) := quaddobl_invfactorial(k);
        if not plus
         then Min(res(k)); 
        end if;
        plus := not plus;
      end if;
    end loop;
    return res;
  end quaddobl_cos_series;

  procedure QuadDobl_Line_Test ( numdeg,dendeg : in integer32 ) is

    use QuadDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : QuadDobl_Complex_Vectors.Vector(0..dim);
    num : QuadDobl_Complex_Vectors.Vector(0..numdeg);
    den : QuadDobl_Complex_Vectors.Vector(0..dendeg);
    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);
    info : integer32;

  begin
    cff := (0..dim => Create(zero));
    cff(1) := Create(one);
    QuadDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
    end if;
  end QuadDobl_Line_Test;

  procedure QuadDobl_log_Test ( numdeg,dendeg : in integer32 ) is

    use QuadDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant QuadDobl_Complex_Vectors.Vector(0..dim)
        := QuadDobl_log_series(dim);
    num : QuadDobl_Complex_Vectors.Vector(0..numdeg);
    den : QuadDobl_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    nbr : constant quad_double := create(0.1);
    pnt : constant Complex_Number := Create(nbr);
    eva : Complex_Number;
    arg : constant quad_double := create(1.1);
    chkpnt : constant quad_double := Quad_Double_Numbers.log(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    QuadDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := QuadDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 1.1      : "); put(eva); new_line;
      put("The value of log(1.1) : "); put(chkpnt); new_line;
    end if;
  end QuadDobl_log_Test;

  procedure QuadDobl_sin_Test ( numdeg,dendeg : in integer32 ) is

    use QuadDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant QuadDobl_Complex_Vectors.Vector(0..dim)
        := quaddobl_sin_series(dim);
    num : QuadDobl_Complex_Vectors.Vector(0..numdeg);
    den : QuadDobl_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    arg : constant quad_double := create(0.1);
    pnt : constant Complex_Number := Create(arg);
    eva : Complex_Number;
    chkpnt : constant quad_double
           := QuadDobl_Mathematical_Functions.SIN(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    QuadDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := QuadDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      : "); put(eva); new_line;
      put("The value of sin(0.1) : "); put(chkpnt); new_line;
    end if;
  end QuadDobl_Sin_Test;

  procedure QuadDobl_exp_Test ( numdeg,dendeg : in integer32 ) is

    use QuadDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant QuadDobl_Complex_Vectors.Vector(0..dim)
        := quaddobl_exp_series(dim);
    num : QuadDobl_Complex_Vectors.Vector(0..numdeg);
    den : QuadDobl_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    arg : constant quad_double := create(0.1);
    pnt : constant Complex_Number := Create(arg);
    eva : Complex_Number;
    chkpnt : constant quad_double := Quad_Double_Numbers.exp(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    QuadDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := QuadDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      : "); put(eva); new_line;
      put("The value of exp(0.1) : "); put(chkpnt); new_line;
    end if;
  end QuadDobl_exp_Test;

  procedure QuadDobl_cos_Test ( numdeg,dendeg : in integer32 ) is

    use QuadDobl_Complex_Numbers;

    dim : constant integer32 := numdeg + dendeg;
    cff : constant QuadDobl_Complex_Vectors.Vector(0..dim)
        := quaddobl_cos_series(dim);
    num : QuadDobl_Complex_Vectors.Vector(0..numdeg);
    den : QuadDobl_Complex_Vectors.Vector(0..dendeg);
    info : integer32;
    arg : constant quad_double := create(0.1);
    pnt : constant Complex_Number := Create(arg);
    eva : Complex_Number;
    chkpnt : constant quad_double
           := QuadDobl_Mathematical_Functions.COS(arg);

  begin
    put_line("The coefficient vector of the series :"); put_line(cff);
    QuadDobl_Rational_Approximations.Pade(numdeg,dendeg,cff,num,den,info);
    if info /= 0 then
      put("info : "); put(info,1); put_line(" error!");
    else
      put_line("The coefficients of the numerator :"); put_line(num);
      put_line("The coefficients of the denominator :"); put_line(den);
      eva := QuadDobl_Rational_Approximations.Evaluate(num,den,pnt);
      put("The value at 0.1      : "); put(eva); new_line;
      put("The value of cos(0.1) : "); put(chkpnt); new_line;
    end if;
  end QuadDobl_cos_Test;

  procedure QuadDobl_Pade_Approximation
              ( nbequ,nbsteps : in integer32;
                srv : in QuadDobl_Complex_Series_Vectors.Vector;
                pv : in QuadDobl_Pade_Approximants.Pade_Vector ) is

    dd_arg : double_double;
    qd_arg : quad_double;
    pnt : QuadDobl_Complex_Numbers.Complex_Number;
    value : QuadDobl_Complex_Vectors.Vector(1..nbequ);
    approx : QuadDobl_Complex_Vectors.Vector(pv'range);
    valsrv : QuadDobl_Complex_Vectors.Vector(srv'range);

  begin
    for k in 1..nbsteps loop
      dd_arg := Double_Double_Numbers.create(double_float(k)*0.04);
      qd_arg := Quad_Double_Numbers.create(dd_arg);
      pnt := QuadDobl_Complex_Numbers.Create(qd_arg);
      valsrv := QuadDobl_CSeries_Vector_Functions.Eval(srv,pnt);
      put_line("The value of the series approximation :");
      put_line(valsrv);
      approx := QuadDobl_Pade_Approximants.Eval(pv,pnt);
      put_line("The value of the rational approximation :");
      put_line(approx);
      value := QuadDobl_Homotopy.Eval(valsrv,pnt);
      put("Series approximation evaluated at "); put(qd_arg,3); put_line(" :");
      put_line(value);
      value := QuadDobl_Homotopy.Eval(approx,pnt);
      put("Pade approximation evaluated at "); put(qd_arg,3); put_line(" :");
      put_line(value);
    end loop;
  end QuadDobl_Pade_Approximation;

  procedure QuadDobl_Test_Homotopy is

    use QuadDobl_Complex_Polynomials;

    start,target : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..1);
    startpol,targetpol : Poly;
    trm : Term;
    tpow : constant natural32 := 1;
    nbr : quad_double := create(1.0);
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Complex_Numbers.Create(nbr);

  begin
    trm.cf := QuadDobl_Complex_Numbers.Create(nbr);
    trm.dg := new Standard_Natural_Vectors.Vector'(1..1 => 0);
    trm.dg(1) := 2;
    startpol := Create(trm);
    trm.dg(1) := 0;
    Sub(startpol,trm);
    start(1) := startpol;
    trm.dg(1) := 2;
    nbr := create(3.0);
    trm.cf := QuadDobl_Complex_Numbers.Create(nbr);
    targetpol := Create(trm);
    trm.dg(1) := 0;
    nbr := create(1.5);
    trm.cf := QuadDobl_Complex_Numbers.Create(nbr);
    Sub(targetpol,trm);
    target(1) := targetpol;
    QuadDobl_Homotopy.Create(target,start,tpow,gamma);
  end QuadDobl_Test_Homotopy;

  procedure QuadDobl_Test_Start_Solutions
              ( sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    sol : QuadDobl_Complex_Solutions.Solution(1);
    nbr : quad_double;
    two,one : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    nbr := create(0.0);
    sol.t := QuadDobl_Complex_Numbers.Create(nbr);
    sol.m := 1;
    nbr := create(-1.0);
    sol.v(1) := QuadDobl_Complex_Numbers.Create(nbr);
    sol.err := create(0.0);
    sol.rco := create(1.0);
    sol.res := create(0.0); 
    two := new QuadDobl_Complex_Solutions.Solution'(sol);
    QuadDobl_Complex_Solutions.Construct(two,sols);
    nbr := create(1.0);
    sol.v(1) := QuadDobl_Complex_Numbers.Create(nbr);
    one := new QuadDobl_Complex_Solutions.Solution'(sol);
    QuadDobl_Complex_Solutions.Construct(one,sols);
  end QuadDobl_Test_Start_Solutions;

  procedure QuadDobl_Test_Case
              ( sols : out QuadDobl_Complex_Solutions.Solution_List ) is
  begin
    QuadDobl_Test_Homotopy;
    QuadDobl_Test_Start_Solutions(sols);
  end QuadDobl_Test_Case;

  procedure QuadDobl_Pade_Homotopy
              ( numdeg,dendeg,nbeq,nbsteps : in integer32;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

    lnk : constant QuadDobl_Complex_Solutions.Link_to_Solution
        := QuadDobl_Complex_Solutions.Head_Of(sols);
    sol : constant QuadDobl_Complex_Solutions.Solution := lnk.all;
    nbt : constant natural32 := natural32(numdeg+dendeg+1);
    nit : constant natural32 := 4*nbt;
    srv : QuadDobl_Complex_Series_Vectors.Vector(sol.v'range);
    eva : QuadDobl_Complex_Series_Vectors.Vector(1..nbeq);
    pv : QuadDobl_Pade_Approximants.Pade_Vector(srv'range);

  begin
    Homotopy_Pade_Approximants.QuadDobl_Pade_Approximant
      (sol.v,nbeq+1,nbeq,numdeg,dendeg,nit,srv,eva,pv);
    put_line("The solution series :");
    QuadDobl_Complex_Series_Vectors_io.put(srv);
    put_line("The evaluated solution series :");
    QuadDobl_Complex_Series_Vectors_io.put(eva);
    QuadDobl_Pade_Approximation(nbeq,nbsteps,srv,pv);
    put_line("The Pade approximant :");
    for i in pv'range loop
      put_line(QuadDobl_Pade_Approximants_io.Write(pv(i)));
    end loop;
    QuadDobl_Pade_Approximants.Clear(pv);
  end QuadDobl_Pade_Homotopy;

  procedure QuadDobl_Homotopy_Test ( numdeg,dendeg : in integer32 ) is

    nbeq : integer32;
    nbsteps : integer32 := 5;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    ans : character;

  begin
    new_line;
    put("Run an example test case ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    if ans = 'y' then
      QuadDobl_Test_Case(sols); nbeq := 1; nbsteps := 50;
    else
      put("Random gamma ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
       -- Homotopy_Series_Readers.QuadDobl_Reader(nbeq,sols,tpow=>1);
        Homotopy_Series_Readers.QuadDobl_Reader(nbeq,sols);
      else
        declare
          one : constant quad_double := create(1.0);
          gamma : constant QuadDobl_Complex_Numbers.Complex_Number
                := QuadDobl_Complex_Numbers.Create(one);
        begin
         -- Homotopy_Series_Readers.QuadDobl_Reader(nbeq,sols,1,gamma);
          Homotopy_Series_Readers.QuadDobl_Reader(nbeq,sols,gamma);
        end;
      end if;
      new_line;
      put("Give the number of steps : "); get(nbsteps);
      new_line;
    end if;
    QuadDobl_Pade_Homotopy(numdeg,dendeg,nbeq,nbsteps,sols);
  end QuadDobl_Homotopy_Test;

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
      when '0' => QuadDobl_Line_Test(degnum,degden);
      when '1' => QuadDobl_log_Test(degnum,degden);
      when '2' => QuadDobl_exp_Test(degnum,degden);
      when '3' => QuadDobl_sin_Test(degnum,degden);
      when '4' => QuadDobl_cos_Test(degnum,degden);
      when '5' => QuadDobl_Homotopy_Test(degnum,degden);
      when others => null;
    end case;
  end Main;

end Test_QuadDobl_Pade_Approximants;
