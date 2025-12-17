with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Double_Leading_Evaluations;
with Double_Ordered_Evaluations;
with Random_Laurent_Homotopy;
with Laurent_Homotopy_Derivatives;

package body Test_Newton_Puiseux is

  procedure Evaluate_and_Differentiate
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                zt0 : in Standard_Complex_Vectors.Vector;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
		cjm : out Standard_Complex_Matrices.Matrix;
		ejm : out Standard_Floating_Matrices.Matrix;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := hcf'last;
    tol : constant double_float := 1.0e-12;
    dfc : double_float;

    use Double_Leading_Evaluations;

  begin
    if vrblvl > 0 then
      put_line("-> in Test_Newton_Puiseux.evaluate_and_differentiate ...");
    end if;
    for i in 1..dim loop -- run over all polynomials
      if vrblvl > 0 then
        put("at monomial 1");
        put(" of polynomial "); put(i,1); put_line(" ...");
      end if;
      ydg(i) := hct(i)(1); -- initialize with first monomial
      ycf(i) := hcf(i)(1)*Leading_Coefficient(hdg(i)(1).all,zt0);
      for k in 1..dim loop
        ejm(i,k) := hct(i)(1);
        cjm(i,k) := hcf(i)(1)*Leading_Coefficient(hdg(i)(1).all,zt0,k);
      end loop;
      put(ycf(i)); put(" t^"); put(ydg(i)); new_line;
      for j in 2..hcf(i)'last loop -- run over all monomials
        if vrblvl > 0 then
          put("at monomial "); put(j,1);
          put(" of polynomial "); put(i,1); put_line(" ...");
        end if;
        dfc := abs(ydg(i) - hct(i)(j));
        if dfc < tol then
          ycf(i) := ycf(i) + hcf(i)(j)*Leading_Coefficient(hdg(i)(j).all,zt0);
          for k in 1..dim loop
            cjm(i,k) := cjm(i,k)
              + hcf(i)(j)*Leading_Coefficient(hdg(i)(j).all,zt0,k);
          end loop;
        elsif hct(i)(j) < ydg(i) then
          ydg(i) := hct(i)(j);
          ycf(i) := hcf(i)(j)*Leading_Coefficient(hdg(i)(j).all,zt0);
          for k in 1..dim loop
            ejm(i,k) := hct(i)(j);
            cjm(i,k) := hcf(i)(j)*Leading_Coefficient(hdg(i)(j).all,zt0,k);
          end loop;
        end if;
        put(ycf(i)); put(" t^"); put(ydg(i)); new_line;
      end loop;
    end loop;
  end Evaluate_and_Differentiate;

  procedure Evaluate_All_Monomials
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                zt0 : in Standard_Complex_Vectors.Vector;
                ycf : in Standard_Complex_VecVecs.VecVec;
                ydg : in Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := hcf'last;

    use Double_Leading_Evaluations;

  begin
    if vrblvl > 0 then
      put_line("-> in Test_Newton_Puiseux.evaluate_all_monomials ...");
    end if;
    for i in 1..dim loop -- run over all polynomials
      for j in hcf(i)'range loop -- run over all monomials
        if vrblvl > 0 then
          put("at monomial "); put(j,1);
          put(" of polynomial "); put(i,1); put_line(" ...");
        end if;
        ydg(i)(j) := hct(i)(j);
        ycf(i)(j) := hcf(i)(j)*Leading_Coefficient(hdg(i)(j).all,zt0);
        put(ycf(i)(j)); put(" t^"); put(ydg(i)(j)); new_line;
      end loop;
    end loop;
  end Evaluate_All_Monomials;

  procedure Sum ( cf : in out Standard_Complex_Vectors.Vector;
                  dg : in Standard_Floating_Vectors.Vector ) is

    tol : constant double_float := 1.0e-12;
    dif : double_float;

  begin
    for i in cf'first..cf'last-1 loop
      dif := abs(dg(i) - dg(i+1));
      if dif < tol then
        cf(i) := cf(i) + cf(i+1);
        cf(i+1) := create(0.0);
      end if;
    end loop;
  end Sum;

  procedure Sum ( cf : in Standard_Complex_VecVecs.VecVec;
                  dg : in Standard_Floating_VecVecs.VecVec ) is
  begin
    for i in cf'range loop
      Sum(cf(i).all,dg(i).all);
    end loop;
  end Sum;

  function Positive_Minimum
             ( v : Standard_Floating_Vectors.Vector ) return double_float is

    res : double_float;
    tol : constant double_float := 1.0E-12;
    idx : integer32;

  begin
    for i in v'range loop -- find first positive number
      if v(i) > tol then
        idx := i;
        res := v(i);
        exit;
      end if;
    end loop;
    for i in idx+1..v'last loop
      if v(i) > tol and then v(i) < res
       then res := v(i);
      end if;
    end loop;
    return res;
  end Positive_Minimum;

  function Positive_Minimum
             ( c : Standard_Complex_Vectors.Vector;
               v : Standard_Floating_Vectors.Vector ) return double_float is

    res : double_float;
    tol : constant double_float := 1.0E-12;
    idx : integer32;

  begin
    for i in v'range loop -- find first positive number
      if AbsVal(c(i)) > tol then
        if v(i) > tol then
          idx := i;
          res := v(i);
          exit;
        end if;
      end if;
    end loop;
    for i in idx+1..v'last loop
      if AbsVal(c(i)) > tol then
        if v(i) > tol and then v(i) < res
         then res := v(i);
        end if;
      end if;
    end loop;
    return res;
  end Positive_Minimum;

  function Coefficient ( c : Standard_Complex_Vectors.Vector;
                         e : Standard_Floating_Vectors.Vector;
                         p : double_float ) return Complex_Number is

    tol : constant double_float := 1.0e-12;

  begin
    for i in e'range loop
      if abs(e(i) - p) < tol
       then return c(i);
      end if;
    end loop;
    return create(0.0);
  end Coefficient;

  procedure Leading_Powers_by_Evaluation
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                lcf : in Standard_Complex_Vectors.Vector;
                lpw : in Standard_Floating_Vectors.Vector;
                psm : out Standard_Floating_Vectors.Vector; 
                cfp : out Standard_Complex_Vectors.Vector; 
                vrblvl : in integer32 := 0 ) is

    cfy : Standard_Complex_VecVecs.VecVec(hcf'range);
    dgy : Standard_Floating_VecVecs.VecVec(hct'range);
    err,sumerr : double_float;

  begin
    for i in hcf'range loop
      cfy(i) := new Standard_Complex_Vectors.Vector(hcf(i)'range);
      dgy(i) := new Standard_Floating_Vectors.Vector(hcf(i)'range);
    end loop;
    Evaluate_All_Monomials(hcf,hct,hdg,lcf,cfy,dgy,vrblvl-1);
    if vrblvl > 0 then
      for i in cfy'range loop
        put("all values of polynomial "); put(i,1); put_line(" :");
        for j in cfy(i)'range loop
          put(cfy(i)(j)); put(" t^"); put(dgy(i)(j)); new_line;
        end loop;
      end loop;
    end if;
    Sum(cfy,dgy);
    if vrblvl > 0 then
      put_line("After adding coefficients with same monomial :");
      for i in cfy'range loop
        put("all values of polynomial "); put(i,1); put_line(" :");
        for j in cfy(i)'range loop
          put(cfy(i)(j)); put(" t^"); put(dgy(i)(j)); new_line;
        end loop;
      end loop;
    end if;
    if vrblvl > 0 then
      put_line("Computing the positive minima :");
      sumerr := 0.0;
    end if;
    for i in cfy'range loop
      psm(i) := Positive_Minimum(cfy(i).all,dgy(i).all);       
      if vrblvl > 0 then
        put("p :"); put(psm(i));
        put(", lpw :"); put(lpw(i)); err := abs(psm(i) - lpw(i));
        put(", err :"); put(err,3); new_line;
        sumerr := sumerr + err;
      end if;
    end loop;
    if vrblvl > 0
     then put("error sum :"); put(sumerr,3); new_line;
    end if;
    for i in cfy'range loop
      cfp(i) := Coefficient(cfy(i).all,dgy(i).all,psm(i));
    end loop;
    for i in cfy'range loop
      Standard_Complex_Vectors.Clear(cfy(i));
      Standard_Floating_Vectors.Clear(dgy(i));
    end loop;
  end Leading_Powers_by_Evaluation;

  procedure Second_Order_Derivatives
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 ) is

    t : constant double_float := abs(Standard_Random_Numbers.Random);
    z,dhzt : Standard_Complex_Vectors.Vector(cff'range);
    idx : Standard_Integer_Vectors.Vector(z'range) := (z'range => 0);

  begin
    if vrblvl > 0
     then put("a random t :"); put(t); new_line;
    end if;
    for i in z'range loop
      z(i) := cff(i)(cff(i)'first);
      z(i) := z(i) + cff(i)(cff(i)'first)*(t**pwr(i)(pwr'first));
    end loop;
    if vrblvl > 0
     then put_line("first order evaluated at t :"); put_line(z);
    end if;
    for i in z'range loop
      idx(i) := 1;
      for j in i..z'last loop
        idx(j) := idx(j) + 1;
        put("idx :"); put(idx); new_line;
        dhzt := Laurent_Homotopy_Derivatives.Diff(hcf,hct,hdg,idx,z,t);
        put("derivatives at"); put(idx); put_line(" :"); put_line(dhzt);
        idx(j) := idx(j) - 1;
      end loop;
      idx(i) := 0;
    end loop;
  end Second_Order_Derivatives;

  procedure Run_Newton_Step
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 ) is

    nbr : constant integer32 := pwr(pwr'first)'last;
    ycf : Standard_Complex_Vectors.Vector(hcf'range);
    ydg : Standard_Floating_Vectors.Vector(hcf'range);
    lcf : Standard_Complex_Vectors.Vector(hcf'range);
    lpw : Standard_Floating_Vectors.Vector(hcf'range);
    cA : Standard_Complex_Matrices.Matrix(hcf'range,hcf'range);
    eA : Standard_Floating_Matrices.Matrix(hcf'range,hcf'range);
    posmin,err,sumerr : double_float;
    psm : Standard_Floating_Vectors.Vector(hcf'range);
    cfp : Standard_Complex_Vectors.Vector(hcf'range);

  begin
    if vrblvl > 0 then
      put("-> in Test_Newton_Puiseux.run_newton_step, deg : ");
      put(nbr,1);  put_line(" ...");
    end if;
    for i in hcf'range loop
      lcf(i) := cff(i)(cff(i)'first); -- get leading coefficients
      lpw(i) := pwr(i)(pwr(i)'first); -- get leading exponents
    end loop;
    sumerr := 0.0;
    for i in hct'range loop
      posmin := Positive_minimum(hct(i).all);
      put("p :"); put(posmin);
      put(", lpw :"); put(lpw(i)); err := abs(posmin - lpw(i));
      put(", err :"); put(err,3); new_line;
      sumerr := sumerr + err;
    end loop;
    put("error sum :"); put(sumerr,3); new_line;
    Evaluate_and_Differentiate(hcf,hct,hdg,lcf,ycf,ydg,cA,eA,1);
    if vrblvl > 0 then
      put_line("the function value :");
      for i in ycf'range loop
        put(ycf(i)); put(" t^"); put(ydg(i)); new_line;
      end loop;
      put_line("the Jacobian :");
      for i in cA'range(1) loop
        for j in cA'range(2) loop
          put("A("); put(i,1); put(","); put(j,1); put(") : ");
          put(cA(i,j)); put(" t^"); put(eA(i,j)); new_line;
        end loop;
      end loop;
    end if;
    Leading_Powers_by_Evaluation(hcf,hct,hdg,lcf,lpw,psm,cfp,vrblvl);
    for i in cfp'range loop -- Jacobian is diagonal for the test example
      cfp(i) := -cfp(i)/cA(i,i);
    end loop;
    for i in cfp'range loop
      put(cfp(i)); put(" t^"); put(psm(i)); new_line;
    end loop;
    sumerr := 0.0;
    for i in cff'range loop
      put(cff(i)(cff(i)'first+1));
      err := AbsVal(cfp(i) - cff(i)(cff(i)'first+1));
      put(", err :"); put(err,3); new_line;
      sumerr := sumerr + err;
    end loop;
    put("error sum :"); put(sumerr,3); new_line;
    Double_Ordered_Evaluations.First_Order_Evaluation
      (hcf,hct,hdg,cff,pwr,psm,cfp,vrblvl-1);
    put_line("smallest positive powers :");
    for i in psm'range loop
      put(i,1); put(" :"); put(psm(i)); new_line;
    end loop;
    put_line("first order terms :");
    for i in cfp'range loop
      put(cfp(i)); put(" t^"); put(psm(i)); new_line;
    end loop;
    if nbr > 1 then
      Double_Ordered_Evaluations.Second_Order_Evaluation
        (hcf,hct,hdg,cff,pwr,psm,cfp,vrblvl-1);
      put_line("smallest positive powers :");
      for i in psm'range loop
        put(i,1); put(" :"); put(psm(i)); new_line;
      end loop;
      put_line("Second order terms :");
      sumerr := 0.0;
      for i in cfp'range loop
        cfp(i) := -cfp(i)/cA(i,i);
        put(cfp(i)); put(" t^"); put(psm(i)); new_line;
        put(cff(i)(2)); put(" t^"); put(pwr(i)(2)); new_line;
        err := AbsVal(cfp(i) - cff(i)(2));
        sumerr := sumerr + err;
        put("error : "); put(err,3); put(" t^");
        err := abs(psm(i) - pwr(i)(2));
        put(err,3); new_line;
        sumerr := sumerr + err;
      end loop;
      put("sum of errors :"); put(sumerr,3); new_line;
     -- Second_Order_Derivatives(hcf,hct,hdg,cff,pwr,vrblvl-1);
      put_line("Computing third order evaluations ...");
      Double_Ordered_Evaluations.Third_Order_Evaluation
        (hcf,hct,hdg,cff,pwr,psm,cfp,vrblvl-1);
      put_line("smallest positive powers :");
      for i in psm'range loop
        put(i,1); put(" :"); put(psm(i)); new_line;
      end loop;
      put_line("second order terms :");
      for i in psm'range loop
        put(cff(i)(2)); put(" t^"); put(pwr(i)(2)); new_line;
      end loop;
    end if;
  end Run_Newton_Step;

  procedure Scale_Homotopy_Powers
              ( hct : in Standard_Floating_VecVecs.VecVec ) is

    powers : Standard_Floating_Vectors.Link_to_Vector;
    minpwr : double_float;

  begin
    for i in hct'range loop
      powers := hct(i);
      minpwr := powers(powers'first);
      for j in powers'first+1..powers'last loop
        if powers(j) < minpwr
         then minpwr := powers(j);
        end if;
      end loop;
      for j in powers'range loop
        powers(j) := powers(j) - minpwr;
      end loop;
    end loop;
  end Scale_Homotopy_Powers;

  procedure Define_Homotopy
              ( dim : in integer32;
                nbm,nbt : in Standard_Integer_Vectors.Vector;
                cff : out Standard_Complex_VecVecs.VecVec;
                pwr : out Standard_Floating_VecVecs.VecVec;
                hdg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : out Standard_Complex_VecVecs.VecVec;
                hct : out Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 ) is

    pdg : Standard_Integer_VecVecs.Array_of_VecVecs(1..dim);
    pcf : Standard_Complex_VecVecs.VecVec(1..dim);
    pct : Standard_Floating_VecVecs.VecVec(1..dim);

  begin
    Random_Laurent_Homotopy.Random_Laurent_System
      (dim,dim,-9,9,nbm,pdg,pcf,pct);
    if vrblvl > 0 then
      for i in 1..dim loop
        put("-> coefficients and degrees of polynomial ");
        put(i,1); put_line(" :");
        for j in 1..nbm(i) loop
          put(pcf(i)(j)); 
          put(" t^"); put(pct(i)(j));
          put("  "); put(pdg(i)(j)); new_line;
        end loop;
      end loop;
    end if;
    Random_Laurent_Homotopy.Random_Power_Series(dim,nbt,cff,pwr);
    if vrblvl > 0 then
      for i in 1..dim loop
        put("-> a random power series "); put(i,1); put_line(" : ");
        put(cff(i)(0)); new_line;
        for j in 1..nbt(i) loop
          put(cff(i)(j)); put("  t^"); put(pwr(i)(j)); new_line;
        end loop;
      end loop;
    end if;
    Random_Laurent_Homotopy.Random_Homotopy(pdg,pcf,pct,cff,pwr,hdg,hcf,hct);
    Scale_Homotopy_Powers(hct);
    if vrblvl > 0 then
      for i in 1..dim loop
        put("-> coefficients and degrees of homotopy ");
        put(i,1); put_line(" :");
        for j in hcf(i)'range loop
          put(hcf(i)(j)); 
          put(" t^"); put(hct(i)(j));
          put("  "); put(hdg(i)(j)); new_line;
        end loop;
      end loop;
    end if;
    Random_Laurent_Homotopy.Test_Random_Homotopy(hdg,hcf,hct,cff,pwr);
  end Define_Homotopy;

  procedure Test ( dim : in integer32 ) is

    nbm : Standard_Integer_Vectors.Vector(1..dim) := (1..dim => 0);
    nbt : Standard_Integer_Vectors.Vector(1..dim) := (1..dim => 0);
    cff : Standard_Complex_VecVecs.VecVec(1..dim);
    pwr : Standard_Floating_VecVecs.VecVec(1..dim);
    hdg : Standard_Integer_VecVecs.Array_of_VecVecs(1..dim);
    hcf : Standard_Complex_VecVecs.VecVec(1..dim);
    hct : Standard_Floating_VecVecs.VecVec(1..dim);

  begin
    put_line("Reading the number of monomials for every polynomial ...");
    for i in 1..dim loop
      put("  Give the number of monomials in polynomial "); put(i,1);
      put(" : "); get(nbm(i));
    end loop;
    put_line("Reading the number of terms in every series ...");
    for i in 1..dim loop
      put("  Give the number of terms in series "); put(i,1);
      put(" : "); get(nbt(i));
    end loop;
    Define_Homotopy(dim,nbm,nbt,cff,pwr,hdg,hcf,hct,1);
    Run_Newton_Step(hcf,hct,hdg,cff,pwr,4);
  end Test;

  procedure Main is

    dim : integer32 := 0;

  begin
    new_line;
    put("Give the number of variables : "); get(dim);
    Test(dim);
  end Main;

end Test_Newton_Puiseux;
