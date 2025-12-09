with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Double_Leading_Evaluations;
with Double_Puiseux_Operations;
with Random_Laurent_Homotopy;

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
      end loop;
    end loop;
  end Evaluate_and_Differentiate;

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

  procedure Run_Newton_Step
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := hcf'last;
    ycf : Standard_Complex_Vectors.Vector(hcf'range);
    ydg : Standard_Floating_Vectors.Vector(hcf'range);
    lcf : Standard_Complex_Vectors.Vector(hcf'range);
    lpw : Standard_Floating_Vectors.Vector(hcf'range);
    cA : Standard_Complex_Matrices.Matrix(hcf'range,hcf'range);
    eA : Standard_Floating_Matrices.Matrix(hcf'range,hcf'range);
    tol : constant double_float := 1.0e-12;
    idx1,idx2 : Standard_Integer_Vectors.Vector(hcf'range);
    work : Standard_Integer_VecVecs.VecVec(0..dim);
    sol : Standard_Floating_Vectors.Vector(hcf'range);
    fail : boolean;
    posmin,err,sumerr : double_float;

  begin
    if vrblvl > 0 then
      put_line("-> in Test_Newton_Puiseux.run_newton_step ...");
    end if;
    for i in hcf'range loop
      lcf(i) := cff(i)(cff(i)'first); -- get leading coefficients
      lpw(i) := pwr(i)(pwr(i)'first); -- get leading exponents
    end loop;
    sumerr := 0.0;
    for i in hct'range loop
      posmin := Positive_minimum(hct(i).all);
      put("p : "); put(posmin);
      put(", lpw : "); put(lpw(i)); err := abs(posmin - lpw(i));
      put(", err : "); put(err,3); new_line;
      sumerr := sumerr + err;
    end loop;
    put("error sum : "); put(sumerr,3); new_line;
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
    put_line("-> computing the tropical Cramer vector ...");
    for i in work'range loop
      work(i) := new Standard_Integer_Vectors.Vector'(1..dim => 0);
    end loop;
    Double_Puiseux_Operations.Leading_Powers
      (dim,tol,eA,ydg,work,sol,idx1,idx2,fail,2);
    put_line("the solution :"); put_line(sol);
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
    Run_Newton_Step(hcf,hct,hdg,cff,pwr,1);
  end Test;

  procedure Main is

    dim : integer32 := 0;

  begin
    new_line;
    put("Give the number of variables : "); get(dim);
    Test(dim);
  end Main;

end Test_Newton_Puiseux;
