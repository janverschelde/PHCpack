with Ada.Text_IO;                       use Ada.Text_IO;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Numbers_Polar;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Random_Numbers;
with Standard_Random_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;
with Double_Leading_Evaluations;
with Double_Real_Powered_Series;
with Random_Laurent_Homotopy;

package body Test_Leading_Evaluations is

  procedure Test_Monomial_Derivative
              ( deg : in Standard_Integer_Vectors.Vector;
                pwr : in Standard_Floating_Vectors.Vector;
                cff : in Standard_Complex_Vectors.Vector;
                idx : in integer32; err : out double_float ) is

    lpr : constant double_float
        := Double_Leading_Evaluations.Leading_Power(deg,pwr,idx);
    lcf : constant Complex_Number
        := Double_Leading_Evaluations.Leading_Coefficient(deg,cff,idx);
   -- variables for the random point test
    t_rnd : constant double_float := abs(Standard_Random_Numbers.Random);
    t_pwr : constant double_float := t_rnd**lpr;
    t_val : constant Complex_Number := t_pwr*lcf;
    t_cff : Standard_Complex_Vectors.Vector(cff'range);
    t_tst,error : Complex_Number;

  begin
    put("-> testing derivative with respect to variable "); put(idx,1);
    put_line(" ...");
    put("leading power : "); put(lpr); new_line;
    put("leading coefficient : "); put(lcf); new_line;
    put("value at random point t = "); put(t_rnd); put_line(" :");
    put(t_val); new_line;
    for i in t_cff'range loop
      t_cff(i) := (t_rnd**pwr(i))*cff(i);
    end loop;
    t_tst := Double_Leading_Evaluations.Leading_Coefficient(deg,t_cff,idx);
    put_line("comparing to alternative evaluation :");
    put(t_tst); new_line;
    error := t_tst - t_val;
    err := AbsVal(error);
    put("error :"); put(err,3);
    if err < 1.0E-10
     then put_line(", okay");
     else put_line(" >= 1.0E-10, bug!");
    end if;
  end Test_Monomial_Derivative;

  procedure Show_Indices ( dim,idxsum : in integer32 ) is

    cnt : integer32 := 0;

    procedure Write ( idx : in Standard_Integer_Vectors.Vector;
                      continue : out boolean ) is
    begin
      cnt := cnt + 1;
      put(cnt,3); put(" :"); put(idx); new_line;
      continue := true;
    end Write;
 
    procedure Write_Indices is
      new Double_Leading_Evaluations.Enumerate_Indices(Write);

  begin
    Write_Indices(dim,idxsum);
  end Show_Indices;

  procedure Test_Indexed_Derivative
              ( deg : in Standard_Integer_Vectors.Vector;
                cff : in Standard_Complex_Vectors.Vector;
                err : out double_float ) is

    x,y : Complex_Number;
    idx : Standard_Integer_Vectors.Vector(deg'range) := (deg'range => 0);
    idxerr : double_float;

    use Double_Leading_Evaluations;

  begin
    err := 0.0;
    put_line("testing first derivatives ...");
    for i in deg'range loop
      x := Leading_Coefficient(deg,cff,i,0);
      idx(i) := 1;
      y := Indexed_Derivative(deg,cff,idx);
      idx(i) := 0;
      idxerr := AbsVal(x - y);
      err := err + idxerr;
      put("x"); put(i,1); put(" : "); put(x); new_line;
      put("y"); put(i,1); put(" : "); put(y);
      put(", err :"); put(idxerr,3); new_line;
    end loop;
    put_line("testing second pure derivatives ...");
    for i in deg'range loop
      x := Second_Derivative(deg,cff,i,0);
      idx(i) := 2;
      y := Indexed_Derivative(deg,cff,idx);
      idx(i) := 0;
      idxerr := AbsVal(x - y);
      err := err + idxerr;
      put("x"); put(i,1); put(" : "); put(x); new_line;
      put("y"); put(i,1); put(" : "); put(y);
      put(", err :"); put(idxerr,3); new_line;
    end loop;
    put_line("testing second mixed derivatives ...");
    for i in deg'range loop
      for j in i+1..deg'last loop
        x := Second_Mixed_Derivative(deg,cff,i,j,0);
        idx(i) := 1; idx(j) := 1;
        y := Indexed_Derivative(deg,cff,idx);
        idx(i) := 0; idx(j) := 0;
        idxerr := AbsVal(x - y);
        err := err + idxerr;
        put("x"); put(i,1); put(","); put(j,1); put(" : "); put(x); new_line;
        put("y"); put(i,1); put(","); put(j,1); put(" : "); put(y);
        put(", err :"); put(idxerr,3); new_line;
      end loop;
    end loop;
    put_line("testing third pure derivatives ...");
    for i in deg'range loop
      x := Third_Derivative(deg,cff,i,0);
      idx(i) := 3;
      y := Indexed_Derivative(deg,cff,idx);
      idx(i) := 0;
      idxerr := AbsVal(x - y);
      err := err + idxerr;
      put("x"); put(i,1); put(" : "); put(x); new_line;
      put("y"); put(i,1); put(" : "); put(y);
      put(", err :"); put(idxerr,3); new_line;
    end loop;
    put_line("testing third semi mixed derivatives ...");
    for i in deg'range loop
      for j in i+1..deg'last loop
        x := Third_Semi_Mixed_Derivative(deg,cff,i,j,0);
        idx(i) := 2; idx(j) := 1;
        y := Indexed_Derivative(deg,cff,idx);
        idx(i) := 0; idx(j) := 0;
        idxerr := AbsVal(x - y);
        err := err + idxerr;
        put("x"); put(i,1); put(","); put(j,1); put(" : "); put(x); new_line;
        put("y"); put(i,1); put(","); put(j,1); put(" : "); put(y);
        put(", err :"); put(idxerr,3); new_line;
      end loop;
    end loop;
    put_line("testing third fully mixed derivatives ...");
    for i in deg'range loop
      for j in i+1..deg'last loop
        for k in j+1..deg'last loop
          x := Third_Fully_Mixed_Derivative(deg,cff,i,j,k,0);
          idx(i) := 1; idx(j) := 1; idx(k) := 1;
          y := Indexed_Derivative(deg,cff,idx);
          idx(i) := 0; idx(j) := 0; idx(k) := 0;
          idxerr := AbsVal(x - y);
          err := err + idxerr;
          put("x"); put(i,1); put(","); put(j,1); put(","); put(k,1);
          put(" : "); put(x); new_line;
          put("y"); put(i,1); put(","); put(j,1); put(","); put(k,1);
          put(" : "); put(y);
          put(", err :"); put(idxerr,3); new_line;
        end loop;
      end loop;
    end loop;
  end Test_Indexed_Derivative;

  procedure Test_Indexed_Monomial_Derivatives ( dim : in integer32 ) is
               
    deg : constant Standard_Integer_Vectors.Vector(1..dim)
        := Random_Laurent_Homotopy.Random_Monomial_Support(dim,-9,9);
    cff : constant Standard_Complex_Vectors.Vector(1..dim)
        := Standard_Random_Vectors.Random_Vector(1,dim);
    err : double_float;

  begin
    Test_Indexed_Derivative(deg,cff,err);
    put("sum of errors :"); put(err,3); new_line;
  end Test_Indexed_Monomial_Derivatives;

  procedure Test_Indexed_Derivatives is

    dim : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    Test_Indexed_Monomial_Derivatives(dim);
    put_line("All indices with sum equal to 4 :");
    Show_Indices(dim,4);
  end Test_Indexed_Derivatives;

  procedure Test_Monomial ( dim : in integer32 ) is
               
    deg : constant Standard_Integer_Vectors.Vector(1..dim)
        := Random_Laurent_Homotopy.Random_Monomial_Support(dim,-9,9);
    pwr : constant Standard_Floating_Vectors.Vector(1..dim)
        := Double_Real_Powered_Series.Random_Leading_Powers(dim);
    cff : constant Standard_Complex_Vectors.Vector(1..dim)
        := Standard_Random_Vectors.Random_Vector(1,dim);
    lpr : constant double_float
        := Double_Leading_Evaluations.Leading_Power(deg,pwr);
    lcf : constant Complex_Number
        := Double_Leading_Evaluations.Leading_Coefficient(deg,cff);
   -- variables for the random point test
    t_rnd : constant double_float := abs(Standard_Random_Numbers.Random);
    t_pwr : constant double_float := t_rnd**lpr;
    t_val : constant Complex_Number := t_pwr*lcf;
    t_cff : Standard_Complex_Vectors.Vector(1..dim);
    t_tst,error : Complex_Number;
    err,diferr : double_float;

  begin
    put("the degrees :"); put(deg); new_line;
    put_line("leading powers of the series :"); put_line(pwr);
    put_line("leading coefficients of the series :"); put_line(cff);
    put("power value : "); put(lpr); new_line;
    put("coeff value : "); put(lcf); new_line;
    put(t_rnd); put("**"); put(lpr); put(" = ");
    put(t_pwr); new_line;
    put("evaluating at t ="); put(t_rnd); put_line(" :");
    put(t_val); new_line;
    put("|value| :");
    put(Standard_Complex_Numbers_Polar.radius(t_val)); new_line;
    put("  angle : ");
    put(Standard_Complex_Numbers_Polar.angle(t_val)); new_line;
    for i in t_cff'range loop
      t_cff(i) := (t_rnd**pwr(i))*cff(i);
    end loop;
    put_line("leading coefficients multiplied by the point :");
    put_line(t_cff);
    t_tst := Double_Leading_Evaluations.Leading_Coefficient(deg,t_cff);
    put_line("value at the random point :");
    put(t_tst); new_line;
    put("|value| :");
    put(Standard_Complex_Numbers_Polar.radius(t_tst)); new_line;
    put("  angle : ");
    put(Standard_Complex_Numbers_Polar.angle(t_tst)); new_line;
    error := t_tst - t_val;
    err := AbsVal(error);
    put("error :"); put(err,3);
    if err < 1.0E-10
     then put_line(", okay");
     else put_line(" >= 1.0E-10, bug!");
    end if;
    for i in 1..dim loop
      if deg(i) /= 0 then
        Test_Monomial_Derivative(deg,pwr,cff,i,diferr);
        err := err + diferr;
      end if;
    end loop;
    put("sum of all errors :"); put(err,3);
    if err < 1.0E-10
     then put_line(", okay");
     else put_line(" >= 1.0E-10, bug!");
    end if;
  end Test_Monomial;

  procedure Test_Polynomial ( nbr,dim : in integer32 ) is

    deg : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Random_Laurent_Homotopy.Random_Polynomial_Support(nbr,dim,-9,9);
    pwr : constant Standard_Floating_Vectors.Vector(1..dim)
        := Double_Real_Powered_Series.Random_Leading_Powers(dim);
    val : Standard_Floating_Vectors.Vector(1..nbr);
   -- lpr : double_float;
    diflpr,df2 : double_float;
    idx : integer32;

  begin
    put_line("The degrees : ");
    for i in deg'range loop
      put(deg(i).all); new_line;
    end loop;
    put_line("leading powers of the series :"); put_line(pwr);
   -- Double_Leading_Evaluations.Leading_Power(deg,pwr,lpr,idx,1);
    Double_Leading_Evaluations.Evaluate_Powers(deg,pwr,val,idx,1);
    put_line("evaluated powers :"); put_line(val);
    put("power value : "); put(val(val'first));
    put(" at index "); put(idx,1); new_line;
    for k in 1..dim loop
      put("-> derivative "); put(k,1); 
      put(", at index "); put(idx,1);
      diflpr := Double_Leading_Evaluations.Leading_Power(deg(idx).all,pwr,k);
      put(" : "); put(diflpr); new_line;
      for j in 1..nbr loop
        if j /= idx then
          df2 := Double_Leading_Evaluations.Leading_Power(deg(j).all,pwr,k);
          put("  at monomial "); put(j,1);
          put(" : "); put(df2);
          if df2 < diflpr
           then put_line(", larger!");
           else put_line(", okay.");
          end if;
        end if;
      end loop;
    end loop;
  end Test_Polynomial;

  procedure Test_System
              ( nbp,dim : in integer32;
                nbm : in Standard_Integer_Vectors.Vector ) is

    pdg : Standard_Integer_VecVecs.Array_of_VecVecs(1..nbp);
    pcf : Standard_Complex_VecVecs.VecVec(1..nbp);
    pct : Standard_Floating_VecVecs.VecVec(1..nbp);
    pwr : constant Standard_Floating_Vectors.Vector(1..dim)
        := Double_Real_Powered_Series.Random_Leading_Powers(dim);
    cff : constant Standard_Complex_Vectors.Vector(1..dim)
        := Standard_Random_Vectors.Random_Vector(1,dim);

  begin
    Random_Laurent_Homotopy.Random_Laurent_System
      (nbp,dim,-9,9,nbm,pdg,pcf,pct);
    for i in 1..nbp loop
      put("-> coefficients and degrees of polynomial ");
      put(i,1); put_line(" :");
      for j in 1..nbm(i) loop
        put(pcf(i)(j)); 
        put(" t^"); put(pct(i)(j));
        put("  "); put(pdg(i)(j)); new_line;
      end loop;
    end loop;
    put_line("Leading coefficients and powers of the series :");
    for i in 1..dim loop
      put(cff(i)); put("  t^"); put(pwr(i)); new_line;
    end loop;
    for i in 1..nbp loop
      declare
        val : Standard_Floating_Vectors.Vector(1..nbm(i));
        ycf : Standard_Complex_Vectors.Vector(1..nbm(i));
        idx : integer32;
      begin
        put("Evaluating polynomial "); put(i,1); put_line(" ...");
        Double_Leading_Evaluations.Evaluate_Polynomial
          (pcf(i).all,pct(i).all,pdg(i).all,cff,pwr,ycf,val,idx,1);
        put_line("evaluated powers :"); put_line(val);
        put("power value : "); put(val(val'first));
        put(" at index "); put(idx,1); new_line;
      end;
    end loop;
    put_line("Evaluating the system ...");
    declare
      ycf : Standard_Complex_VecVecs.VecVec(1..nbp);
      ydg : Standard_Floating_VecVecs.VecVec(1..nbp);
    begin
      for i in 1..nbp loop
        ycf(i) := new Standard_Complex_Vectors.Vector(1..nbm(i));
        ydg(i) := new Standard_Floating_Vectors.Vector(1..nbm(i));
      end loop;
      Double_Leading_Evaluations.Evaluate_System
        (pcf,pct,pdg,cff,pwr,ycf,ydg,1);
      for i in 1..nbp loop
        put("evaluated powers at "); put(i,1); put_line(" :");
        put_line(ydg(i).all);
      end loop;
    end;   
  end Test_System;

  procedure Test_Homotopy
              ( dim : in integer32;
                nbm,nbt : in Standard_Integer_Vectors.Vector ) is

    pdg : Standard_Integer_VecVecs.Array_of_VecVecs(1..dim);
    pcf : Standard_Complex_VecVecs.VecVec(1..dim);
    pct : Standard_Floating_VecVecs.VecVec(1..dim);
    cff : Standard_Complex_VecVecs.VecVec(1..dim);
    pwr : Standard_Floating_VecVecs.VecVec(1..dim);
    hdg : Standard_Integer_VecVecs.Array_of_VecVecs(1..dim);
    hcf : Standard_Complex_VecVecs.VecVec(1..dim);
    hct : Standard_Floating_VecVecs.VecVec(1..dim);

  begin
    Random_Laurent_Homotopy.Random_Laurent_System
      (dim,dim,-9,9,nbm,pdg,pcf,pct);
    for i in 1..dim loop
      put("-> coefficients and degrees of polynomial ");
      put(i,1); put_line(" :");
      for j in 1..nbm(i) loop
        put(pcf(i)(j)); 
        put(" t^"); put(pct(i)(j));
        put("  "); put(pdg(i)(j)); new_line;
      end loop;
    end loop;
    Double_Real_Powered_Series.Random_Power_Series(dim,nbt,cff,pwr);
    for i in 1..dim loop
      put("-> a random power series "); put(i,1); put_line(" : ");
      put(cff(i)(0)); new_line;
      for j in 1..nbt(i) loop
        put(cff(i)(j)); put("  t^"); put(pwr(i)(j)); new_line;
      end loop;
    end loop;
    Random_Laurent_Homotopy.Random_Homotopy(pdg,pcf,pct,cff,pwr,hdg,hcf,hct);
    for i in 1..dim loop
      put("-> coefficients and degrees of homotopy ");
      put(i,1); put_line(" :");
      for j in hcf(i)'range loop
        put(hcf(i)(j)); 
        put(" t^"); put(hct(i)(j));
        put("  "); put(hdg(i)(j)); new_line;
      end loop;
    end loop;
    Random_Laurent_Homotopy.Test_Random_Homotopy(hdg,hcf,hct,cff,pwr);
  end Test_Homotopy;

  procedure Main is

    nbr,dim,npol : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the number of variables : "); get(dim);
    put("Give the number of polynomials : "); get(npol);
    if npol = 1 and dim > 1 then
      put("Give the number of monomials : "); get(nbr);
      if nbr = 1
       then Test_Monomial(dim);
       else Test_Polynomial(nbr,dim);
      end if;
    else
      put_line("Reading the number of monomials for every polynomial ...");
      declare
        nbm : Standard_Integer_Vectors.Vector(1..npol) := (1..npol => 0);
        nbt : Standard_Integer_Vectors.Vector(1..dim) := (1..dim => 0);
      begin
        for i in 1..npol loop
          put("  Give the number of monomials in polynomial "); put(i,1);
          put(" : "); get(nbm(i));
        end loop;
        if npol /= dim then
          Test_System(npol,dim,nbm);
        else
          put("Test homotopy ? (y/n) ");
          Ask_Yes_or_No(ans);
          if ans = 'n' then
            Test_System(npol,dim,nbm);
          else
            for i in 1..dim loop
              put("  Give the number of terms in series "); put(i,1);
              put(" : "); get(nbt(i));
            end loop;
            Test_Homotopy(dim,nbm,nbt);
          end if;
        end if;
      end;
    end if;
  end Main;

end Test_Leading_Evaluations;
