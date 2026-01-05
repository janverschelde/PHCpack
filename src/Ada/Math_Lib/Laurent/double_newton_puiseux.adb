with Ada.Text_IO;                       use Ada.Text_IO;
with Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Double_Leading_Evaluations;
with Double_Ordered_Evaluations;
with Double_Real_Powered_Series;
with Laurent_Homotopy_Derivatives;

package body Double_Newton_Puiseux is

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
      put_line("-> in Double_Newton_Puiseux.evaluate_and_differentiate ...");
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
      if vrblvl > 0
       then put(ycf(i)); put(" t^"); put(ydg(i)); new_line;
      end if;
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
        if vrblvl > 0
         then put(ycf(i)); put(" t^"); put(ydg(i)); new_line;
        end if;
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
      put_line("-> in Double_Newton_Puiseux.evaluate_all_monomials ...");
    end if;
    for i in 1..dim loop -- run over all polynomials
      for j in hcf(i)'range loop -- run over all monomials
        if vrblvl > 0 then
          put("at monomial "); put(j,1);
          put(" of polynomial "); put(i,1); put_line(" ...");
        end if;
        ydg(i)(j) := hct(i)(j);
        ycf(i)(j) := hcf(i)(j)*Leading_Coefficient(hdg(i)(j).all,zt0);
        if vrblvl > 0
         then put(ycf(i)(j)); put(" t^"); put(ydg(i)(j)); new_line;
        end if;
      end loop;
    end loop;
  end Evaluate_All_Monomials;

  procedure Leading_Powers_by_Evaluation
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                lcf : in Standard_Complex_Vectors.Vector;
                psm : out Standard_Floating_Vectors.Vector; 
                cfp : out Standard_Complex_Vectors.Vector; 
                vrblvl : in integer32 := 0 ) is

    cfy : Standard_Complex_VecVecs.VecVec(hcf'range);
    dgy : Standard_Floating_VecVecs.VecVec(hct'range);

  begin
    if vrblvl > 0 then
      put("-> in Double_Newton_Puiseux.");
      put_line("leading_powers_by_evaluations ...");
    end if;
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
    Double_Real_Powered_Series.Normalize(cfy,dgy);
    if vrblvl > 0 then
      put_line("After adding coefficients with same monomial :");
      for i in cfy'range loop
        put("all values of polynomial "); put(i,1); put_line(" :");
        for j in cfy(i)'range loop
          put(cfy(i)(j)); put(" t^"); put(dgy(i)(j)); new_line;
        end loop;
      end loop;
    end if;
    psm := Double_Real_Powered_Series.Positive_Minima(cfy,dgy);       
    cfp := Double_Real_Powered_Series.Coefficients(cfy,dgy,psm);
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
    if vrblvl > 0 then
      put_line("-> in Double_Newton_Puiseux.second_order_derivatives ...");
    end if;
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

  procedure Diagonal_Leading_Terms
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0 : in Standard_Complex_Vectors.Vector;
                cA : out Standard_Complex_Matrices.Matrix;
                eA : out Standard_Floating_Matrices.Matrix;
                cf1 : out Standard_Complex_Vectors.Vector;
                pw1 : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    ycf : Standard_Complex_Vectors.Vector(hcf'range);
    ydg : Standard_Floating_Vectors.Vector(hcf'range);

  begin
    if vrblvl > 0
     then put_line("-> in Double_Newton_Puiseux.diagonal_leading_terms ...");
    end if;
    Evaluate_and_Differentiate(hcf,hct,hdg,cf0,ycf,ydg,cA,eA,vrblvl-1);
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
    Leading_Powers_by_Evaluation(hcf,hct,hdg,cf0,pw1,cf1,vrblvl-1);
    for i in cf1'range loop -- Jacobian is diagonal for the test example
      cf1(i) := -cf1(i)/cA(i,i);
    end loop;
    if vrblvl > 0 then
      for i in cf1'range loop
        put(cf1(i)); put(" t^"); put(pw1(i)); new_line;
      end loop;
    end if;
  end Diagonal_Leading_Terms;

  procedure Diagonal_Second_Terms
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0,cf1 : in Standard_Complex_Vectors.Vector;
                pw1 : in Standard_Floating_Vectors.Vector;
                cA : in Standard_Complex_Matrices.Matrix;
                cf2 : out Standard_Complex_Vectors.Vector;
                pw2 : out Standard_Floating_Vectors.Vector;
                tol : in double_float := 1.0E-12;
                vrblvl : in integer32 := 0 ) is

    cf2a : Standard_Complex_Vectors.Vector(cf1'range);
    pw2a : Standard_Floating_Vectors.Vector(pw1'range);
    dif,sumdif : double_float;

  begin
    if vrblvl > 0
     then put_line("-> in Double_Newton_Puiseux.diagonal_second_terms ...");
    end if;
    Double_Ordered_Evaluations.First_Derivative_First_Order
      (hcf,hct,hdg,cf0,cf1,pw1,cf2a,pw2a,vrblvl-1);
    for i in cf1'range loop -- Jacobian is diagonal for the test example
      cf2a(i) := -cf2a(i)/cA(i,i);
    end loop;
    Double_Ordered_Evaluations.Second_Derivative_First_Order
      (hcf,hct,hdg,cf0,cf1,pw1,cf2,pw2,vrblvl-1);
    for i in cf1'range loop -- Jacobian is diagonal for the test example
      cf2(i) := -cf2(i)/cA(i,i);
    end loop;
    if vrblvl > 0
     then put_line("first and second derivative evaluations :");
    end if;
    sumdif := 0.0;
    for i in cf1'range loop
      if vrblvl > 0 then
        put(cf2a(i)); put(" t^"); put(pw2a(i)); new_line;
        put(cf2(i)); put(" t^"); put(pw2(i)); new_line;
      end if;
      dif := AbsVal(cf2a(i) - cf2(i)) + abs(pw2a(i)-pw2(i));
      if vrblvl > 0
       then put("difference :"); put(dif,3); new_line;
      end if;
      sumdif := sumdif + dif;
    end loop;
    if vrblvl > 0
     then put("sum of differences :"); put(sumdif,3); new_line;
    end if;
    if sumdif > tol then
      cf2a := cf2; pw2a := pw2;
      Double_Ordered_Evaluations.Third_Derivative_First_Order
        (hcf,hct,hdg,cf0,cf1,pw1,cf2,pw2,vrblvl-1);
      for i in cf1'range loop -- Jacobian is diagonal for the test example
        cf2(i) := -cf2(i)/cA(i,i);
      end loop;
      if vrblvl > 0
       then put_line("second and third derivative evaluations :");
      end if;
      sumdif := 0.0;
      for i in cf1'range loop
        if vrblvl > 0 then
          put(cf2a(i)); put(" t^"); put(pw2a(i)); new_line;
          put(cf2(i)); put(" t^"); put(pw2(i)); new_line;
        end if;
        dif := AbsVal(cf2a(i) - cf2(i)) + abs(pw2a(i)-pw2(i));
        if vrblvl > 0
         then put("difference :"); put(dif,3); new_line;
        end if;
        sumdif := sumdif + dif;
      end loop;
      if vrblvl > 0
       then put("sum of differences :"); put(sumdif,3); new_line;
      end if;
    end if;
    for difsum in 4..8 loop
      exit when (sumdif < tol);
      cf2a := cf2; pw2a := pw2;
      if vrblvl > 0 then
        put("computing derivative "); put(integer32(difsum),1);
        put_line(" first order evaluations ...");
      end if;
      Double_Ordered_Evaluations.First_Order_Evaluation
        (hcf,hct,hdg,cf0,cf1,pw1,integer32(difsum),cf2,pw2,vrblvl-1);
      for i in cf1'range loop -- Jacobian is diagonal for the test example
        cf2(i) := -cf2(i)/cA(i,i);
      end loop;
      if vrblvl > 0 then
        put("derivative "); put(integer32(difsum)-1,1);
        put(" and "); put(integer32(difsum),1); put_line(" evaluations :");
      end if;
      sumdif := 0.0;
      for i in cf1'range loop
        if vrblvl > 0 then
          put(cf2a(i)); put(" t^"); put(pw2a(i)); new_line;
          put(cf2(i)); put(" t^"); put(pw2(i)); new_line;
        end if;
        dif := AbsVal(cf2a(i) - cf2(i)) + abs(pw2a(i)-pw2(i));
        if vrblvl > 0
         then put("difference :"); put(dif,3); new_line;
        end if;
        sumdif := sumdif + dif;
      end loop;
      if vrblvl > 0
       then put("sum of differences :"); put(sumdif,3); new_line;
      end if;
    end loop;
  end Diagonal_Second_Terms;

  procedure Diagonal_Third_Terms
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0,cf1,cf2 : in Standard_Complex_Vectors.Vector;
                pw1,pw2 : in Standard_Floating_Vectors.Vector;
                cA : in Standard_Complex_Matrices.Matrix;
                cf3 : out Standard_Complex_Vectors.Vector;
                pw3 : out Standard_Floating_Vectors.Vector;
                tol : in double_float := 1.0E-12;
                vrblvl : in integer32 := 0 ) is

    cf3a : Standard_Complex_Vectors.Vector(cf1'range);
    pw3a : Standard_Floating_Vectors.Vector(pw1'range);
    dif,sumdif : double_float;

  begin
    if vrblvl > 0
     then put_line("-> in Double_Newton_Puiseux.diagonal_third_terms ...");
    end if;
    Double_Ordered_Evaluations.First_Derivative_Second_Order
      (hcf,hct,hdg,cf0,cf1,cf2,pw1,pw2,cf3a,pw3a,vrblvl-1);
    for i in cf3a'range loop -- Jacobian is diagonal for the test example
      cf3a(i) := -cf3a(i)/cA(i,i);
    end loop;
    Double_Ordered_Evaluations.Second_Derivative_Second_Order
      (hcf,hct,hdg,cf0,cf1,cf2,pw1,pw2,cf3,pw3,vrblvl-1);
    for i in cf3'range loop -- Jacobian is diagonal for the test example
      cf3(i) := -cf3(i)/cA(i,i);
    end loop;
    if vrblvl > 0
     then put_line("first and second derivative evaluations :");
    end if;
    sumdif := 0.0;
    for i in cf1'range loop
      if vrblvl > 0 then
        put(cf3a(i)); put(" t^"); put(pw3a(i)); new_line;
        put(cf3(i)); put(" t^"); put(pw3(i)); new_line;
      end if;
      dif := AbsVal(cf3a(i) - cf3(i)) + abs(pw3a(i)-pw3(i));
      if vrblvl > 0
       then put("difference :"); put(dif,3); new_line;
      end if;
      sumdif := sumdif + dif;
    end loop;
    if vrblvl > 0 then
      put("sum of differences :"); put(sumdif,3); new_line;
      if sumdif > tol
       then put_line("higher derivatives are needed ...");
      end if;
    end if;
    if sumdif > tol then
      cf3a := cf3; pw3a := pw3;
      Double_Ordered_Evaluations.Third_Derivative_Second_Order
        (hcf,hct,hdg,cf0,cf1,cf2,pw1,pw2,cf3,pw3,vrblvl-1);
      for i in cf3'range loop -- Jacobian is diagonal for the test example
        cf3(i) := -cf3(i)/cA(i,i);
      end loop;
      if vrblvl > 0
       then put_line("second and third derivative evaluations :");
      end if;
      sumdif := 0.0;
      for i in cf1'range loop
        if vrblvl > 0 then
          put(cf3a(i)); put(" t^"); put(pw3a(i)); new_line;
          put(cf3(i)); put(" t^"); put(pw3(i)); new_line;
        end if;
        dif := AbsVal(cf3a(i) - cf3(i)) + abs(pw3a(i)-pw3(i));
        if vrblvl > 0
         then put("difference :"); put(dif,3); new_line;
        end if;
        sumdif := sumdif + dif;
      end loop;
      if vrblvl > 0 then
        put("sum of differences :"); put(sumdif,3); new_line;
        if sumdif > tol
         then put_line("higher derivatives are needed ...");
        end if;
      end if;
    end if;
    for difsum in 4..5 loop
      exit when (sumdif < tol);
      cf3a := cf3; pw3a := pw3;
      if vrblvl > 0 then
        put("computing derivative "); put(integer32(difsum),1);
        put_line(" second order evaluations ...");
      end if;
      Double_Ordered_Evaluations.Second_Order_Evaluation
        (hcf,hct,hdg,cf0,cf1,cf2,pw1,pw2,integer32(difsum),cf3,pw3,vrblvl-1);
      for i in cf1'range loop -- Jacobian is diagonal for the test example
        cf3(i) := -cf3(i)/cA(i,i);
      end loop;
      if vrblvl > 0 then
        put("derivative "); put(integer32(difsum)-1,1);
        put(" and "); put(integer32(difsum),1); put_line(" evaluations :");
      end if;
      sumdif := 0.0;
      for i in cf1'range loop
        if vrblvl > 0 then
          put(cf3a(i)); put(" t^"); put(pw3a(i)); new_line;
          put(cf3(i)); put(" t^"); put(pw3(i)); new_line;
        end if;
        dif := AbsVal(cf3a(i) - cf3(i)) + abs(pw3a(i)-pw3(i));
        if vrblvl > 0
         then put("difference :"); put(dif,3); new_line;
        end if;
        sumdif := sumdif + dif;
      end loop;
      if vrblvl > 0
       then put("sum of differences :"); put(sumdif,3); new_line;
      end if;
    end loop;
  end Diagonal_Third_Terms;

  procedure Diagonal_Fourth_Terms
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0,cf1,cf2,cf3 : in Standard_Complex_Vectors.Vector;
                pw1,pw2,pw3 : in Standard_Floating_Vectors.Vector;
                cA : in Standard_Complex_Matrices.Matrix;
                cf4 : out Standard_Complex_Vectors.Vector;
                pw4 : out Standard_Floating_Vectors.Vector;
                tol : in double_float := 1.0E-12;
                vrblvl : in integer32 := 0 ) is

    cf4a : Standard_Complex_Vectors.Vector(cf1'range);
    pw4a : Standard_Floating_Vectors.Vector(pw1'range);
    dif,sumdif : double_float;

  begin
    if vrblvl > 0
     then put_line("-> in Double_Newton_Puiseux.diagonal_fourth_terms ...");
    end if;
    Double_Ordered_Evaluations.First_Derivative_Third_Order
      (hcf,hct,hdg,cf0,cf1,cf2,cf3,pw1,pw2,pw3,cf4a,pw4a,vrblvl-1);
    for i in cf4a'range loop -- Jacobian is diagonal for the test example
      cf4a(i) := -cf4a(i)/cA(i,i);
    end loop;
    Double_Ordered_Evaluations.Second_Derivative_Third_Order
      (hcf,hct,hdg,cf0,cf1,cf2,cf3,pw1,pw2,pw3,cf4,pw4,vrblvl-1);
    for i in cf3'range loop -- Jacobian is diagonal for the test example
      cf4(i) := -cf4(i)/cA(i,i);
    end loop;
    if vrblvl > 0
     then put_line("first and second derivative evaluations :");
    end if;
    sumdif := 0.0;
    for i in cf1'range loop
      if vrblvl > 0 then
        put(cf4a(i)); put(" t^"); put(pw4a(i)); new_line;
        put(cf4(i)); put(" t^"); put(pw4(i)); new_line;
      end if;
      dif := AbsVal(cf4a(i) - cf4(i)) + abs(pw4a(i) - pw4(i));
      if vrblvl > 0
       then put("difference :"); put(dif,3); new_line;
      end if;
      sumdif := sumdif + dif;
    end loop;
    if vrblvl > 0 then
      put("sum of differences :"); put(sumdif,3); new_line;
      if sumdif > tol
       then put_line("higher derivatives are needed ...");
      end if;
    end if;
    if sumdif > tol then
      cf4a := cf4; pw4a := pw4;
      Double_Ordered_Evaluations.Third_Derivative_Third_Order
        (hcf,hct,hdg,cf0,cf1,cf2,cf3,pw1,pw2,pw3,cf4,pw4,vrblvl-1);
      for i in cf4'range loop -- Jacobian is diagonal for the test example
        cf4(i) := -cf4(i)/cA(i,i);
      end loop;
      if vrblvl > 0
       then put_line("second and third derivative evaluations :");
      end if;
      sumdif := 0.0;
      for i in cf1'range loop
        if vrblvl > 0 then
          put(cf4a(i)); put(" t^"); put(pw4a(i)); new_line;
          put(cf4(i)); put(" t^"); put(pw4(i)); new_line;
        end if;
        dif := AbsVal(cf4a(i) - cf4(i)) + abs(pw4a(i) - pw4(i));
        if vrblvl > 0
         then put("difference :"); put(dif,3); new_line;
        end if;
        sumdif := sumdif + dif;
      end loop;
      if vrblvl > 0 then
        put("sum of differences :"); put(sumdif,3); new_line;
        if sumdif > tol
         then put_line("higher derivatives are needed ...");
        end if;
      end if;
    end if;
    for difsum in 4..5 loop
      exit when (sumdif < tol);
      cf4a := cf4; pw4a := pw4;
      if vrblvl > 0 then
        put("computing derivative "); put(integer32(difsum),1);
        put_line(" third order evaluations ...");
      end if;
      Double_Ordered_Evaluations.Third_Order_Evaluation
        (hcf,hct,hdg,cf0,cf1,cf2,cf3,pw1,pw2,pw3,integer32(difsum),
         cf4,pw4,vrblvl-1);
      for i in cf1'range loop -- Jacobian is diagonal for the test example
        cf4(i) := -cf4(i)/cA(i,i);
      end loop;
      if vrblvl > 0 then
        put("derivative "); put(integer32(difsum)-1,1);
        put(" and "); put(integer32(difsum),1); put_line(" evaluations :");
      end if;
      sumdif := 0.0;
      for i in cf1'range loop
        if vrblvl > 0 then
          put(cf4a(i)); put(" t^"); put(pw4a(i)); new_line;
          put(cf4(i)); put(" t^"); put(pw4(i)); new_line;
        end if;
        dif := AbsVal(cf4a(i) - cf4(i)) + abs(pw4a(i) - pw4(i));
        if vrblvl > 0
         then put("difference :"); put(dif,3); new_line;
        end if;
        sumdif := sumdif + dif;
      end loop;
      if vrblvl > 0
       then put("sum of differences :"); put(sumdif,3); new_line;
      end if;
    end loop;
  end Diagonal_Fourth_Terms;

  function Error_Sum ( idx : in integer32;
                       cff : in Standard_Complex_VecVecs.VecVec;
                       pwr : in Standard_Floating_VecVecs.VecVec;
                       cfp : in Standard_Complex_Vectors.Vector;
                       psm : in Standard_Floating_Vectors.Vector;
                       vrblvl : in integer32 := 0 )
                     return double_float is

  -- DESCRIPTION :
  --   Given the generated power series in (cff, pwr) and the
  --   computed new coefficients and powers in (cfp, pwr),
  --   returns the sum of the errors.
  --   If vrblvl > 0, then the series and errors are written.
                     
    res : double_float := 0.0;
    err : double_float;

  begin
    for i in cff'range loop
      if vrblvl > 0 then
        put("x"); put(i,1); put(" :");
        put(cff(i)(idx)); put(" t^"); put(pwr(i)(idx)); new_line;
        put("y"); put(i,1); put(" :");
        put(cfp(i)); put(" t^"); put(psm(i)); new_line;
      end if;
      err := AbsVal(cfp(i) - cff(i)(idx)); res := res + err;
      if vrblvl > 0 then
        put("error :"); put(err,3);
        put(" t^");
      end if;
      err := abs(psm(i) - pwr(i)(idx)); res := res + err;
      if vrblvl > 0 then
        put(err,3); new_line;
      end if;
    end loop;
    return res;
  end Error_Sum;

  procedure Diagonal_Newton_Steps
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cf0 : in Standard_Complex_Vectors.Vector;
                nbr : in integer32;
                cf1,cf2,cf3,cf4 : out Standard_Complex_Vectors.Vector; 
                pw1,pw2,pw3,pw4 : out Standard_Floating_Vectors.Vector;
                tol : in double_float := 1.0E-12;
                vrblvl : in integer32 := 0 ) is

    cA : Standard_Complex_Matrices.Matrix(hcf'range,hcf'range);
    eA : Standard_Floating_Matrices.Matrix(hcf'range,hcf'range);

  begin
    if vrblvl > 0
     then put_line("-> in Double_Newton_Puiseux.four_newton_steps ...");
    end if;
    Diagonal_Leading_Terms(hcf,hct,hdg,cf0,cA,eA,cf1,pw1,vrblvl-1);
    if nbr > 1 then
      Diagonal_Second_Terms(hcf,hct,hdg,cf0,cf1,pw1,cA,cf2,pw2,tol,vrblvl-1);
      if nbr > 2 then
        Diagonal_Third_Terms
          (hcf,hct,hdg,cf0,cf1,cf2,pw1,pw2,cA,cf3,pw3,tol,vrblvl-1);
        if nbr > 3 then
          Diagonal_Fourth_Terms
            (hcf,hct,hdg,cf0,cf1,cf2,cf3,pw1,pw2,pw3,cA,cf4,pw4,tol,vrblvl-1);
        end if;
      end if;
    end if;
  end Diagonal_Newton_Steps;

  procedure Run_Newton_Step
              ( hcf : in Standard_Complex_VecVecs.VecVec;
                hct : in Standard_Floating_VecVecs.VecVec;
                hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                tol : in double_float := 1.0E-12;
                vrblvl : in integer32 := 0 ) is

    nbr : constant integer32 := pwr(pwr'first)'last;
    lcf : Standard_Complex_Vectors.Vector(hcf'range);
    cA : Standard_Complex_Matrices.Matrix(hcf'range,hcf'range);
    eA : Standard_Floating_Matrices.Matrix(hcf'range,hcf'range);
    sumerr : double_float;
    psm,pw2,pw3,pw4 : Standard_Floating_Vectors.Vector(hcf'range);
    cfp,cf2,cf3,cf4 : Standard_Complex_Vectors.Vector(hcf'range);
    ans : character;

  begin
    if vrblvl > 0 then
      put("-> in Double_Newton_Puiseux.run_newton_step, deg : ");
      put(nbr,1);  put_line(" ...");
    end if;
    for i in hcf'range loop
      lcf(i) := cff(i)(cff(i)'first); -- get leading coefficients
    end loop;
    Diagonal_Leading_Terms(hcf,hct,hdg,lcf,cA,eA,cfp,psm,1);
    put_line("Computing error sum ...");
    sumerr := Error_Sum(1,cff,pwr,cfp,psm,1);
    put("error sum :"); put(sumerr,3); new_line;
    if nbr > 1 then
      put("Continue ? (y/n) "); Communications_with_User.Ask_Yes_or_No(ans);
      if ans /= 'y'
       then return;
      end if;
      Diagonal_Second_Terms(hcf,hct,hdg,lcf,cfp,psm,cA,cf2,pw2,tol,1);
      put_line("Computing error sum ...");
      sumerr := Error_Sum(2,cff,pwr,cf2,pw2,1);
      put("error sum :"); put(sumerr,3); new_line;
    end if;
    if nbr > 2 then
      put("Continue ? (y/n) "); Communications_with_User.Ask_Yes_or_No(ans);
      if ans /= 'y'
       then return;
      end if;
      Diagonal_Third_Terms(hcf,hct,hdg,lcf,cfp,cf2,psm,pw2,cA,cf3,pw3,tol,1);
      put_line("Computing error sum ...");
      sumerr := Error_Sum(3,cff,pwr,cf3,pw3,1);
      put("error sum :"); put(sumerr,3); new_line;
   end if;
   if nbr > 3 then
      put("Continue ? (y/n) "); Communications_with_User.Ask_Yes_or_No(ans);
      if ans /= 'y'
       then return;
      end if;
      Diagonal_Fourth_Terms
        (hcf,hct,hdg,lcf,cfp,cf2,cf3,psm,pw2,pw3,cA,cf4,pw4,tol,1);
      put_line("Computing error sum ...");
      sumerr := Error_Sum(4,cff,pwr,cf4,pw4,1);
      put("error sum :"); put(sumerr,3); new_line;
      put_line("first order terms :");
      for i in psm'range loop
        put(cff(i)(1)); put(" t^"); put(pwr(i)(1)); new_line;
      end loop;
      if nbr > 1 then
        put_line("second order terms :");
        for i in psm'range loop
          put(cff(i)(2)); put(" t^"); put(pwr(i)(2)); new_line;
        end loop;
        if nbr > 2 then
          put_line("third order terms :");
          for i in psm'range loop
            put(cff(i)(3)); put(" t^"); put(pwr(i)(3)); new_line;
          end loop;
          if nbr > 3 then
            put_line("fourth order terms :");
            for i in psm'range loop
              put(cff(i)(4)); put(" t^"); put(pwr(i)(4)); new_line;
            end loop;
          end if;
        end if;
      end if;
    end if;
  end Run_Newton_Step;

end Double_Newton_Puiseux;
