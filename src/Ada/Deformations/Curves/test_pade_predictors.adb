with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs_io;        use QuadDobl_Complex_VecVecs_io;
with QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with Standard_Complex_Series_Vectors_io;
with DoblDobl_Complex_Series_Vectors_io;
with QuadDobl_Complex_Series_Vectors_io;
with Series_and_Homotopies;
with Test_Series_Predictors;
with Standard_Pade_Approximants_io;
with DoblDobl_Pade_Approximants_io;
with QuadDobl_Pade_Approximants_io;
with Homotopy_Pade_Approximants;

package body Test_Pade_Predictors is

  procedure Standard_Check_Prediction
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                pv : in Standard_Pade_Approximants.Pade_Vector;
                step : in double_float ) is

    phm : Standard_Complex_Poly_Systems.Poly_Sys(hom'range);
    valxpv : Standard_Complex_Vectors.Vector(phm'range);
    nrmxpv : double_float;
    xpv : constant Standard_Complex_Vectors.Vector(pv'range)
        := Standard_Pade_Approximants.Eval(pv,step);

  begin
    put_line("The evaluated Pade approximation :"); put_line(xpv);
    phm := Series_and_Homotopies.Eval(hom,step);
    valxpv := Standard_Complex_Poly_SysFun.Eval(phm,xpv);
    nrmxpv := Standard_Complex_Vector_Norms.Max_Norm(valxpv);
    put("The actual error of the Pade predictor : "); put(nrmxpv,3); new_line;
    Standard_Complex_Poly_Systems.Clear(phm);
  end Standard_Check_Prediction;

  procedure DoblDobl_Check_Prediction
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector;
                step : in double_double ) is

    phm : DoblDobl_Complex_Poly_Systems.Poly_Sys(hom'range);
    valxpv : DoblDobl_Complex_Vectors.Vector(phm'range);
    nrmxpv : double_double;
    xpv : constant DoblDobl_Complex_Vectors.Vector(pv'range)
        := DoblDobl_Pade_Approximants.Eval(pv,step);

  begin
    put_line("The evaluated Pade approximant :"); put_line(xpv);
    phm := Series_and_Homotopies.Eval(hom,step);
    valxpv := DoblDobl_Complex_Poly_SysFun.Eval(phm,xpv);
    nrmxpv := DoblDobl_Complex_Vector_Norms.Max_Norm(valxpv);
    put("The actual error of the Pade predictor : "); put(nrmxpv,3); new_line;
    DoblDobl_Complex_Poly_Systems.Clear(phm);
  end DoblDobl_Check_Prediction;

  procedure QuadDobl_Check_Prediction
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                pv : in QuadDobl_Pade_Approximants.Pade_Vector;
                step : in quad_double ) is

    phm : QuadDobl_Complex_Poly_Systems.Poly_Sys(hom'range);
    valxpv : QuadDobl_Complex_Vectors.Vector(phm'range);
    nrmxpv : quad_double;
    xpv : constant QuadDobl_Complex_Vectors.Vector(pv'range)
        := QuadDobl_Pade_Approximants.Eval(pv,step);

  begin
    put_line("The evaluated Pade approximant :"); put_line(xpv);
    phm := Series_and_Homotopies.Eval(hom,step);
    valxpv := QuadDobl_Complex_Poly_SysFun.Eval(phm,xpv);
    nrmxpv := QuadDobl_Complex_Vector_Norms.Max_Norm(valxpv);
    put("The actual error of the Pade predictor : "); put(nrmxpv,3); new_line;
    QuadDobl_Complex_Poly_Systems.Clear(phm);
  end QuadDobl_Check_Prediction;

  procedure Standard_Step_Prediction
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                srv,eva : in Standard_Complex_Series_Vectors.Vector;
                pv : in Standard_Pade_Approximants.Pade_Vector ) is

    step : double_float := 0.0;
    ans : character;

  begin
    new_line;
    put_line("Interactive loop to experiment with the step size ...");
    loop
      put("Give the step size : "); get(step);
      Test_Series_Predictors.Standard_Check_Prediction(hom,srv,eva,step);
      Standard_Check_Prediction(hom,pv,step);
      put("Check another step size ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Step_Prediction;

  procedure DoblDobl_Step_Prediction
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                srv,eva : in DoblDobl_Complex_Series_Vectors.Vector;
                pv : in DoblDobl_Pade_Approximants.Pade_Vector ) is

    step : double_double := create(0.0);
    ans : character;

  begin
    new_line;
    put_line("Interactive loop to experiment with the step size ...");
    skip_line; -- needed for the prompting of the double double step
    loop
      put("Give the step size : "); get(step);
      Test_Series_Predictors.DoblDobl_Check_Prediction(hom,srv,eva,step);
      DoblDobl_Check_Prediction(hom,pv,step);
      put("Check another step size ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end DoblDobl_Step_Prediction;

  procedure QuadDobl_Step_Prediction
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                srv,eva : in QuadDobl_Complex_Series_Vectors.Vector;
                pv : QuadDobl_Pade_Approximants.Pade_Vector ) is

    step : quad_double := create(0.0);
    ans : character;

  begin
    new_line;
    put_line("Interactive loop to experiment with the step size ...");
    skip_line; -- needed for prompting of the quad double step
    loop
      put("Give the step size : "); get(step);
      Test_Series_Predictors.QuadDobl_Check_Prediction(hom,srv,eva,step);
      QuadDobl_Check_Prediction(hom,pv,step);
      put("Check another step size ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end QuadDobl_Step_Prediction;

  procedure Forward_Pole_Radius
              ( neq,degnum,degden : in integer32;
                solv : in Standard_Complex_Vectors.Vector;
                srv : out Standard_Complex_Series_Vectors.Vector;
                eva : out Standard_Complex_Series_Vectors.Vector;
                pv : out Standard_Pade_Approximants.Pade_Vector;
                poles : out Standard_Complex_VecVecs.VecVec;
                fpr : out double_float; verbose : in boolean := false ) is

    nbt : constant natural32 := natural32(degnum+degden+1);
    lead,idx : integer32;

  begin
    Homotopy_Pade_Approximants.Standard_Pade_Approximant
      (solv,neq+1,neq,degnum,degden,nbt,srv,eva,pv,verbose);
    if verbose then
      put_line("The solution series :");
      Standard_Complex_Series_Vectors_io.put(srv);
      put_line("The evaluated solution series :");
      Standard_Complex_Series_Vectors_io.put(eva);
      put_line("The Pade approximant :");
      for i in pv'range loop
        put_line(Standard_Pade_Approximants_io.Write(pv(i)));
      end loop;
    end if;
    poles := Homotopy_Pade_Approximants.Standard_Poles(pv);
    if verbose
     then put_line("The poles : "); put_line(poles);
    end if;
    Homotopy_Pade_Approximants.Closest_Pole(poles,lead,idx,fpr);
    if verbose then
      put("The radius of the smallest pole : ");
      put(fpr,3); new_line;
    end if;
  end Forward_Pole_Radius;

  procedure Forward_Pole_Radius
              ( neq,degnum,degden : in integer32;
                solv : in DoblDobl_Complex_Vectors.Vector;
                srv : out DoblDobl_Complex_Series_Vectors.Vector;
                eva : out DoblDobl_Complex_Series_Vectors.Vector;
                pv : out DoblDobl_Pade_Approximants.Pade_Vector;
                poles : out DoblDobl_Complex_VecVecs.VecVec;
                fpr : out double_double; verbose : in boolean := false ) is

    nbt : constant natural32 := natural32(degnum+degden+1);
    lead,idx : integer32;

  begin
    Homotopy_Pade_Approximants.DoblDobl_Pade_Approximant
      (solv,neq+1,neq,degnum,degden,nbt,srv,eva,pv);
    if verbose then
      put_line("The solution series :");
      DoblDobl_Complex_Series_Vectors_io.put(srv);
      put_line("The evaluated solution series :");
      DoblDobl_Complex_Series_Vectors_io.put(eva);
      put_line("The Pade approximant :");
      for i in pv'range loop
        put_line(DoblDobl_Pade_Approximants_io.Write(pv(i)));
      end loop;
    end if;
    poles := Homotopy_Pade_Approximants.DoblDobl_Poles(pv);
    if verbose
     then put_line("The poles : "); put_line(poles);
    end if;
    Homotopy_Pade_Approximants.Closest_Pole(poles,lead,idx,fpr);
    if verbose then
      put("The radius of the smallest pole : ");
      put(fpr,3); new_line;
    end if;
  end Forward_Pole_Radius;

  procedure Forward_Pole_Radius
              ( neq,degnum,degden : in integer32;
                solv : in QuadDobl_Complex_Vectors.Vector;
                srv : out QuadDobl_Complex_Series_Vectors.Vector;
                eva : out QuadDobl_Complex_Series_Vectors.Vector;
                pv : out QuadDobl_Pade_Approximants.Pade_Vector;
                poles : out QuadDobl_Complex_VecVecs.VecVec;
                fpr : out quad_double; verbose : in boolean := false ) is

    nbt : constant natural32 := natural32(degnum+degden+1);
    lead,idx : integer32;

  begin
    Homotopy_Pade_Approximants.QuadDobl_Pade_Approximant
      (solv,neq+1,neq,degnum,degden,nbt,srv,eva,pv);
    if verbose then
      put_line("The solution series :");
      QuadDobl_Complex_Series_Vectors_io.put(srv);
      put_line("The evaluated solution series :");
      QuadDobl_Complex_Series_Vectors_io.put(eva);
      put_line("The Pade approximant :");
      for i in pv'range loop
        put_line(QuadDobl_Pade_Approximants_io.Write(pv(i)));
      end loop;
    end if;
    poles := Homotopy_Pade_Approximants.QuadDobl_Poles(pv);
    if verbose
     then put_line("The poles : "); put_line(poles);
    end if;
    Homotopy_Pade_Approximants.Closest_Pole(poles,lead,idx,fpr);
    if verbose then
      put("The radius of the smallest pole : ");
      put(fpr,3); new_line;
    end if;
  end Forward_Pole_Radius;

  function Forward_Pole_Radius
              ( neq,degnum,degden : integer32;
                solv : Standard_Complex_Vectors.Vector;
                verbose : boolean := false ) return double_float is

    res : double_float;
    srv : Standard_Complex_Series_Vectors.Vector(solv'range);
    eva : Standard_Complex_Series_Vectors.Vector(1..neq);
    pv : Standard_Pade_Approximants.Pade_Vector(solv'range);
    poles : Standard_Complex_VecVecs.VecVec(solv'range);

  begin
    Forward_Pole_Radius
      (neq,degnum,degden,solv,srv,eva,pv,poles,res,verbose);
    Standard_Pade_Approximants.Clear(pv);
    Standard_Complex_Series_Vectors.Clear(srv);
    Standard_Complex_Series_Vectors.Clear(eva);
    Standard_Complex_VecVecs.Clear(poles);
    return res;
  end Forward_Pole_Radius;

  function Forward_Pole_Radius
              ( neq,degnum,degden : integer32;
                solv : DoblDobl_Complex_Vectors.Vector;
                verbose : boolean := false ) return double_double is

    res : double_double;
    srv : DoblDobl_Complex_Series_Vectors.Vector(solv'range);
    eva : DoblDobl_Complex_Series_Vectors.Vector(1..neq);
    pv : DoblDobl_Pade_Approximants.Pade_Vector(solv'range);
    poles : DoblDobl_Complex_VecVecs.VecVec(solv'range);

  begin
    Forward_Pole_Radius
      (neq,degnum,degden,solv,srv,eva,pv,poles,res,verbose);
    DoblDobl_Pade_Approximants.Clear(pv);
    DoblDobl_Complex_Series_Vectors.Clear(srv);
    DoblDobl_Complex_Series_Vectors.Clear(eva);
    DoblDobl_Complex_VecVecs.Clear(poles);
    return res;
  end Forward_Pole_Radius;

  function Forward_Pole_Radius
              ( neq,degnum,degden : integer32;
                solv : QuadDobl_Complex_Vectors.Vector;
                verbose : boolean := false ) return quad_double is

    res : quad_double;
    srv : QuadDobl_Complex_Series_Vectors.Vector(solv'range);
    eva : QuadDobl_Complex_Series_Vectors.Vector(1..neq);
    pv : QuadDobl_Pade_Approximants.Pade_Vector(solv'range);
    poles : QuadDobl_Complex_VecVecs.VecVec(solv'range);

  begin
    Forward_Pole_Radius
      (neq,degnum,degden,solv,srv,eva,pv,poles,res,verbose);
    QuadDobl_Pade_Approximants.Clear(pv);
    QuadDobl_Complex_Series_Vectors.Clear(srv);
    QuadDobl_Complex_Series_Vectors.Clear(eva);
    QuadDobl_Complex_VecVecs.Clear(poles);
    return res;
  end Forward_Pole_Radius;

  procedure Forward_Pole_Radius
              ( neq,degnum,degden : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List;
                radii : out Standard_Floating_Vectors.Vector;
                fpr : out double_float; verbose : in boolean := false ) is

    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    idx : integer32 := 0;

  begin
    fpr := -1.0;
    while not Standard_Complex_Solutions.Is_Null(tmp) loop
      idx := idx + 1;
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      radii(idx) := Forward_Pole_Radius(neq,degnum,degden,ls.v,verbose);
      if fpr = -1.0 then
        fpr := radii(idx);
      elsif radii(idx) >= 0.0 then
        if radii(idx) < fpr
         then fpr := radii(idx);
        end if;
      end if;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Forward_Pole_Radius;

  procedure Forward_Pole_Radius
              ( neq,degnum,degden : in integer32;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                radii : out Double_Double_Vectors.Vector;
                fpr : out double_double; verbose : in boolean := false ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    idx : integer32 := 0;
    minone : constant double_double := create(-1.0);

  begin
    fpr := minone;
    while not DoblDobl_Complex_Solutions.Is_Null(tmp) loop
      idx := idx + 1;
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      radii(idx) := Forward_Pole_Radius(neq,degnum,degden,ls.v,verbose);
      if fpr = minone then
        fpr := radii(idx);
      elsif radii(idx) >= 0.0 then
        if radii(idx) < fpr
         then fpr := radii(idx);
        end if;
      end if;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Forward_Pole_Radius;

  procedure Forward_Pole_Radius
              ( neq,degnum,degden : in integer32;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                radii : out Quad_Double_Vectors.Vector;
                fpr : out quad_double; verbose : in boolean := false ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    idx : integer32 := 0;
    minone : constant quad_double := create(-1.0);

  begin
    fpr := minone;
    while not QuadDobl_Complex_Solutions.Is_Null(tmp) loop
      idx := idx + 1;
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      radii(idx) := Forward_Pole_Radius(neq,degnum,degden,ls.v,verbose);
      if fpr = minone then
        fpr := radii(idx);
      elsif radii(idx) >= 0.0 then
        if radii(idx) < fpr
         then fpr := radii(idx);
        end if;
      end if;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Forward_Pole_Radius;

  procedure Standard_Test_Pade_Prediction
              ( neq : in integer32;
                sol : in Standard_Complex_Solutions.Solution;
                pv : out Standard_Pade_Approximants.Pade_Vector ) is

    degnum,degden : integer32 := 0;
    srv : Standard_Complex_Series_Vectors.Vector(sol.v'range);
    eva : Standard_Complex_Series_Vectors.Vector(1..neq);
    poles : Standard_Complex_VecVecs.VecVec(sol.v'range);
    minpole : double_float;

  begin
    new_line;
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    Forward_Pole_Radius
      (neq,degnum,degden,sol.v,srv,eva,pv,poles,minpole,true);
  end Standard_Test_Pade_Prediction;

  procedure DoblDobl_Test_Pade_Prediction
              ( neq : in integer32;
                sol : in DoblDobl_Complex_Solutions.Solution;
                pv : out DoblDobl_Pade_Approximants.Pade_Vector ) is

    degnum,degden : integer32 := 0;
    srv : DoblDobl_Complex_Series_Vectors.Vector(sol.v'range);
    eva : DoblDobl_Complex_Series_Vectors.Vector(1..neq);
    poles : DoblDobl_Complex_VecVecs.VecVec(pv'range);
    minpole : double_double;

  begin
    new_line;
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    Forward_Pole_Radius
      (neq,degnum,degden,sol.v,srv,eva,pv,poles,minpole,true);
  end DoblDobl_Test_Pade_Prediction;

  procedure QuadDobl_Test_Pade_Prediction
              ( neq : in integer32;
                sol : in QuadDobl_Complex_Solutions.Solution;
                pv : out QuadDobl_Pade_Approximants.Pade_Vector ) is

    degnum,degden : integer32 := 0;
    srv : QuadDobl_Complex_Series_Vectors.Vector(sol.v'range);
    eva : QuadDobl_Complex_Series_Vectors.Vector(1..neq);
    poles : QuadDobl_Complex_VecVecs.VecVec(pv'range);
    minpole : quad_double;

  begin
    new_line;
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    Forward_Pole_Radius
      (neq,degnum,degden,sol.v,srv,eva,pv,poles,minpole,true);
  end QuadDobl_Test_Pade_Prediction;

end Test_Pade_Predictors;
