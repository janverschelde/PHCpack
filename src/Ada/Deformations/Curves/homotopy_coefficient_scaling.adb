with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors_io;
with Standard_Coefficient_Homotopy;
with Standard_Complex_Series_Vectors_io;
with Standard_CSeries_Vector_Functions;
with Series_and_Solutions;
with Hyperplane_Solution_Scaling;

package body Homotopy_Coefficient_Scaling is

  procedure Last_Coefficients
              ( file : in file_type;
                fcf : in Standard_Complex_Series_Vectors.Link_to_Vector;
                t : in double_float;
                gamma : in Standard_Complex_Numbers.Complex_Number ) is

    nbreqs : constant integer32 
           := Standard_Coefficient_Homotopy.Number_of_Equations;
    hcp,hcq : Standard_Complex_Vectors.Link_to_Vector;

  begin
    put_line(file,"The last of wrk_fcf :");
    Standard_Complex_Series_Vectors_io.put_line(file,fcf);
    put(file,"Number of equations in the coefficient homotopy : ");
    put(file,nbreqs,1); new_line(file);
    if nbreqs > 0 then
      hcp := Standard_Coefficient_Homotopy.Target_Coefficients(nbreqs);
      hcq := Standard_Coefficient_Homotopy.Start_Coefficients(nbreqs);
      put_line(file,"Coefficients of the last target equation :");
      Standard_Complex_Vectors_io.put_line(file,hcp);
      put_line(file,"Coefficients of the last start equation :");
      Standard_Complex_Vectors_io.put_line(file,hcq);
      declare
        wrk : Standard_Complex_Vectors.Vector(hcp'range);
        onemint : constant double_float := 1.0 - t;
        use Standard_Complex_Numbers;
      begin
        for i in wrk'range loop
          wrk(i) := t*hcp(i);
        end loop;
       -- the constant is the last coefficient in the vector
        wrk(wrk'last) := wrk(wrk'last) + gamma*onemint*hcq(hcq'last);
       -- the Z0 coefficient is the next to last coefficient
        wrk(wrk'last-1) := wrk(wrk'last-1) + gamma*onemint*hcq(hcq'last-1);
        put_line(file,"last coefficients with t value multiplied in :");
        Standard_Complex_Vectors_io.put_line(file,wrk);
        put_line(file,"the coefficients in fcf :");
        for i in fcf'range loop
          Standard_Complex_Numbers_io.put(file,fcf(i).cff(0));
          new_line(file);
        end loop;
      end;
    end if;
  end Last_Coefficients;

  procedure Scale_Solution_Coefficients
              ( file : in file_type;
                fhm : Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys;
                hcf : in Standard_Complex_Series_VecVecs.VecVec;
                sol : in out Standard_Complex_Vectors.Vector;
                t : in double_float;
                gamma : in Standard_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   Scales the solution and the coefficients in the homotopy.

    use Standard_Complex_Numbers;

    prd : Complex_Number;
    hcp,hcq : Standard_Complex_Vectors.Link_to_Vector;
    nbreqs : constant integer32 
           := Standard_Coefficient_Homotopy.Number_of_Equations;
    srv : Standard_Complex_Series_Vectors.Vector(sol'range);
    eva : Standard_Complex_Series_Vectors.Vector(fhm'range);
    evh : Standard_Complex_Vectors.Vector(fhm'range);
    fcf : constant Standard_Complex_Series_Vectors.Link_to_Vector
        := hcf(hcf'last);
    cmt : constant Complex_Number := Create(t);

  begin
    if nbreqs > 0 then
      srv := Series_and_Solutions.Create(sol,0);
      eva := Standard_CSeries_Poly_SysFun.Eval(fhm,hcf,srv);
      evh := Standard_Coefficient_Homotopy.Eval(sol,cmt);
      put_line(file,"The evaluated solution before scaling :");
      Standard_Complex_Vectors_io.put_line(file,evh);
      put_line(file,"The gamma constant : ");
      put(file,gamma); new_line(file);
      put_line(file,"The evaluated solution series before scaling :");
      Standard_Complex_Series_Vectors_io.put_line(file,eva);
      Hyperplane_Solution_Scaling.Scale(sol);
      hcq := Standard_Coefficient_Homotopy.Start_Coefficients(nbreqs);
      hcp := Standard_Coefficient_Homotopy.Target_Coefficients(nbreqs);
      hcq(hcq'last) := -sol(sol'last);
      if t = 0.0 then
        fcf(fcf'last).cff(0) := -gamma*sol(sol'last);
        fcf(fcf'last).cff(1) := gamma*sol(sol'last) + hcp(hcp'last);
      else
        prd := hcp(hcp'first)*sol(sol'first);
        for i in sol'first+1..sol'last loop
          prd := prd + hcp(i)*sol(i);
        end loop;
        hcp(hcp'last) := -prd;
        for i in fcf'range loop
          fcf(i).cff(0) := Standard_Complex_Numbers.create(0.0);
          fcf(i).cff(1) := hcp(i);
        end loop;
        fcf(fcf'last-1).cff(0) := gamma;
        fcf(fcf'last-1).cff(1) := fcf(fcf'last-1).cff(1) - gamma;
        fcf(fcf'last).cff(0) := -gamma*sol(sol'last);
        fcf(fcf'last).cff(1) := fcf(fcf'last).cff(1) + gamma*sol(sol'last);
        Standard_CSeries_Vector_Functions.Shift(fcf.all,-t);
      end if;
      srv := Series_and_Solutions.Create(sol,0);
      eva := Standard_CSeries_Poly_SysFun.Eval(fhm,hcf,srv);
      put_line(file,"The evaluated solution series after scaling :");
      Standard_Complex_Series_Vectors_io.put_line(file,eva);
      evh := Standard_Coefficient_Homotopy.Eval(sol,cmt);
      put_line(file,"The evaluated solution after scaling :");
      Standard_Complex_Vectors_io.put_line(file,evh);
    end if;
  end Scale_Solution_Coefficients;

end Homotopy_Coefficient_Scaling;
