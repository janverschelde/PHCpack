with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with Standard_Complex_Series_Functions;
with Standard_Complex_Series_Vectors_io;
with Standard_CSeries_Vector_Functions;
with Standard_CSeries_Poly_SysFun;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_Series_Functions;
with DoblDobl_Complex_Series_Vectors_io;
with DoblDobl_CSeries_Vector_Functions;
with DoblDobl_CSeries_Poly_SysFun;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Series_Functions;
with QuadDobl_Complex_Series_Vectors_io;
with QuadDobl_CSeries_Vector_Functions;
with QuadDobl_CSeries_Poly_SysFun;
with Complex_Series_and_Polynomials_io;
with Power_Series_Methods;               use Power_Series_Methods;
with Series_and_Solutions;

package body Series_and_Predictors is

  procedure Newton_Prediction
              ( maxdeg,nit : in integer32;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in Standard_Complex_Vectors.Vector;
                srv : out Standard_Complex_Series_Vectors.Vector;
                eva : out Standard_Complex_Series_Vectors.Vector ) is

    info : integer32;

  begin
    srv := Series_and_Solutions.Create(sol,0);
    if hom'last = sol'last
     then Run_LU_Newton(maxdeg,nit,hom,srv,info);
     else Run_QR_Newton(maxdeg,nit,hom,srv);
    end if;
    eva := Standard_CSeries_Poly_SysFun.Eval(hom,srv);
  end Newton_Prediction;

  procedure Newton_Prediction
              ( maxdeg,nit : in integer32;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Vectors.Vector;
                srv : out DoblDobl_Complex_Series_Vectors.Vector;
                eva : out DoblDobl_Complex_Series_Vectors.Vector ) is

    info : integer32;

  begin
    srv := Series_and_Solutions.Create(sol,0);
    if hom'last = sol'last
     then Run_LU_Newton(maxdeg,nit,hom,srv,info);
     else Run_QR_Newton(maxdeg,nit,hom,srv);
    end if;
    eva := DoblDobl_CSeries_Poly_SysFun.Eval(hom,srv);
  end Newton_Prediction;

  procedure Newton_Prediction
              ( maxdeg,nit : in integer32;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in QuadDobl_Complex_Vectors.Vector;
                srv : out QuadDobl_Complex_Series_Vectors.Vector;
                eva : out QuadDobl_Complex_Series_Vectors.Vector ) is

    info : integer32;

  begin
    srv := Series_and_Solutions.Create(sol,0);
    if hom'last = sol'last
     then Run_LU_Newton(maxdeg,nit,hom,srv,info);
     else Run_QR_Newton(maxdeg,nit,hom,srv);
    end if;
    eva := QuadDobl_CSeries_Poly_SysFun.Eval(hom,srv);
  end Newton_Prediction;

  procedure Newton_Prediction
              ( file : in file_type; maxdeg,nit : in integer32;
                hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sol : in Standard_Complex_Vectors.Vector;
                srv : out Standard_Complex_Series_Vectors.Vector;
                eva : out Standard_Complex_Series_Vectors.Vector;
                verbose : in boolean := false ) is

    info : integer32;

  begin
    srv := Series_and_Solutions.Create(sol,0);
    if verbose then
      put_line(file,"The solution :");
      Standard_Complex_Vectors_io.put_line(file,sol);
      new_line(file);
      put_line(file,"The series for the solution :");
      Standard_Complex_Series_Vectors_io.put(file,srv);
    end if;
    if hom'last = sol'last then
      if not verbose then
        Run_LU_Newton(maxdeg,nit,hom,srv,info);
      else
        put_line(file,"Applying LU Newton ...");
        Run_LU_Newton(file,maxdeg,nit,hom,srv,info,true);
      end if;
    else
      if not verbose then
        Run_QR_Newton(maxdeg,nit,hom,srv);
      else
        put_line(file,"Applying QR Newton ...");
        Run_QR_Newton(file,maxdeg,nit,hom,srv,true);
      end if;
    end if;
    eva := Standard_CSeries_Poly_SysFun.Eval(hom,srv);
    if verbose then
      put_line(file,"The evaluated series : ");
      Complex_Series_and_Polynomials_io.put(file,eva);
    end if;
  end Newton_Prediction;

  procedure Newton_Prediction
              ( file : in file_type; maxdeg,nit : in integer32;
                hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in DoblDobl_Complex_Vectors.Vector;
                srv : out DoblDobl_Complex_Series_Vectors.Vector;
                eva : out DoblDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false ) is

    info : integer32;

  begin
    srv := Series_and_Solutions.Create(sol,0);
    if verbose then
      put_line(file,"The solution :");
      DoblDobl_Complex_Vectors_io.put_line(file,sol);
      new_line(file);
      put_line(file,"The series for the solution :");
      DoblDobl_Complex_Series_Vectors_io.put(file,srv);
    end if;
    if hom'last = sol'last then
      if not verbose then
        Run_LU_Newton(maxdeg,nit,hom,srv,info);
      else
        put_line(file,"Applying LU Newton ...");
        Run_LU_Newton(file,maxdeg,nit,hom,srv,info,true);
      end if;
    else
      if not verbose then
        Run_QR_Newton(maxdeg,nit,hom,srv);
      else
        put_line(file,"Applying QR Newton ...");
        Run_QR_Newton(file,maxdeg,nit,hom,srv,true);
      end if;
    end if;
    eva := DoblDobl_CSeries_Poly_SysFun.Eval(hom,srv);
    if verbose then
      put_line(file,"The evaluated series : ");
      Complex_Series_and_Polynomials_io.put(file,eva);
    end if;
  end Newton_Prediction;

  procedure Newton_Prediction
              ( file : in file_type; maxdeg,nit : in integer32;
                hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sol : in QuadDobl_Complex_Vectors.Vector;
                srv : out QuadDobl_Complex_Series_Vectors.Vector;
                eva : out QuadDobl_Complex_Series_Vectors.Vector;
                verbose : in boolean := false ) is

    info : integer32;

  begin
    srv := Series_and_Solutions.Create(sol,0);
    if verbose then
      put_line(file,"The solution :");
      QuadDobl_Complex_Vectors_io.put_line(file,sol);
      new_line(file);
      put_line(file,"The series for the solution :");
      QuadDobl_Complex_Series_Vectors_io.put(file,srv);
    end if;
    if hom'last = sol'last then
      if not verbose then
        Run_LU_Newton(maxdeg,nit,hom,srv,info);
      else
        put_line(file,"Applying LU Newton ...");
        Run_LU_Newton(file,maxdeg,nit,hom,srv,info,true);
      end if;
    else
      if not verbose then
        Run_QR_Newton(maxdeg,nit,hom,srv);
      else
        put_line(file,"Applying QR Newton ...");
        Run_QR_Newton(file,maxdeg,nit,hom,srv,true);
      end if;
    end if;
    eva := QuadDobl_CSeries_Poly_SysFun.Eval(hom,srv);
    if verbose then
      put_line(file,"The evaluated series : ");
      Complex_Series_and_Polynomials_io.put(file,eva);
    end if;
  end Newton_Prediction;

  function Predicted_Error
             ( evls : Standard_Complex_Series_Vectors.Vector;
               step : double_float ) return double_float is

    eva : constant Standard_Complex_Vectors.Vector
        := Standard_CSeries_Vector_Functions.Eval(evls,step);
    res : constant double_float
        := Standard_Complex_Vector_Norms.Max_Norm(eva);

  begin
    return res;
  end Predicted_Error;

  function Predicted_Error
             ( evls : DoblDobl_Complex_Series_Vectors.Vector;
               step : double_double ) return double_double is

    eva : constant DoblDobl_Complex_Vectors.Vector
        := DoblDobl_CSeries_Vector_Functions.Eval(evls,step);
    res : constant double_double
        := DoblDobl_Complex_Vector_Norms.Max_Norm(eva);

  begin
    return res;
  end Predicted_Error;

  function Predicted_Error
             ( evls : QuadDobl_Complex_Series_Vectors.Vector;
               step : quad_double ) return quad_double is

    eva : constant QuadDobl_Complex_Vectors.Vector
        := QuadDobl_CSeries_Vector_Functions.Eval(evls,step);
    res : constant quad_double
        := QuadDobl_Complex_Vector_Norms.Max_Norm(eva);

  begin
    return res;
  end Predicted_Error;

  procedure Order ( v : in Standard_Complex_Series_Vectors.Vector;
                    tol : in double_float; vk,ord : out integer32 ) is

    res : integer32
        := Standard_Complex_Series_Functions.Order(v(v'first),tol);
    vkord : integer32;

  begin
    vk := v'first;
    for k in v'first+1..v'last loop
      exit when (res = 0);
      vkord := Standard_Complex_Series_Functions.Order(v(k),tol);
      if vkord < res
       then res := vkord; vk := k;
      end if;
    end loop;
    ord := res;
  end Order;

  procedure Order ( v : in DoblDobl_Complex_Series_Vectors.Vector;
                    tol : in double_float; vk,ord : out integer32 ) is

    res : integer32
        := DoblDobl_Complex_Series_Functions.Order(v(v'first),tol);
    vkord : integer32;

  begin
    vk := v'first;
    for k in v'first+1..v'last loop
      exit when (res = 0);
      vkord := DoblDobl_Complex_Series_Functions.Order(v(k),tol);
      if vkord < res
       then res := vkord; vk := k;
      end if;
    end loop;
    ord := res;
  end Order;

  procedure Order ( v : in QuadDobl_Complex_Series_Vectors.Vector;
                    tol : in double_float; vk,ord : out integer32 ) is

    res : integer32
        := QuadDobl_Complex_Series_Functions.Order(v(v'first),tol);
    vkord : integer32;

  begin
    vk := v'first;
    for k in v'first+1..v'last loop
      exit when (res = 0);
      vkord := QuadDobl_Complex_Series_Functions.Order(v(k),tol);
      if vkord < res
       then res := vkord; vk := k;
      end if;
    end loop;
    ord := res;
  end Order;

  function Set_Step_Size
             ( v : Standard_Complex_Series_Vectors.Vector;
               tolcff,tolres : double_float ) return double_float is

    vk,ord : integer32;
    res,valcff,pwr,arg : double_float;

    use Standard_Mathematical_Functions;

  begin
    Order(v,tolcff,vk,ord);
    valcff := Standard_Complex_Numbers.AbsVal(v(vk).cff(ord));
    arg := tolres/valcff;
    if ord = 0 then
      res := arg;
    else
      pwr := 1.0/double_float(ord);
      res := arg**pwr;
    end if;
    return res;
  end Set_Step_Size;

  function Set_Step_Size
             ( v : DoblDobl_Complex_Series_Vectors.Vector;
               tolcff,tolres : double_float ) return double_float is

    vk,ord : integer32;
    dd_valcff : double_double;
    res,valcff,pwr,arg : double_float;

    use Standard_Mathematical_Functions;

  begin
    Order(v,tolcff,vk,ord);
    dd_valcff := DoblDobl_Complex_Numbers.AbsVal(v(vk).cff(ord));
    valcff := hi_part(dd_valcff);
    arg := tolres/valcff;
    if ord = 0 then
      res := arg;
    else
      pwr := 1.0/double_float(ord);
      res := arg**pwr;
    end if;
    return res;
  end Set_Step_Size;

  function Set_Step_Size
             ( v : QuadDobl_Complex_Series_Vectors.Vector;
               tolcff,tolres : double_float ) return double_float is

    vk,ord : integer32;
    qd_valcff : quad_double;
    res,valcff,pwr,arg : double_float;

    use Standard_Mathematical_Functions;

  begin
    Order(v,tolcff,vk,ord);
    qd_valcff := QuadDobl_Complex_Numbers.AbsVal(v(vk).cff(ord));
    valcff := hihi_part(qd_valcff);
    arg := tolres/valcff;
    if ord = 0 then
      res := arg;
    else
      pwr := 1.0/double_float(ord);
      res := arg**pwr;
    end if;
    return res;
  end Set_Step_Size;

  function Set_Step_Size
             ( file : file_type;
               v : Standard_Complex_Series_Vectors.Vector;
               tolcff,tolres : double_float;
               verbose : boolean := false ) return double_float is

    vk,ord : integer32;
    res,valcff,pwr,arg : double_float;

    use Standard_Mathematical_Functions;

  begin
    Order(v,tolcff,vk,ord);
    valcff := Standard_Complex_Numbers.AbsVal(v(vk).cff(ord));
    arg := tolres/valcff;
    if verbose then
      put(file,"order : "); put(file,ord,1);
      put(file," at component : "); put(file,vk,1);
      put(file," arg = "); put(file,arg); new_line(file);
    end if;
    if ord = 0 then
      res := arg;
    else
      pwr := 1.0/double_float(ord);
      res := arg**pwr;
    end if;
    return res;
  end Set_Step_Size;

  function Set_Step_Size
             ( file : file_type;
               v : DoblDobl_Complex_Series_Vectors.Vector;
               tolcff,tolres : double_float;
               verbose : boolean := false ) return double_float is

    vk,ord : integer32;
    dd_valcff : double_double;
    res,valcff,pwr,arg : double_float;

    use Standard_Mathematical_Functions;

  begin
    Order(v,tolcff,vk,ord);
    dd_valcff := DoblDobl_Complex_Numbers.AbsVal(v(vk).cff(ord));
    valcff := hi_part(dd_valcff);
    arg := tolres/valcff;
    if verbose then
      put(file,"order : "); put(file,ord,1);
      put(file," at component : "); put(file,vk,1);
      put(file," arg = "); put(file,arg); new_line(file);
    end if;
    if ord = 0 then
      res := arg;
    else
      pwr := 1.0/double_float(ord);
      res := arg**pwr;
    end if;
    return res;
  end Set_Step_Size;

  function Set_Step_Size
             ( file : file_type;
               v : QuadDobl_Complex_Series_Vectors.Vector;
               tolcff,tolres : double_float;
               verbose : boolean := false ) return double_float is

    vk,ord : integer32;
    qd_valcff : quad_double;
    res,valcff,pwr,arg : double_float;

    use Standard_Mathematical_Functions;

  begin
    Order(v,tolcff,vk,ord);
    qd_valcff := QuadDobl_Complex_Numbers.AbsVal(v(vk).cff(ord));
    valcff := hihi_part(qd_valcff);
    arg := tolres/valcff;
    if verbose then
      put(file,"order : "); put(file,ord,1);
      put(file," at component : "); put(file,vk,1);
      put(file," arg = "); put(file,arg); new_line(file);
    end if;
    if ord = 0 then
      res := arg;
    else
      pwr := 1.0/double_float(ord);
      res := arg**pwr;
    end if;
    return res;
  end Set_Step_Size;

  function Predicted_Solution
             ( srv : Standard_Complex_Series_Vectors.Vector;
               step : double_float )
             return Standard_Complex_Vectors.Vector is

    res : constant Standard_Complex_Vectors.Vector
        := Standard_CSeries_Vector_Functions.Eval(srv,step);

  begin
    return res;
  end Predicted_Solution;

  function Predicted_Solution
             ( srv : DoblDobl_Complex_Series_Vectors.Vector;
               step : double_double )
             return DoblDobl_Complex_Vectors.Vector is

    res : constant DoblDobl_Complex_Vectors.Vector
        := DoblDobl_CSeries_Vector_Functions.Eval(srv,step);

  begin
    return res;
  end Predicted_Solution;

  function Predicted_Solution
             ( srv : QuadDobl_Complex_Series_Vectors.Vector;
               step : quad_double )
             return QuadDobl_Complex_Vectors.Vector is

    res : constant QuadDobl_Complex_Vectors.Vector
        := QuadDobl_CSeries_Vector_Functions.Eval(srv,step);

  begin
    return res;
  end Predicted_Solution;

end Series_and_Predictors;
