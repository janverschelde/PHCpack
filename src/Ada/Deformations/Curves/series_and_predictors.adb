with Standard_Mathematical_Functions;
with Standard_Complex_Numbers;
with Standard_Complex_Vector_Norms;
with Standard_Series_Vector_Functions;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vector_Norms;
with DoblDobl_Series_Vector_Functions;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vector_Norms;
with QuadDobl_Series_Vector_Functions;

package body Series_and_Predictors is

  function Predicted_Error
             ( evls : Standard_Dense_Series_Vectors.Vector;
               step : double_float ) return double_float is

    eva : constant Standard_Complex_Vectors.Vector
        := Standard_Series_Vector_Functions.Eval(evls,step);
    res : constant double_float
        := Standard_Complex_Vector_Norms.Max_Norm(eva);

  begin
    return res;
  end Predicted_Error;

  function Predicted_Error
             ( evls : DoblDobl_Dense_Series_Vectors.Vector;
               step : double_double ) return double_double is

    eva : constant DoblDobl_Complex_Vectors.Vector
        := DoblDobl_Series_Vector_Functions.Eval(evls,step);
    res : constant double_double
        := DoblDobl_Complex_Vector_Norms.Max_Norm(eva);

  begin
    return res;
  end Predicted_Error;

  function Predicted_Error
             ( evls : QuadDobl_Dense_Series_Vectors.Vector;
               step : quad_double ) return quad_double is

    eva : constant QuadDobl_Complex_Vectors.Vector
        := QuadDobl_Series_Vector_Functions.Eval(evls,step);
    res : constant quad_double
        := QuadDobl_Complex_Vector_Norms.Max_Norm(eva);

  begin
    return res;
  end Predicted_Error;

  function Least_Order
             ( s : Standard_Dense_Series.Series; tol : double_float )
             return integer32 is
  begin
    for k in 0..s.order loop
      if Standard_Complex_Numbers.AbsVal(s.cff(k)) > tol
       then return k;
      end if;
    end loop;
    return s.order+1;
  end Least_Order;

  function Least_Order
             ( s : DoblDobl_Dense_Series.Series; tol : double_float )
             return integer32 is
  begin
    for k in 0..s.order loop
      if DoblDobl_Complex_Numbers.AbsVal(s.cff(k)) > tol
       then return k;
      end if;
    end loop;
    return s.order+1;
  end Least_Order;

  function Least_Order
             ( s : QuadDobl_Dense_Series.Series; tol : double_float )
             return integer32 is
  begin
    for k in 0..s.order loop
      if QuadDobl_Complex_Numbers.AbsVal(s.cff(k)) > tol
       then return k;
      end if;
    end loop;
    return s.order+1;
  end Least_Order;

  procedure Least_Order
             ( v : in Standard_Dense_Series_Vectors.Vector;
               tol : in double_float; vk,ord : out integer32 ) is

    res : integer32 := Least_Order(v(v'first),tol);
    vkord : integer32;

  begin
    vk := v'first;
    for k in v'first+1..v'last loop
      exit when (res = 0);
      vkord := Least_Order(v(k),tol);
      if vkord < res
       then res := vkord; vk := k;
      end if;
    end loop;
    ord := res;
  end Least_Order;

  procedure Least_Order
             ( v : in DoblDobl_Dense_Series_Vectors.Vector;
               tol : in double_float; vk,ord : out integer32 ) is

    res : integer32 := Least_Order(v(v'first),tol);
    vkord : integer32;

  begin
    vk := v'first;
    for k in v'first+1..v'last loop
      exit when (res = 0);
      vkord := Least_Order(v(k),tol);
      if vkord < res
       then res := vkord; vk := k;
      end if;
    end loop;
    ord := res;
  end Least_Order;

  procedure Least_Order
             ( v : in QuadDobl_Dense_Series_Vectors.Vector;
               tol : in double_float; vk,ord : out integer32 ) is

    res : integer32 := Least_Order(v(v'first),tol);
    vkord : integer32;

  begin
    vk := v'first;
    for k in v'first+1..v'last loop
      exit when (res = 0);
      vkord := Least_Order(v(k),tol);
      if vkord < res
       then res := vkord; vk := k;
      end if;
    end loop;
    ord := res;
  end Least_Order;

  function Set_Step_Size
             ( v : Standard_Dense_Series_Vectors.Vector;
               tolcff,tolres : double_float ) return double_float is

    vk,ord : integer32;
    res,valcff,pwr : double_float;

    use Standard_Mathematical_Functions;

  begin
    Least_Order(v,tolcff,vk,ord);
    valcff := Standard_Complex_Numbers.AbsVal(v(vk).cff(ord));
    pwr := 1.0/double_float(ord);
    res := (tolres/valcff)**pwr;
    return res;
  end Set_Step_Size;

  function Set_Step_Size
             ( v : DoblDobl_Dense_Series_Vectors.Vector;
               tolcff,tolres : double_float ) return double_float is

    vk,ord : integer32;
    dd_valcff : double_double;
    res,valcff,pwr : double_float;

    use Standard_Mathematical_Functions;

  begin
    Least_Order(v,tolcff,vk,ord);
    dd_valcff := DoblDobl_Complex_Numbers.AbsVal(v(vk).cff(ord));
    valcff := hi_part(dd_valcff);
    pwr := 1.0/double_float(ord);
    res := (tolres/valcff)**pwr;
    return res;
  end Set_Step_Size;

  function Set_Step_Size
             ( v : QuadDobl_Dense_Series_Vectors.Vector;
               tolcff,tolres : double_float ) return double_float is

    vk,ord : integer32;
    qd_valcff : quad_double;
    res,valcff,pwr : double_float;

    use Standard_Mathematical_Functions;

  begin
    Least_Order(v,tolcff,vk,ord);
    qd_valcff := QuadDobl_Complex_Numbers.AbsVal(v(vk).cff(ord));
    valcff := hihi_part(qd_valcff);
    pwr := 1.0/double_float(ord);
    res := (tolres/valcff)**pwr;
    return res;
  end Set_Step_Size;

  function Predicted_Solution
             ( srv : Standard_Dense_Series_Vectors.Vector;
               step : double_float )
             return Standard_Complex_Vectors.Vector is

    res : constant Standard_Complex_Vectors.Vector
        := Standard_Series_Vector_Functions.Eval(srv,step);

  begin
    return res;
  end Predicted_Solution;

  function Predicted_Solution
             ( srv : DoblDobl_Dense_Series_Vectors.Vector;
               step : double_double )
             return DoblDobl_Complex_Vectors.Vector is

    res : constant DoblDobl_Complex_Vectors.Vector
        := DoblDobl_Series_Vector_Functions.Eval(srv,step);

  begin
    return res;
  end Predicted_Solution;

  function Predicted_Solution
             ( srv : QuadDobl_Dense_Series_Vectors.Vector;
               step : quad_double )
             return QuadDobl_Complex_Vectors.Vector is

    res : constant QuadDobl_Complex_Vectors.Vector
        := QuadDobl_Series_Vector_Functions.Eval(srv,step);

  begin
    return res;
  end Predicted_Solution;

end Series_and_Predictors;
