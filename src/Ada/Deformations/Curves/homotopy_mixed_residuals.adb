with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers_Polar;
with DoblDobl_Complex_Numbers_Polar;
with QuadDobl_Complex_Numbers_Polar;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Homotopy;
with Standard_Coefficient_Homotopy;
with DoblDobl_Homotopy;
with DoblDobl_Coefficient_Homotopy;
with QuadDobl_Homotopy;
with QuadDobl_Coefficient_Homotopy;
with Standard_Mixed_Residuals;
with DoblDobl_Mixed_Residuals;
with QuadDobl_Mixed_Residuals;

package body Homotopy_Mixed_Residuals is

  function Standard_AbsVal_Homotopy
             return Standard_Complex_Poly_SysFun.Eval_Poly_Sys is

    hom : constant Standard_Complex_Poly_Systems.Poly_Sys
        := Standard_Homotopy.Homotopy_System;
    neq : constant integer32 := hom'last;
    abh : Standard_Complex_Poly_Systems.Poly_Sys(1..neq)
        := Standard_Mixed_Residuals.AbsVal(hom);
    res : constant Standard_Complex_Poly_SysFun.Eval_Poly_Sys(1..neq)
        := Standard_Complex_Poly_SysFun.Create(abh);

  begin
    Standard_Complex_Poly_Systems.Clear(abh);
    return res;
  end Standard_AbsVal_Homotopy;

  function DoblDobl_AbsVal_Homotopy
             return DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys is

    hom : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
        := DoblDobl_Homotopy.Homotopy_System;
    neq : constant integer32 := hom'last;
    abh : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..neq)
        := DoblDobl_Mixed_Residuals.AbsVal(hom);
    res : constant DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(1..neq)
        := DoblDobl_Complex_Poly_SysFun.Create(abh);

  begin
    DoblDobl_Complex_Poly_Systems.Clear(abh);
    return res;
  end DoblDobl_AbsVal_Homotopy;

  function QuadDobl_AbsVal_Homotopy
             return QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys is

    hom : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
        := QuadDobl_Homotopy.Homotopy_System;
    neq : constant integer32 := hom'last;
    abh : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..neq)
        := QuadDobl_Mixed_Residuals.AbsVal(hom);
    res : constant QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys(1..neq)
        := QuadDobl_Complex_Poly_SysFun.Create(abh);

  begin
    QuadDobl_Complex_Poly_Systems.Clear(abh);
    return res;
  end QuadDobl_AbsVal_Homotopy;

  function Residual ( abh : Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : Standard_Complex_Vectors.Vector;
                      t : Standard_Complex_Numbers.Complex_Number )
                    return double_float is

    val : Standard_Complex_Vectors.Vector(abh'range);
    abz : constant Standard_Complex_Vectors.Vector(z'range)
        := Standard_Mixed_Residuals.AbsVal(z);
    azt : Standard_Complex_Vectors.Vector(z'first..z'last+1);
    avl : Standard_Complex_Vectors.Vector(abh'range);
    len : constant double_float := double_float(val'last);
    vpz,vap : double_float;
    res : double_float := 0.0;

  begin
    if Standard_Coefficient_Homotopy.Number_of_Equations = -1
     then val := Standard_Homotopy.Eval(z,t);
     else val := Standard_Coefficient_Homotopy.Eval(z,t);
    end if;
    azt(abz'range) := abz;
    azt(azt'last) := t;
    avl := Standard_Complex_Poly_SysFun.Eval(abh,azt);
    for i in val'range loop
      vpz := Standard_Complex_Numbers_Polar.Radius(val(i));
      vap := Standard_Complex_Numbers_Polar.Radius(avl(i));
      res := res + vpz/(vap+1.0);
    end loop;
    res := res/len;
    return res;
  end Residual;

  function Residual ( abh : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : DoblDobl_Complex_Vectors.Vector;
                      t : DoblDobl_Complex_Numbers.Complex_Number )
                    return double_double is

    val : DoblDobl_Complex_Vectors.Vector(abh'range);
    abz : constant DoblDobl_Complex_Vectors.Vector(z'range)
        := DoblDobl_Mixed_Residuals.AbsVal(z);
    azt : DoblDobl_Complex_Vectors.Vector(z'first..z'last+1);
    avl : DoblDobl_Complex_Vectors.Vector(abh'range);
    len : constant double_float := double_float(val'last);
    vpz,vap : double_double;
    res : double_double := Double_Double_Numbers.create(0.0);

  begin
    if DoblDobl_Coefficient_Homotopy.Number_of_Equations = -1
     then val := DoblDobl_Homotopy.Eval(z,t);
     else val := DoblDobl_Coefficient_Homotopy.Eval(z,t);
    end if;
    azt(abz'range) := abz;
    azt(azt'last) := t;
    avl := DoblDobl_Complex_Poly_SysFun.Eval(abh,azt);
    for i in val'range loop
      vpz := DoblDobl_Complex_Numbers_Polar.Radius(val(i));
      vap := DoblDobl_Complex_Numbers_Polar.Radius(avl(i));
      res := res + vpz/(vap+1.0);
    end loop;
    res := res/len;
    return res;
  end Residual;

  function Residual ( abh : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : QuadDobl_Complex_Vectors.Vector;
                      t : QuadDobl_Complex_Numbers.Complex_Number )
                    return quad_double is

    val : QuadDobl_Complex_Vectors.Vector(abh'range);
    abz : constant QuadDobl_Complex_Vectors.Vector(z'range)
        := QuadDobl_Mixed_Residuals.AbsVal(z);
    azt : QuadDobl_Complex_Vectors.Vector(z'first..z'last+1);
    avl : QuadDobl_Complex_Vectors.Vector(abh'range);
    len : constant double_float := double_float(val'last);
    vpz,vap : quad_double;
    res : quad_double := Quad_Double_Numbers.create(0.0);

  begin
    if QuadDobl_Coefficient_Homotopy.Number_of_Equations = -1
     then val := QuadDobl_Homotopy.Eval(z,t);
     else val := QuadDobl_Coefficient_Homotopy.Eval(z,t);
    end if;
    azt(abz'range) := abz;
    azt(azt'last) := t;
    avl := QuadDobl_Complex_Poly_SysFun.Eval(abh,azt);
    for i in val'range loop
      vpz := QuadDobl_Complex_Numbers_Polar.Radius(val(i));
      vap := QuadDobl_Complex_Numbers_Polar.Radius(avl(i));
      res := res + vpz/(vap+1.0);
    end loop;
    res := res/len;
    return res;
  end Residual;

  function Residual ( file : file_type;
                      abh : Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : Standard_Complex_Vectors.Vector;
                      t : Standard_Complex_Numbers.Complex_Number )
                    return double_float is

    val : Standard_Complex_Vectors.Vector(abh'range);
    abz : constant Standard_Complex_Vectors.Vector(z'range)
        := Standard_Mixed_Residuals.AbsVal(z);
    azt : Standard_Complex_Vectors.Vector(z'first..z'last+1);
    avl : Standard_Complex_Vectors.Vector(abh'range);
    len : constant double_float := double_float(val'last);
    vaz,vpz,vap,res : double_float;
    sumres : double_float := 0.0;

  begin
    if Standard_Coefficient_Homotopy.Number_of_Equations = -1
     then val := Standard_Homotopy.Eval(z,t);
     else val := Standard_Coefficient_Homotopy.Eval(z,t);
    end if;
    vaz := Standard_Complex_Vector_Norms.Max_Norm(abz);
    azt(abz'range) := abz;
    azt(azt'last) := t;
    avl := Standard_Complex_Poly_SysFun.Eval(abh,azt);
    for i in val'range loop
      vpz := Standard_Complex_Numbers_Polar.Radius(val(i));
      vap := Standard_Complex_Numbers_Polar.Radius(avl(i));
      res := vpz/(vap+1.0);
      Standard_Mixed_Residuals.Write_Residuals(file,vaz,vpz,vap,res);
      sumres := sumres + res;
    end loop;
    sumres := sumres/len;
    return sumres;
  end Residual;

  function Residual ( file : file_type;
                      abh : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : DoblDobl_Complex_Vectors.Vector;
                      t : DoblDobl_Complex_Numbers.Complex_Number )
                    return double_double is

    val : DoblDobl_Complex_Vectors.Vector(abh'range);
    abz : constant DoblDobl_Complex_Vectors.Vector(z'range)
        := DoblDobl_Mixed_Residuals.AbsVal(z);
    azt : DoblDobl_Complex_Vectors.Vector(z'first..z'last+1);
    avl : DoblDobl_Complex_Vectors.Vector(abh'range);
    len : constant double_float := double_float(val'last);
    vaz,vpz,vap,res : double_double;
    sumres : double_double := Double_Double_Numbers.create(0.0);

  begin
    if DoblDobl_Coefficient_Homotopy.Number_of_Equations = -1
     then val := DoblDobl_Homotopy.Eval(z,t);
     else val := DoblDobl_Coefficient_Homotopy.Eval(z,t);
    end if;
    vaz := DoblDobl_Complex_Vector_Norms.Max_Norm(abz);
    azt(abz'range) := abz;
    azt(azt'last) := t;
    avl := DoblDobl_Complex_Poly_SysFun.Eval(abh,azt);
    for i in val'range loop
      vpz := DoblDobl_Complex_Numbers_Polar.Radius(val(i));
      vap := DoblDobl_Complex_Numbers_Polar.Radius(avl(i));
      res := vpz/(vap + 1.0);
      DoblDobl_Mixed_Residuals.Write_Residuals(file,vaz,vpz,vap,res);
      sumres := sumres + res;
    end loop;
    sumres := sumres/len;
    return sumres;
  end Residual;

  function Residual ( file : file_type;
                      abh : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : QuadDobl_Complex_Vectors.Vector;
                      t : QuadDobl_Complex_Numbers.Complex_Number )
                    return quad_double is

    val : QuadDobl_Complex_Vectors.Vector(abh'range);
    abz : constant QuadDobl_Complex_Vectors.Vector(z'range)
        := QuadDobl_Mixed_Residuals.AbsVal(z);
    azt : QuadDobl_Complex_Vectors.Vector(z'first..z'last+1);
    avl : QuadDobl_Complex_Vectors.Vector(abh'range);
    len : constant double_float := double_float(val'last);
    vaz,vpz,vap,res : quad_double;
    sumres : quad_double := Quad_Double_Numbers.create(0.0);

  begin
    if QuadDobl_Coefficient_Homotopy.Number_of_Equations = -1
     then val := QuadDobl_Homotopy.Eval(z,t);
     else val := QuadDobl_Coefficient_Homotopy.Eval(z,t);
    end if;
    vaz := QuadDobl_Complex_Vector_Norms.Max_Norm(abz);
    azt(abz'range) := abz;
    azt(azt'last) := t;
    avl := QuadDobl_Complex_Poly_SysFun.Eval(abh,azt);
    for i in val'range loop
      vpz := QuadDobl_Complex_Numbers_Polar.Radius(val(i));
      vap := QuadDobl_Complex_Numbers_Polar.Radius(avl(i));
      res := vpz/(vap+1.0);
      QuadDobl_Mixed_Residuals.Write_Residuals(file,vaz,vpz,vap,res);
      sumres := sumres + res;
    end loop;
    sumres := sumres/len;
    return sumres;
  end Residual;

end Homotopy_Mixed_Residuals;
