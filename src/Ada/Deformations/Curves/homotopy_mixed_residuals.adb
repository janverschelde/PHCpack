with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Complex_Numbers_Polar;
with DoblDobl_Complex_Numbers_Polar;
with QuadDobl_Complex_Numbers_Polar;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Standard_Mixed_Residuals;
with DoblDobl_Mixed_Residuals;
with QuadDobl_Mixed_Residuals;

package body Homotopy_Mixed_Residuals is

  function Residual ( abh : Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : Standard_Complex_Vectors.Vector;
                      t : Standard_Complex_Numbers.Complex_Number )
                    return double_float is

    val : constant Standard_Complex_Vectors.Vector(abh'range)
        := Standard_Homotopy.Eval(z,t);
    abz : constant Standard_Complex_Vectors.Vector(z'range)
        := Standard_Mixed_Residuals.AbsVal(z);
    azt : Standard_Complex_Vectors.Vector(z'first..z'last+1);
    avl : Standard_Complex_Vectors.Vector(abh'range);
    len : constant double_float := double_float(val'last);
    vpz,vap : double_float;
    res : double_float := 0.0;

  begin
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

    val : constant DoblDobl_Complex_Vectors.Vector(abh'range)
        := DoblDobl_Homotopy.Eval(z,t);
    abz : constant DoblDobl_Complex_Vectors.Vector(z'range)
        := DoblDobl_Mixed_Residuals.AbsVal(z);
    azt : DoblDobl_Complex_Vectors.Vector(z'first..z'last+1);
    avl : DoblDobl_Complex_Vectors.Vector(abh'range);
    len : constant double_float := double_float(val'last);
    vpz,vap : double_double;
    res : double_double := Double_Double_Numbers.create(0.0);

  begin
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

    val : constant QuadDobl_Complex_Vectors.Vector(abh'range)
        := QuadDobl_Homotopy.Eval(z,t);
    abz : constant QuadDobl_Complex_Vectors.Vector(z'range)
        := QuadDobl_Mixed_Residuals.AbsVal(z);
    azt : QuadDobl_Complex_Vectors.Vector(z'first..z'last+1);
    avl : QuadDobl_Complex_Vectors.Vector(abh'range);
    len : constant double_float := double_float(val'last);
    vpz,vap : quad_double;
    res : quad_double := Quad_Double_Numbers.create(0.0);

  begin
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

    val : constant Standard_Complex_Vectors.Vector(abh'range)
        := Standard_Homotopy.Eval(z,t);
    abz : constant Standard_Complex_Vectors.Vector(z'range)
        := Standard_Mixed_Residuals.AbsVal(z);
    azt : Standard_Complex_Vectors.Vector(z'first..z'last+1);
    avl : Standard_Complex_Vectors.Vector(abh'range);
    len : constant double_float := double_float(val'last);
    vaz,vpz,vap,res : double_float;
    sumres : double_float := 0.0;

  begin
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

    val : constant DoblDobl_Complex_Vectors.Vector(abh'range)
        := DoblDobl_Homotopy.Eval(z,t);
    abz : constant DoblDobl_Complex_Vectors.Vector(z'range)
        := DoblDobl_Mixed_Residuals.AbsVal(z);
    azt : DoblDobl_Complex_Vectors.Vector(z'first..z'last+1);
    avl : DoblDobl_Complex_Vectors.Vector(abh'range);
    len : constant double_float := double_float(val'last);
    vaz,vpz,vap,res : double_double;
    sumres : double_double := Double_Double_Numbers.create(0.0);

  begin
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

    val : constant QuadDobl_Complex_Vectors.Vector(abh'range)
        := QuadDobl_Homotopy.Eval(z,t);
    abz : constant QuadDobl_Complex_Vectors.Vector(z'range)
        := QuadDobl_Mixed_Residuals.AbsVal(z);
    azt : QuadDobl_Complex_Vectors.Vector(z'first..z'last+1);
    avl : QuadDobl_Complex_Vectors.Vector(abh'range);
    len : constant double_float := double_float(val'last);
    vaz,vpz,vap,res : quad_double;
    sumres : quad_double := Quad_Double_Numbers.create(0.0);

  begin
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
