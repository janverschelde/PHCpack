with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_Polar;     use QuadDobl_Complex_Numbers_Polar;

package body QuadDobl_Mixed_Residuals is

  function AbsVal ( z : Vector ) return Vector is

    res : Vector(z'range);
    rad : quad_double;

  begin
    for i in z'range loop
      rad := Radius(z(i));
      res(i) := Create(rad);
    end loop;
    return res;
  end AbsVal;

  function AbsVal ( t : QuadDobl_Complex_Polynomials.Term )
                  return QuadDobl_Complex_Polynomials.Term is

    res : QuadDobl_Complex_Polynomials.Term;
    rad : quad_double;

  begin
    QuadDobl_Complex_Polynomials.Copy(t,res);
    rad := Radius(t.cf);
    res.cf := Create(rad);
    return res;
  end AbsVal;

  function AbsVal ( t : QuadDobl_Complex_Laurentials.Term )
                  return QuadDobl_Complex_Laurentials.Term is

    res : QuadDobl_Complex_Laurentials.Term;
    rad : quad_double;

  begin
    QuadDobl_Complex_Laurentials.Copy(t,res);
    rad := Radius(t.cf);
    res.cf := Create(rad);
    return res;
  end AbsVal;

  function AbsVal ( p : QuadDobl_Complex_Polynomials.Poly )
                  return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;

    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Term; cont : out boolean ) is 

      abt : Term := AbsVal(t);

    begin
      Add(res,abt);
      QuadDobl_Complex_Polynomials.Clear(abt);
      cont := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end AbsVal;

  function AbsVal ( p : QuadDobl_Complex_Laurentials.Poly )
                  return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Laurentials;

    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Term; cont : out boolean ) is 

      abt : Term := AbsVal(t);

    begin
      Add(res,abt);
      QuadDobl_Complex_Laurentials.Clear(abt);
      cont := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end AbsVal;

  function AbsVal ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                  return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := AbsVal(p(i));
    end loop;
    return res;
  end AbsVal;

  function AbsVal ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                  return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := AbsVal(p(i));
    end loop;
    return res;
  end AbsVal;

  function Residual ( pol,abp : QuadDobl_Complex_Polynomials.Poly;
                      z : Vector ) return quad_double is

    val : constant Complex_Number
        := QuadDobl_Complex_Poly_Functions.Eval(pol,z);
    abz : constant Vector(z'range) := AbsVal(z);
    avl : constant Complex_Number
        := QuadDobl_Complex_Poly_Functions.Eval(abp,abz);
    res : constant quad_double := Radius(val)/(Radius(avl) + 1.0);

  begin
    return res;
  end Residual;

  function Residual ( pol,abp : QuadDobl_Complex_Poly_Functions.Eval_Poly;
                      z : Vector ) return quad_double is

    val : constant Complex_Number
        := QuadDobl_Complex_Poly_Functions.Eval(pol,z);
    abz : constant Vector(z'range) := AbsVal(z);
    avl : constant Complex_Number
        := QuadDobl_Complex_Poly_Functions.Eval(abp,abz);
    res : constant quad_double := Radius(val)/(Radius(avl) + 1.0);

  begin
    return res;
  end Residual;

  function Residual ( pol,abp : QuadDobl_Complex_Laurentials.Poly;
                      z : Vector ) return quad_double is

    val : constant Complex_Number
        := QuadDobl_Complex_Laur_Functions.Eval(pol,z);
    abz : constant Vector(z'range) := AbsVal(z);
    avl : constant Complex_Number
        := QuadDobl_Complex_Laur_Functions.Eval(abp,abz);
    res : constant quad_double := Radius(val)/(Radius(avl) + 1.0);

  begin
    return res;
  end Residual;

  function Residual ( pol,abp : QuadDobl_Complex_Laur_Functions.Eval_Poly;
                      z : Vector ) return quad_double is

    val : constant Complex_Number
        := QuadDobl_Complex_Laur_Functions.Eval(pol,z);
    abz : constant Vector(z'range) := AbsVal(z);
    avl : constant Complex_Number
        := QuadDobl_Complex_Laur_Functions.Eval(abp,abz);
    res : constant quad_double := Radius(val)/(Radius(avl) + 1.0);

  begin
    return res;
  end Residual;

  function Residual ( pol,abp : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                      z : Vector ) return quad_double is

    len : constant double_float := double_float(pol'last);
    res : quad_double := create(0.0);

  begin
    for i in pol'range loop
      res := res + Residual(pol(i),abp(i),z);
    end loop;
    res := res/len;
    return res;
  end Residual;

  function Residual ( pol,abp : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : Vector ) return quad_double is

    len : constant double_float := double_float(pol'last);
    res : quad_double := create(0.0);

  begin
    for i in pol'range loop
      res := res + Residual(pol(i),abp(i),z);
    end loop;
    res := res/len;
    return res;
  end Residual;

  function Residual ( pol,abp : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                      z : Vector ) return quad_double is

    len : constant double_float := double_float(pol'last);
    res : quad_double := create(0.0);

  begin
    for i in pol'range loop
      res := res + Residual(pol(i),abp(i),z);
    end loop;
    res := res/len;
    return res;
  end Residual;

  function Residual ( pol,abp : QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                      z : Vector ) return quad_double is

    len : constant double_float := double_float(pol'last);
    res : quad_double := create(0.0);

  begin
    for i in pol'range loop
      res := res + Residual(pol(i),abp(i),z);
    end loop;
    res := res/len;
    return res;
  end Residual;

end QuadDobl_Mixed_Residuals;
