with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_Polar;     use Standard_Complex_Numbers_Polar;
with Standard_Complex_Vector_Norms;

package body Standard_Mixed_Residuals is

  function AbsVal ( z : Vector ) return Vector is

    res : Vector(z'range);
    rad : double_float;

  begin
    for i in z'range loop
      rad := Radius(z(i));
      res(i) := Create(rad);
    end loop;
    return res;
  end AbsVal;

  function AbsVal ( t : Standard_Complex_Polynomials.Term )
                  return Standard_Complex_Polynomials.Term is

    res : Standard_Complex_Polynomials.Term;
    rad : double_float;

  begin
    Standard_Complex_Polynomials.Copy(t,res);
    rad := Radius(t.cf);
    res.cf := Create(rad);
    return res;
  end AbsVal;

  function AbsVal ( t : Standard_Complex_Laurentials.Term )
                  return Standard_Complex_Laurentials.Term is

    res : Standard_Complex_Laurentials.Term;
    rad : double_float;

  begin
    Standard_Complex_Laurentials.Copy(t,res);
    rad := Radius(t.cf);
    res.cf := Create(rad);
    return res;
  end AbsVal;

  function AbsVal ( p : Standard_Complex_Polynomials.Poly )
                  return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Term; cont : out boolean ) is 

      abt : Term := AbsVal(t);

    begin
      Add(res,abt);
      Standard_Complex_Polynomials.Clear(abt);
      cont := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end AbsVal;

  function AbsVal ( p : Standard_Complex_Laurentials.Poly )
                  return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;

    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Term; cont : out boolean ) is 

      abt : Term := AbsVal(t);

    begin
      Add(res,abt);
      Standard_Complex_Laurentials.Clear(abt);
      cont := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end AbsVal;

  function AbsVal ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                  return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := AbsVal(p(i));
    end loop;
    return res;
  end AbsVal;

  function AbsVal ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                  return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := AbsVal(p(i));
    end loop;
    return res;
  end AbsVal;

  function Residual ( val,avl : Complex_Number ) return double_float is

    vpz : constant double_float := Radius(val);
    vap : constant double_float := Radius(avl);
    res : constant double_float := vpz/(vap+1.0);

  begin
    return res;
  end Residual;

  procedure Residual ( val,avl : in Complex_Number;
                       vpz,vap,res : out double_float ) is
  begin
    vpz := Radius(val);
    vap := Radius(avl);
    res := vpz/(vap+1.0);
  end Residual;

  function Residual ( val,avl : Vector ) return double_float is

    res : double_float := 0.0;
    len : constant double_float := double_float(val'length);

  begin
    for i in val'range loop
      res := res + Residual(val(i),avl(i));
    end loop;
    return (res/len);
  end Residual;

  procedure Residual ( pol,abp : in Standard_Complex_Polynomials.Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out double_float ) is

    val : constant Complex_Number
        := Standard_Complex_Poly_Functions.Eval(pol,z);
    avl : Complex_Number;

  begin
    abz := AbsVal(z);
    avl := Standard_Complex_Poly_Functions.Eval(abp,abz);
    vaz := Standard_Complex_Vector_Norms.Max_Norm(abz);
    Residual(val,avl,vpz,vap,res);
  end Residual;

  procedure Residual ( pol,abp : in Standard_Complex_Poly_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out double_float ) is

    val : constant Complex_Number
        := Standard_Complex_Poly_Functions.Eval(pol,z);
    avl : Complex_Number;

  begin
    abz := AbsVal(z);
    avl := Standard_Complex_Poly_Functions.Eval(abp,abz);
    vaz := Standard_Complex_Vector_Norms.Max_Norm(abz);
    Residual(val,avl,vpz,vap,res);
  end Residual;

  procedure Residual ( pol,abp : in Standard_Complex_Laurentials.Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out double_float ) is

    val : constant Complex_Number
        := Standard_Complex_Laur_Functions.Eval(pol,z);
    avl : Complex_Number;

  begin
    abz := AbsVal(z);
    avl := Standard_Complex_Laur_Functions.Eval(abp,abz);
    vaz := Standard_Complex_Vector_Norms.Max_Norm(abz);
    Residual(val,avl,vpz,vap,res);
  end Residual;

  procedure Residual ( pol,abp : in Standard_Complex_Laur_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out double_float ) is

    val : constant Complex_Number
        := Standard_Complex_Laur_Functions.Eval(pol,z);
    avl : Complex_Number;

  begin
    abz := AbsVal(z);
    avl := Standard_Complex_Laur_Functions.Eval(abp,abz);
    vaz := Standard_Complex_Vector_Norms.Max_Norm(abz);
    Residual(val,avl,vpz,vap,res);
  end Residual;

  procedure Write_Residuals
              ( file : in file_type;
                vaz,vpz,vap,res : in double_float ) is
  begin
    put(file,"  vaz : "); put(file,vaz,3);
    put(file,"  vpz : "); put(file,vpz,3);
    put(file,"  vap : "); put(file,vap,3);
    put(file,"  res : "); put(file,res,3); new_line(file);
  end Write_Residuals;

  procedure Residual ( file : in file_type;
                       pol,abp : in Standard_Complex_Polynomials.Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out double_float ) is
  begin
    Residual(pol,abp,z,abz,vaz,vpz,vap,res);
    Write_Residuals(file,vaz,vpz,vap,res);
  end Residual;

  procedure Residual ( file : in file_type;
                       pol,abp : in Standard_Complex_Poly_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out double_float ) is
  begin
    Residual(pol,abp,z,abz,vaz,vpz,vap,res);
    Write_Residuals(file,vaz,vpz,vap,res);
  end Residual;

  procedure Residual ( file : in file_type;
                       pol,abp : in Standard_Complex_Laurentials.Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out double_float ) is
  begin
    Residual(pol,abp,z,abz,vaz,vpz,vap,res);
    Write_Residuals(file,vaz,vpz,vap,res);
  end Residual;

  procedure Residual ( file : in file_type;
                       pol,abp : in Standard_Complex_Laur_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out double_float ) is
  begin
    Residual(pol,abp,z,abz,vaz,vpz,vap,res);
    Write_Residuals(file,vaz,vpz,vap,res);
  end Residual;

  function Residual ( pol,abp : Standard_Complex_Polynomials.Poly;
                      z : Vector ) return double_float is

    val : constant Complex_Number
        := Standard_Complex_Poly_Functions.Eval(pol,z);
    abz : constant Vector(z'range) := AbsVal(z);
    avl : constant Complex_Number
        := Standard_Complex_Poly_Functions.Eval(abp,abz);
    res : constant double_float := Radius(val)/(Radius(avl) + 1.0);

  begin
    return res;
  end Residual;

  function Residual ( file : file_type;
                      pol,abp : Standard_Complex_Polynomials.Poly;
                      z : Vector ) return double_float is

    abz : Vector(z'range);
    res,vaz,vpz,vap : double_float;

  begin
    Residual(file,pol,abp,z,abz,vaz,vpz,vap,res);
    return res;
  end Residual;

  function Residual ( pol,abp : Standard_Complex_Poly_Functions.Eval_Poly;
                      z : Vector ) return double_float is

    val : constant Complex_Number
        := Standard_Complex_Poly_Functions.Eval(pol,z);
    abz : constant Vector(z'range) := AbsVal(z);
    avl : constant Complex_Number
        := Standard_Complex_Poly_Functions.Eval(abp,abz);
    res : constant double_float := Radius(val)/(Radius(avl) + 1.0);

  begin
    return res;
  end Residual;

  function Residual ( file : file_type;
                      pol,abp : Standard_Complex_Poly_Functions.Eval_Poly;
                      z : Vector ) return double_float is

    abz : Vector(z'range);
    res,vaz,vpz,vap : double_float;

  begin
    Residual(file,pol,abp,z,abz,vaz,vpz,vap,res);
    return res;
  end Residual;

  function Residual ( pol,abp : Standard_Complex_Laurentials.Poly;
                      z : Vector ) return double_float is

    val : constant Complex_Number
        := Standard_Complex_Laur_Functions.Eval(pol,z);
    abz : constant Vector(z'range) := AbsVal(z);
    avl : constant Complex_Number
        := Standard_Complex_Laur_Functions.Eval(abp,abz);
    res : constant double_float := Radius(val)/(Radius(avl) + 1.0);

  begin
    return res;
  end Residual;

  function Residual ( file : file_type;
                      pol,abp : Standard_Complex_Laurentials.Poly;
                      z : Vector ) return double_float is

    abz : Vector(z'range);
    res,vaz,vpz,vap : double_float;

  begin
    Residual(file,pol,abp,z,abz,vaz,vpz,vap,res);
    return res;
  end Residual;

  function Residual ( pol,abp : Standard_Complex_Laur_Functions.Eval_Poly;
                      z : Vector ) return double_float is

    val : constant Complex_Number
        := Standard_Complex_Laur_Functions.Eval(pol,z);
    abz : constant Vector(z'range) := AbsVal(z);
    avl : constant Complex_Number
        := Standard_Complex_Laur_Functions.Eval(abp,abz);
    res : constant double_float := Radius(val)/(Radius(avl) + 1.0);

  begin
    return res;
  end Residual;

  function Residual ( file : file_type;
                      pol,abp : Standard_Complex_Laur_Functions.Eval_Poly;
                      z : Vector ) return double_float is

    abz : Vector(z'range);
    res,vaz,vpz,vap : double_float;

  begin
    Residual(file,pol,abp,z,abz,vaz,vpz,vap,res);
    return res;
  end Residual;

  function Residual ( pol,abp : Standard_Complex_Poly_Systems.Poly_Sys;
                      z : Vector ) return double_float is

    len : constant double_float := double_float(pol'last);
    res : double_float := 0.0;

  begin
    for i in pol'range loop
      res := res + Residual(pol(i),abp(i),z);
    end loop;
    res := res/len;
    return res;
  end Residual;

  function Residual ( file : file_type;
                      pol,abp : Standard_Complex_Poly_Systems.Poly_Sys;
                      z : Vector ) return double_float is

    len : constant double_float := double_float(pol'last);
    res : double_float := 0.0;

  begin
    for i in pol'range loop
      res := res + Residual(file,pol(i),abp(i),z);
    end loop;
    res := res/len;
    return res;
  end Residual;

  function Residual ( pol,abp : Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : Vector ) return double_float is

    len : constant double_float := double_float(pol'last);
    res : double_float := 0.0;

  begin
    for i in pol'range loop
      res := res + Residual(pol(i),abp(i),z);
    end loop;
    res := res/len;
    return res;
  end Residual;

  function Residual ( file : file_type;
                      pol,abp : Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : Vector ) return double_float is

    len : constant double_float := double_float(pol'last);
    res : double_float := 0.0;

  begin
    for i in pol'range loop
      res := res + Residual(file,pol(i),abp(i),z);
    end loop;
    res := res/len;
    return res;
  end Residual;

  function Residual ( pol,abp : Standard_Complex_Laur_Systems.Laur_Sys;
                      z : Vector ) return double_float is

    len : constant double_float := double_float(pol'last);
    res : double_float := 0.0;

  begin
    for i in pol'range loop
      res := res + Residual(pol(i),abp(i),z);
    end loop;
    res := res/len;
    return res;
  end Residual;

  function Residual ( file : file_type;
                      pol,abp : Standard_Complex_Laur_Systems.Laur_Sys;
                      z : Vector ) return double_float is

    len : constant double_float := double_float(pol'last);
    res : double_float := 0.0;

  begin
    for i in pol'range loop
      res := res + Residual(file,pol(i),abp(i),z);
    end loop;
    res := res/len;
    return res;
  end Residual;

  function Residual ( pol,abp : Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                      z : Vector ) return double_float is

    len : constant double_float := double_float(pol'last);
    res : double_float := 0.0;

  begin
    for i in pol'range loop
      res := res + Residual(pol(i),abp(i),z);
    end loop;
    res := res/len;
    return res;
  end Residual;

  function Residual ( file : file_type;
                      pol,abp : Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                      z : Vector ) return double_float is

    len : constant double_float := double_float(pol'last);
    res : double_float := 0.0;

  begin
    for i in pol'range loop
      res := res + Residual(file,pol(i),abp(i),z);
    end loop;
    res := res/len;
    return res;
  end Residual;

  function Mixed_Residual
              ( valres,absres : in Standard_Complex_Vectors.Vector )
              return double_float is

    res : double_float := 0.0;
    len : constant double_float := double_float(valres'last);

  begin
    for k in valres'range loop 
      res := res + Radius(valres(k))/(Radius(absres(k)) + 1.0);
    end loop;
    return (res/len);
  end Mixed_Residual;

end Standard_Mixed_Residuals;
