with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers_Polar;     use QuadDobl_Complex_Numbers_Polar;
with QuadDobl_Complex_Vector_Norms;

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

  function Residual ( val,avl : Complex_Number ) return quad_double is

    vpz : constant quad_double := Radius(val);
    vap : constant quad_double := Radius(avl);
    res : constant quad_double := vpz/(vap+1.0);

  begin
    return res;
  end Residual;

  procedure Residual ( val,avl : in Complex_Number;
                       vpz,vap,res : out quad_double ) is
  begin
    vpz := Radius(val);
    vap := Radius(avl);
    res := vpz/(vap+1.0);
  end Residual;

  function Residual ( val,avl : Vector ) return quad_double is

    res : quad_double := create(0.0);
    len : constant double_float := double_float(val'length);

  begin
    for i in val'range loop
      res := res + Residual(val(i),avl(i));
    end loop;
    return (res/len);
  end Residual;

  procedure Residual ( pol,abp : in QuadDobl_Complex_Polynomials.Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out quad_double ) is

    val : constant Complex_Number
        := QuadDobl_Complex_Poly_Functions.Eval(pol,z);
    avl : Complex_Number;

  begin
    abz := AbsVal(z);
    avl := QuadDobl_Complex_Poly_Functions.Eval(abp,abz);
    vaz := QuadDobl_Complex_Vector_Norms.Max_Norm(abz);
    Residual(val,avl,vpz,vap,res);
  end Residual;

  procedure Residual ( pol,abp : in QuadDobl_Complex_Poly_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out quad_double ) is

    val : constant Complex_Number
        := QuadDobl_Complex_Poly_Functions.Eval(pol,z);
    avl : Complex_Number;

  begin
    abz := AbsVal(z);
    avl := QuadDobl_Complex_Poly_Functions.Eval(abp,abz);
    vaz := QuadDobl_Complex_Vector_Norms.Max_Norm(abz);
    Residual(val,avl,vpz,vap,res);
  end Residual;

  procedure Residual ( pol,abp : in QuadDobl_Complex_Laurentials.Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out quad_double ) is

    val : constant Complex_Number
        := QuadDobl_Complex_Laur_Functions.Eval(pol,z);
    avl : Complex_Number;

  begin
    abz := AbsVal(z);
    avl := QuadDobl_Complex_Laur_Functions.Eval(abp,abz);
    vaz := QuadDobl_Complex_Vector_Norms.Max_Norm(abz);
    Residual(val,avl,vpz,vap,res);
  end Residual;

  procedure Residual ( pol,abp : in QuadDobl_Complex_Laur_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out quad_double ) is

    val : constant Complex_Number
        := QuadDobl_Complex_Laur_Functions.Eval(pol,z);
    avl : Complex_Number;

  begin
    abz := AbsVal(z);
    avl := QuadDobl_Complex_Laur_Functions.Eval(abp,abz);
    vaz := QuadDobl_Complex_Vector_Norms.Max_Norm(abz);
    Residual(val,avl,vpz,vap,res);
  end Residual;

  procedure Write_Residuals
              ( file : in file_type;
                vaz,vpz,vap,res : in quad_double ) is
  begin
    put(file,"  vaz : "); put(file,vaz,3);
    put(file,"  vpz : "); put(file,vpz,3);
    put(file,"  vap : "); put(file,vap,3);
    put(file,"  res : "); put(file,res,3); new_line(file);
  end Write_Residuals;

  procedure Residual ( file : in file_type;
                       pol,abp : in QuadDobl_Complex_Polynomials.Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out quad_double ) is
  begin
    Residual(pol,abp,z,abz,vaz,vpz,vap,res);
    Write_Residuals(file,vaz,vpz,vap,res);
  end Residual;

  procedure Residual ( file : in file_type;
                       pol,abp : in QuadDobl_Complex_Poly_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out quad_double ) is
  begin
    Residual(pol,abp,z,abz,vaz,vpz,vap,res);
    Write_Residuals(file,vaz,vpz,vap,res);
  end Residual;

  procedure Residual ( file : in file_type;
                       pol,abp : in QuadDobl_Complex_Laurentials.Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out quad_double ) is
  begin
    Residual(pol,abp,z,abz,vaz,vpz,vap,res);
    Write_Residuals(file,vaz,vpz,vap,res);
  end Residual;

  procedure Residual ( file : in file_type;
                       pol,abp : in QuadDobl_Complex_Laur_Functions.Eval_Poly;
                       z : in Vector; abz : out Vector;
                       vaz,vpz,vap,res : out quad_double ) is
  begin
    Residual(pol,abp,z,abz,vaz,vpz,vap,res);
    Write_Residuals(file,vaz,vpz,vap,res);
  end Residual;

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

  function Residual ( file : file_type;
                      pol,abp : QuadDobl_Complex_Polynomials.Poly;
                      z : Vector ) return quad_double is

    abz : Vector(z'range);
    res,vaz,vpz,vap : quad_double;

  begin
    Residual(file,pol,abp,z,abz,vaz,vpz,vap,res);
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

  function Residual ( file : file_type;
                      pol,abp : QuadDobl_Complex_Poly_Functions.Eval_Poly;
                      z : Vector ) return quad_double is

    abz : Vector(z'range);
    res,vaz,vpz,vap : quad_double;

  begin
    Residual(file,pol,abp,z,abz,vaz,vpz,vap,res);
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

  function Residual ( file : file_type;
                      pol,abp : QuadDobl_Complex_Laurentials.Poly;
                      z : Vector ) return quad_double is

    abz : Vector(z'range);
    res,vaz,vpz,vap : quad_double;

  begin
    Residual(file,pol,abp,z,abz,vaz,vpz,vap,res);
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

  function Residual ( file : file_type;
                      pol,abp : QuadDobl_Complex_Laur_Functions.Eval_Poly;
                      z : Vector ) return quad_double is

    abz : Vector(z'range);
    res,vaz,vpz,vap : quad_double;

  begin
    Residual(file,pol,abp,z,abz,vaz,vpz,vap,res);
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

  function Residual ( file : file_type;
                      pol,abp : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                      z : Vector ) return quad_double is

    len : constant double_float := double_float(pol'last);
    res : quad_double := create(0.0);

  begin
    for i in pol'range loop
      res := res + Residual(file,pol(i),abp(i),z);
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

  function Residual ( file : file_type;
                      pol,abp : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                      z : Vector ) return quad_double is

    len : constant double_float := double_float(pol'last);
    res : quad_double := create(0.0);

  begin
    for i in pol'range loop
      res := res + Residual(file,pol(i),abp(i),z);
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

  function Residual ( file : file_type;
                      pol,abp : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                      z : Vector ) return quad_double is

    len : constant double_float := double_float(pol'last);
    res : quad_double := create(0.0);

  begin
    for i in pol'range loop
      res := res + Residual(file,pol(i),abp(i),z);
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

  function Residual ( file : file_type;
                      pol,abp : QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                      z : Vector ) return quad_double is

    len : constant double_float := double_float(pol'last);
    res : quad_double := create(0.0);

  begin
    for i in pol'range loop
      res := res + Residual(file,pol(i),abp(i),z);
    end loop;
    res := res/len;
    return res;
  end Residual;

  function Mixed_Residual
              ( valres,absres : in QuadDobl_Complex_Vectors.Vector )
              return quad_double is

    res : quad_double := Create(0.0);
    len : constant quad_double := create(integer(valres'last));

  begin
    for k in valres'range loop 
      res := res + Radius(valres(k))/(Radius(absres(k)) + 1.0);
    end loop;
    return (res/len);
  end Mixed_Residual;

end QuadDobl_Mixed_Residuals;
