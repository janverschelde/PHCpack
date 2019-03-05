with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_Polar;     use Standard_Complex_Numbers_Polar;
with Standard_Complex_Poly_Functions;
with Standard_Complex_Laur_Functions;

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

end Standard_Mixed_Residuals;
