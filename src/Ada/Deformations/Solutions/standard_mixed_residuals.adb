with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_Polar;     use Standard_Complex_Numbers_Polar;
with Standard_Complex_Poly_Functions;

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

  function AbsVal ( t : Term ) return Term is

    res : Term;
    rad : double_float;

  begin
    Copy(t,res);
    rad := Radius(t.cf);
    res.cf := Create(rad);
    return res;
  end AbsVal;

  function AbsVal ( p : Poly ) return Poly is

    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Term; cont : out boolean ) is 

      abt : Term := AbsVal(t);

    begin
      Add(res,abt);
      cont := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end AbsVal;

  function Residual ( pol,abp : Poly; z : Vector ) return double_float is

    val : constant Complex_Number
        := Standard_Complex_Poly_Functions.Eval(pol,z);
    abz : constant Vector(z'range) := AbsVal(z);
    avl : constant Complex_Number
        := Standard_Complex_Poly_Functions.Eval(abp,abz);
    res : constant double_float := Radius(val)/(Radius(avl) + 1.0);

  begin
    return res;
  end Residual;

end Standard_Mixed_Residuals;
