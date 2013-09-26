with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;

package body Standard_Aitken_Extrapolation is

  function Conjugate_Product ( x,y : Vector ) return Complex_Number is

    res : Complex_Number := Create(0.0);

  begin
    for i in x'range loop
      res := res + x(i)*Conjugate(y(i));
    end loop;
    return res;
  end Conjugate_Product;

  function Norm2 ( x : Vector ) return Complex_Number is
  begin
    return Conjugate_Product(x,x);
  end Norm2;

  function Accelerate ( v0,v1,v2 : Vector ) return Vector is

    res,d1,d2,dd,decrement : Vector(v0'range);
    numerator,denominator,s : Complex_Number;

  begin
    d2 := v2 - v1;
    d1 := v1 - v0;
    dd := d2 - d1;
    numerator := Conjugate_Product(d1,d1);
    denominator := Conjugate_Product(d1,dd);
    if ((REAL_PART(denominator) = 0.0) and (IMAG_PART(denominator) = 0.0))
     then res := v1;
     else s := numerator/denominator;
          decrement := s*d2;
          res := v1 - decrement;
    end if;
    return res;
  end Accelerate;

  function Accelerate ( file : file_type;
                        v0,v1,v2 : Vector ) return Vector is

    res,d1,d2,dd,decrement : Vector(v0'range);
    numerator,denominator,s : Complex_Number;

  begin
    d2 := v2 - v1; 
    d1 := v1 - v0;
    dd := d2 - d1;
    put(file,"  |d1| ="); put(file,Norm2(d1),3);
    put(file,"  |d2| ="); put(file,Norm2(d2),3);
    put(file,"  |dd| ="); put(file,Norm2(dd),3);
    numerator := Conjugate_Product(d1,d1);
    denominator := Conjugate_Product(d1,dd);
    if ((REAL_PART(denominator) = 0.0) and (IMAG_PART(denominator) = 0.0))
     then res := v1; put_line(file,"  |d| = 0");
     else s := numerator/denominator;
          decrement := s*d2;
          res := v1 - decrement;
          put(file,"  |d| ="); put(file,Norm2(decrement),3); new_line(file);
    end if;
    return res;
  end Accelerate;

  function Extrapolate ( v : VecVec ) return VecVec is

    res : VecVec(v'first..v'last-2);
    acc : Vector(v(v'first)'range);

  begin
    for i in v'first..v'last-2 loop
      acc := Accelerate(v(i).all,v(i+1).all,v(i+2).all);
      res(i) := new Vector'(acc);
    end loop;
    return res;
  end Extrapolate;

  function Extrapolate ( file : file_type; v : VecVec ) return VecVec is

    res : VecVec(v'first..v'last-2);
    acc : Vector(v(v'first)'range);

  begin
    for i in v'first..v'last-2 loop
      acc := Accelerate(file,v(i).all,v(i+1).all,v(i+2).all);
      res(i) := new Vector'(acc);
    end loop;
    return res;
  end Extrapolate;

  function Extrapolate ( v : List ) return List is

    res,res_last : List;
    first : List := v;
    second,third : List;
    v0,v1,v2 : Link_to_Vector;

  begin
    if not Is_Null(first) then
      second := Tail_Of(first);
      if not Is_Null(second) then
        third := Tail_Of(second);
        while not Is_Null(third) loop
          v0 := Head_Of(first);
          v1 := Head_Of(second);
          v2 := Head_Of(third);
          Append(res,res_last,Accelerate(v0.all,v1.all,v2.all));
          first := second;
          second := third;
          third := Tail_Of(third);
        end loop;
      end if;
    end if;
    return res;
  end Extrapolate;

  function Extrapolate ( file : file_type; v : List ) return List is

    res,res_last : List;
    first : List := v;
    second,third : List;
    v0,v1,v2 : Link_to_Vector;

  begin
    if not Is_Null(first) then
      second := Tail_Of(first);
      if not Is_Null(second) then
        third := Tail_Of(second);
        while not Is_Null(third) loop
          v0 := Head_Of(first);
          v1 := Head_Of(second);
          v2 := Head_Of(third);
          Append(res,res_last,Accelerate(file,v0.all,v1.all,v2.all));
          first := second;
          second := third;
          third := Tail_Of(third);
        end loop;
      end if;
    end if;
    return res;
  end Extrapolate;

  function Extrapolate ( v : List; k : positive ) return List is

    res,ev : List;

  begin
    if k = 1 then 
      res := Extrapolate(v);
    else
      ev := Extrapolate(v);
      if Length_Of(ev) > 2 then
        res := Extrapolate(ev,k-1);
        Clear(ev);
      else
        res := ev;
      end if;
    end if;
    return res;
  end Extrapolate;

end Standard_Aitken_Extrapolation;
