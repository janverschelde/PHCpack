with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;

package body Multprec_Aitken_Extrapolation is

  function Conjugate_Product ( x,y : Vector ) return Complex_Number is

    res : Complex_Number := Create(integer(0));
    acc : Complex_Number;

  begin
    for i in x'range loop           -- res := res + x(i)*Conjugate(y(i));
      acc := Conjugate(y(i));
      Mul(acc,x(i));
      Add(res,acc);
      Clear(acc);
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
    re,im : Floating_Number;

  begin
    d2 := v2 - v1;
    d1 := v1 - v0;
    dd := d2 - d1;
    numerator := Conjugate_Product(d1,d1);
    denominator := Conjugate_Product(d1,dd);
    re := REAL_PART(denominator);
    im := REAL_PART(denominator);
    if Equal(re,0.0) and Equal(im,0.0) then
      Copy(v1,res);
    else
      s := numerator/denominator;
      decrement := s*d2;
      res := v1 - decrement;
      Clear(s); Clear(decrement);
    end if;
    Clear(d1); Clear(d2); Clear(dd);
    Clear(numerator); Clear(denominator);
    Clear(re); Clear(im);
    return res;
  end Accelerate;

  function Accelerate ( file : file_type;
                        v0,v1,v2 : Vector ) return Vector is

    res,d1,d2,dd,decrement : Vector(v0'range);
    numerator,denominator,s : Complex_Number;
    re,im,nrm : Floating_Number;

  begin
    d2 := v2 - v1;
    d1 := v1 - v0;
    dd := d2 - d1;
    put(file,"  |d1| ="); nrm := Norm2(d1); put(file,nrm,3); Clear(nrm);
    put(file,"  |d2| ="); nrm := Norm2(d2); put(file,nrm,3); Clear(nrm);
    put(file,"  |dd| ="); nrm := Norm2(dd); put(file,nrm,3); Clear(nrm);
    numerator := Conjugate_Product(d1,d1);
    denominator := Conjugate_Product(d1,dd);
    re := REAL_PART(denominator);
    im := IMAG_PART(denominator);
    if Equal(re,0.0) and Equal(im,0.0)
     then Copy(v1,res); put_line(file,"  |d| = 0");
     else s := numerator/denominator;
          decrement := s*d2;
          res := v1 - decrement;
          nrm := Norm2(decrement);
          put(file,"  |d| = "); put(file,nrm,3); new_line(file);
          Clear(nrm); Clear(s); Clear(decrement);
    end if;
    Clear(d1); Clear(d2); Clear(dd);
    Clear(numerator); Clear(denominator);
    Clear(re); Clear(im);
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

end Multprec_Aitken_Extrapolation;
