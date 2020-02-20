with Standard_Mathematical_Functions;
with Standard_Complex_Numbers_Polar;
with Binomial_Coefficients;

package body Standard_Complex_Series_Functions is

  function Eval ( s : Series; t : double_float ) return Complex_Number is

    res : Complex_Number := s.cff(0);
    pwt : double_float := t;

  begin
    for i in 1..(s.deg-1) loop
      res := res + s.cff(i)*pwt;
      pwt := pwt*t;
    end loop;
    res := res + s.cff(s.deg)*pwt;
    return res;
  end Eval;

  function Eval ( s : Series; t : Complex_Number ) return Complex_Number is

    res : Complex_Number := s.cff(0);
    pwt : Complex_Number := t;

  begin
    for i in 1..(s.deg-1) loop
      res := res + s.cff(i)*pwt;
      pwt := pwt*t;
    end loop;
    res := res + s.cff(s.deg)*pwt;
    return res;
  end Eval;

  function Eval ( s : Link_to_Series;
                  t : double_float ) return Complex_Number is
  begin
    if s = null
     then return Create(0.0);
     else return Eval(s.all,t);
    end if;
  end Eval;

  function Eval ( s : Link_to_Series;
                  t : Complex_Number ) return Complex_Number is
  begin
    if s = null
     then return Create(0.0);
     else return Eval(s.all,t);
    end if;
  end Eval;

  function Eval ( s : Series; t : double_float;
                  a,b : integer32 ) return Complex_Number is

    use Standard_Mathematical_Functions;

    pow : double_float := double_float(a)/double_float(b);
    pwt : double_float := t**pow;
    res : Complex_Number := s.cff(0)*pwt;

  begin
    for i in 1..s.deg loop
      pow := double_float(a+i)/double_float(b);
      pwt := t**pow;
      res := res + s.cff(i)*pwt;
    end loop;
    return res;
  end Eval;

  function Eval ( s : Series; t : Complex_Number;
                  a,b : integer32 ) return Complex_Number is

    use Standard_Complex_Numbers_Polar;

    pow : double_float := double_float(a)/double_float(b);
    pwt : Complex_Number := Polar_Exponentiation(t,pow);
    res : Complex_Number := s.cff(0)*pwt;

  begin
    for i in 1..s.deg loop
      pow := double_float(a+i)/double_float(b);
      pwt := Polar_Exponentiation(t,pow);
      res := res + s.cff(i)*pwt;
    end loop;
    return res;
  end Eval;

  function Eval ( s : Link_to_Series; t : double_float;
                  a,b : integer32 ) return Complex_Number is
  begin
    if s = null
     then return Create(0.0);
     else return Eval(s.all,t,a,b);
    end if;
  end Eval;

  function Eval ( s : Link_to_Series; t : Complex_Number;
                  a,b : integer32 ) return Complex_Number is
  begin
    if s = null
     then return Create(0.0);
     else return Eval(s.all,t,a,b);
    end if;
  end Eval;

-- ORDER and FILTER :

  function Order ( s : Series; tol : double_float := 0.0 ) return integer32 is
  begin
    for k in 0..s.deg loop
      if Standard_Complex_Numbers.AbsVal(s.cff(k)) > tol
       then return k;
      end if;
    end loop;
    return s.deg+1;
  end Order;

  function Order ( s : Link_to_Series;
                   tol : double_float := 0.0 ) return integer32 is
  begin
    if s = null
     then return -1;
     else return Order(s.all,tol);
    end if;
  end Order;

  procedure Filter ( s : in out Series; tol : in double_float ) is
  begin
    for i in 0..s.deg loop
      if AbsVal(s.cff(i)) < tol
       then s.cff(i) := Create(0.0);
      end if;
    end loop;
  end Filter;

  procedure Filter ( s : in Link_to_Series; tol : in double_float ) is
  begin
    if s /= null
     then Filter(s.all,tol);
    end if;
  end Filter;

-- SHIFT OPERATORS :

  function Shift ( s : Series; c : double_float ) return Series is

    res : Series(s.deg);
    bcf : double_float;
    sgn : integer32;

    use Binomial_Coefficients;

  begin
    for i in 0..s.deg loop
      res.cff(i) := Create(0.0);
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
       -- bcf := double_float(sgn*binomial(i,j)); -- causes overflow
        bcf := double_float(sgn)*binomial(i,j);
        bcf := bcf*(c**(natural(i-j)));
        res.cff(j) := res.cff(j) + s.cff(i)*bcf;
        sgn := -sgn;
      end loop;
    end loop;
    return res;
  end Shift;

  function Shift ( s : Series; c : Complex_Number ) return Series is

    res : Series(s.deg);
    bcf : double_float;
    rcf : Complex_Number;
    sgn : integer32;

    use Binomial_Coefficients;

  begin
    for i in 0..s.deg loop
      res.cff(i) := Create(0.0);
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
       -- bcf := double_float(sgn*binomial(i,j)); -- causes overflow
        bcf := double_float(sgn)*binomial(i,j);
        rcf := bcf*(c**(natural(i-j)));
        res.cff(j) := res.cff(j) + s.cff(i)*rcf;
        sgn := -sgn;
      end loop;
    end loop;
    return res;
  end Shift;

  function Shift ( s : Link_to_Series;
                   c : double_float ) return Link_to_Series is
  begin
    if s = null
     then return null;
     else return new Series'(Shift(s.all,c));
    end if;
  end Shift;

  function Shift ( s : Link_to_Series;
                   c : Complex_Number ) return Link_to_Series is
  begin
    if s = null
     then return null;
     else return new Series'(Shift(s.all,c));
    end if;
  end Shift;

  procedure Shift ( s : in out Series; c : in double_float ) is
 
    r : constant Series(s.deg) := Shift(s,c);
   
  begin
    s := r;
  end Shift;

  procedure Shift ( s : in out Series; c : in Complex_Number ) is
 
    r : constant Series(s.deg) := Shift(s,c);
   
  begin
    s := r;
  end Shift;

  procedure Shift ( s : in Link_to_Series; c : in double_float ) is
  begin
    if s /= null then
      declare
        r : constant Series(s.deg) := Shift(s.all,c);
      begin
        s.cff := r.cff;
      end;
    end if;
  end Shift;

  procedure Shift ( s : in Link_to_Series; c : in Complex_Number ) is
  begin
    if s /= null then 
      declare
        r : constant Series(s.deg) := Shift(s.all,c);
      begin
        s.cff := r.cff;
      end;
    end if;
  end Shift;

end Standard_Complex_Series_Functions;
