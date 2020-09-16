with PentDobl_Complex_Numbers_Polar;
with Binomial_Coefficients;

package body PentDobl_Complex_Series_Functions is

  function Eval ( s : Series; t : penta_double ) return Complex_Number is

    res : Complex_Number := s.cff(0);
    pwt : penta_double := t;

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
                  t : penta_double ) return Complex_Number is

    zero : constant penta_double := create(0.0);

  begin
    if s = null
     then return Create(zero);
     else return Eval(s.all,t);
    end if;
  end Eval;

  function Eval ( s : Link_to_Series;
                  t : Complex_Number ) return Complex_Number is

    zero : constant penta_double := create(0.0);

  begin
    if s = null
     then return Create(zero);
     else return Eval(s.all,t);
    end if;
  end Eval;

  function Eval ( s : Series; t : penta_double;
                  a,b : integer32 ) return Complex_Number is

    dd_a : penta_double := create(a);
    dd_b : constant penta_double := create(b);
    pow : penta_double := dd_a/dd_b;
    pwt : penta_double := t**pow;
    res : Complex_Number := s.cff(0)*pwt;

  begin
    for i in 1..s.deg loop
      dd_a := create(a+i);
      pow := dd_a/dd_b;
      pwt := t**pow;
      res := res + s.cff(i)*pwt;
    end loop;
    return res;
  end Eval;

  function Eval ( s : Series; t : Complex_Number;
                  a,b : integer32 ) return Complex_Number is

    use PentDobl_Complex_Numbers_Polar;

    dd_a : penta_double := create(a);
    dd_b : constant penta_double := create(b);
    pow : penta_double := dd_a/dd_b;
    pwt : Complex_Number := Polar_Exponentiation(t,pow);
    res : Complex_Number := s.cff(0)*pwt;

  begin
    for i in 1..s.deg loop
      dd_a := create(a+i);
      pow := dd_a/dd_b;
      pwt := Polar_Exponentiation(t,pow);
      res := res + s.cff(i)*pwt;
    end loop;
    return res;
  end Eval;

  function Eval ( s : Link_to_Series; t : penta_double;
                  a,b : integer32 ) return Complex_Number is

    zero : constant penta_double := create(0.0);

  begin
    if s = null
     then return Create(zero);
     else return Eval(s.all,t,a,b);
    end if;
  end Eval;

  function Eval ( s : Link_to_Series; t : Complex_Number;
                  a,b : integer32 ) return Complex_Number is

    zero : constant penta_double := create(0.0);

  begin
    if s = null
     then return Create(zero);
     else return Eval(s.all,t,a,b);
    end if;
  end Eval;

-- ORDER and FILTER :

  function Order ( s : Series; tol : double_float := 0.0 ) return integer32 is
  begin
    for k in 0..s.deg loop
      if PentDobl_Complex_Numbers.AbsVal(s.cff(k)) > tol
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

    zero : constant penta_double := create(0.0);

  begin
    for i in 0..s.deg loop
      if AbsVal(s.cff(i)) < tol
       then s.cff(i) := Create(zero);
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

  function Shift ( s : Series; c : penta_double ) return Series is

    res : Series(s.deg);
    bcf : penta_double;
    sgn : integer32;

    use Binomial_Coefficients;

  begin
    for i in 0..s.deg loop
      res.cff(i) := Create(integer32(0));
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
       -- bcf := penta_double_numbers.create(sgn*binomial(i,j));
        bcf := binomial(i,j);
        bcf := penta_double_numbers.create(sgn)*bcf;
        bcf := bcf*(c**(natural(i-j)));
        res.cff(j) := res.cff(j) + s.cff(i)*bcf;
        sgn := -sgn;
      end loop;
    end loop;
    return res;
  end Shift;

  function Shift ( s : Series; c : Complex_Number ) return Series is

    res : Series(s.deg);
    bcf : penta_double;
    rcf : Complex_Number;
    sgn : integer32;

    use Binomial_Coefficients;

  begin
    for i in 0..s.deg loop
      res.cff(i) := Create(integer32(0));
      if i mod 2 = 0
       then sgn := 1;
       else sgn := -1;
      end if;
      for j in 0..i loop
       -- bcf := penta_double_numbers.create(sgn*binomial(i,j));
        bcf := binomial(i,j);
        bcf := penta_double_numbers.create(sgn)*bcf;
        rcf := bcf*(c**(natural(i-j)));
        res.cff(j) := res.cff(j) + s.cff(i)*rcf;
        sgn := -sgn;
      end loop;
    end loop;
    return res;
  end Shift;

  function Shift ( s : Link_to_Series;
                   c : penta_double ) return Link_to_Series is
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

  procedure Shift ( s : in out Series; c : in penta_double ) is
 
    r : constant Series := Shift(s,c);
   
  begin
    s := r;
  end Shift;

  procedure Shift ( s : in out Series; c : in Complex_Number ) is
 
    r : constant Series := Shift(s,c);
   
  begin
    s := r;
  end Shift;

  procedure Shift ( s : in Link_to_Series; c : in penta_double ) is
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

end PentDobl_Complex_Series_Functions;
