with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;

package body Standard_Dense_Series is

-- CREATORS :

  function Create ( i : integer ) return Series is
  begin
    return Create(i,0);
  end Create;

  function Create ( f : double_float ) return Series is
  begin
    return Create(f,0);
  end Create;

  function Create ( c : Complex_Number ) return Series is
  begin
    return Create(c,0);
  end Create;

  function Create ( i : integer; order : integer32 ) return Series is

    res : Series;

  begin
    res.order := order;
    res.cff(0) := Create(i);
    res.cff(1..order)
      := Standard_Complex_Vectors.Vector'(1..order => Create(0.0));
    return res;
  end Create;

  function Create ( f : double_float; order : integer32 ) return Series is

    res : Series;

  begin
    res.order := order;
    res.cff(0) := Create(f);
    res.cff(1..order)
      := Standard_Complex_Vectors.Vector'(1..order => Create(0.0));
    return res;
  end Create;

  function Create ( c : Complex_Number; order : integer32 ) return Series is

    res : Series;

  begin
    res.order := order;
    res.cff(0) := c;
    res.cff(1..order)
      := Standard_Complex_Vectors.Vector'(1..order => Create(0.0));
    return res;
  end Create;

  function Create ( c : Standard_Complex_Vectors.Vector ) return Series is

    res : Series;

  begin
    if c'last <= max_order then
      res.order := c'last;
      res.cff(0..c'last) := c;
    else
      res.order := max_order;
      res.cff := c(0..max_order);
    end if;
    return res;
  end Create;

  function Create ( s : Series; order : integer32 ) return Series is

    res : Series;
    zero : constant Complex_Number := Create(0.0);

  begin
    res.order := order;
    if order <= s.order then
      for i in 0..res.order loop
        res.cff(i) := s.cff(i);
      end loop;
    else
      for i in 0..s.order loop
        res.cff(i) := s.cff(i);
      end loop;
      for i in s.order+1..order loop
        res.cff(i) := zero;
      end loop;
    end if;
    return res;
  end Create;

-- EQUALITY AND COPY :

  function Equal ( s,t : Series ) return boolean is

    zero : constant Complex_Number := Create(0.0);

  begin
    if s.order <= t.order then
      for i in 0..s.order loop
        if not Standard_Complex_Numbers.Equal(s.cff(i),t.cff(i))
         then return false;
        end if;
      end loop;
      for i in s.order+1..t.order loop
        if not Standard_Complex_Numbers.Equal(t.cff(i),zero)
         then return false;
        end if;
      end loop;
      return true;
    else
      return Standard_Dense_Series.Equal(t,s);
    end if;
  end Equal;

  procedure Copy ( s : in Series; t : in out Series ) is

    zero : constant Complex_Number := Create(0.0);

  begin
    t.order := s.order;
    for i in 0..s.order loop
      t.cff(i) := s.cff(i);
    end loop;
  end Copy;

-- ARITHMETICAL OPERATORS :

  function "+" ( s : Series; c : Complex_Number ) return Series is

    res : Series := s;

  begin
    res.cff(0) := s.cff(0) + c;
    return res;
  end "+";

  function "+" ( c : Complex_Number; s : Series ) return Series is

    res : Series := s;
  
  begin
    res.cff(0) := c + s.cff(0);
    return res;
  end "+";

  procedure Add ( s : in out Series; c : in Complex_Number ) is
  begin
    s.cff(0) := s.cff(0) + c;
  end Add;

  function "+" ( s : Series ) return Series is

    res : constant Series := s;

  begin
    return res;
  end "+";

  function "+" ( s,t : Series ) return Series is

    res : Series;

  begin
    if s.order = t.order then
      res.order := s.order;
      for i in 0..res.order loop
        res.cff(i) := s.cff(i) + t.cff(i);
      end loop;
    elsif s.order < t.order then -- order of result is t.order
      res.order := t.order;
      for i in 0..s.order loop
        res.cff(i) := s.cff(i) + t.cff(i);
      end loop;
      for i in s.order+1..t.order loop -- copy remaining terms from t
        res.cff(i) := t.cff(i);
      end loop;
    else -- s.order > t.order and the order of result is s.order
      res.order := s.order;
      for i in 0..t.order loop
        res.cff(i) := s.cff(i) + t.cff(i);
      end loop;
      for i in t.order+1..s.order loop -- copy remaining terms from s
        res.cff(i) := s.cff(i);
      end loop;
    end if;
    return res;
  end "+";

  procedure Add ( s : in out Series; t : in Series ) is
  begin
    if t.order >= s.order then  -- do not ignore terms of order > s.order!
      for i in 0..s.order loop
        s.cff(i) := s.cff(i) + t.cff(i);
      end loop;
      if t.order > s.order then 
        for i in s.order+1..t.order loop -- copy higher order terms
          s.cff(i) := t.cff(i);
        end loop;
        s.order := t.order; -- adjust the order of s
      end if;
    else          -- add only the coefficients of index <= t.order
      for i in 0..t.order loop
        s.cff(i) := s.cff(i) + t.cff(i);
      end loop;
    end if;
  end Add;

  function "-" ( s : Series; c : Complex_Number ) return Series is

    res : Series := s;

  begin
    res.cff(0) := s.cff(0) - c;
    return res;
  end "-";

  function "-" ( c : Complex_Number; s : Series ) return Series is

    res : Series;

  begin
    res.order := s.order;
    res.cff(0) := c - s.cff(0);
    for k in 1..res.order loop
      res.cff(k) := -s.cff(k);
    end loop;
    return res;
  end "-";

  procedure Sub ( s : in out Series; c : in Complex_Number ) is
  begin
    s.cff(0) := s.cff(0) - c;
  end Sub;

  function "-" ( s : Series ) return Series is

    res : Series;

  begin
    res.order := s.order;
    for i in 0..res.order loop
      res.cff(i) := -s.cff(i);
    end loop;
    return res;
  end "-";

  procedure Min ( s : in out Series ) is
  begin
    for i in 0..s.order loop
      s.cff(i) := -s.cff(i);
    end loop;
  end Min;

  function "-" ( s,t : Series ) return Series is

    res : Series;

  begin
    if s.order = t.order then
      res.order := s.order;
      for i in 0..t.order loop
        res.cff(i) := s.cff(i) - t.cff(i);
      end loop;
    elsif s.order < t.order then -- the order of the result is t.order
      res.order := t.order;
      for i in 0..s.order loop
        res.cff(i) := s.cff(i) - t.cff(i);
      end loop;
      for i in s.order+1..t.order loop -- copy remaining terms
        res.cff(i) := -t.cff(i);       -- of t with minus sign
      end loop;
    else
      res.order := s.order;
      for i in 0..s.order loop
        res.cff(i) := s.cff(i) - t.cff(i);
      end loop;
      for i in t.order+1..s.order loop -- copy remaining terms of s
        res.cff(i) := s.cff(i);
      end loop;
    end if;
    return res;
  end "-";

  procedure Sub ( s : in out Series; t : in Series ) is
  begin
    if t.order >= s.order then -- do not ignore terms of t of index > s.order!
      for i in 0..s.order loop
        s.cff(i) := s.cff(i) - t.cff(i);
      end loop;
      if t.order > s.order then
        for i in s.order+1..t.order loop -- subtract higher order terms
          s.cff(i) := -t.cff(i);
        end loop;
        s.order := t.order; -- s become a series of higher order
      end if;
    else
      for i in 0..t.order loop
        s.cff(i) := s.cff(i) - t.cff(i);
      end loop;
    end if;
  end Sub;

  function "*" ( s : Series; c : Complex_Number ) return Series is

    res : Series;

  begin
    res.order := s.order;
    for k in 0..s.order loop
      res.cff(k) := s.cff(k)*c;
    end loop;
    return res;
  end "*";

  function "*" ( c : Complex_Number; s : Series ) return Series is

    res : Series;

  begin
    res.order := s.order;
    for k in 0..s.order loop
      res.cff(k) := c*s.cff(k);
    end loop;
    return res;
  end "*";

  procedure Mul ( s : in out Series; c : in Complex_Number ) is
  begin
    for i in 0..s.order loop
      s.cff(i) := s.cff(i)*c;
    end loop;
  end Mul;

  function "*" ( s,t : Series ) return Series is

    res : Series;

  begin
    if s.order = t.order then
      res.order := s.order;
      for i in 0..res.order loop
        res.cff(i) := s.cff(0)*t.cff(i);
        for j in 1..i loop
          res.cff(i) := res.cff(i) + s.cff(j)*t.cff(i-j);
        end loop;
      end loop;
    elsif s.order < t.order then
      res.order := t.order;
      for i in 0..res.order loop
        res.cff(i) := s.cff(0)*t.cff(i);
        for j in 1..i loop
          exit when j > s.order;
          res.cff(i) := res.cff(i) + s.cff(j)*t.cff(i-j);
        end loop;
      end loop;
    else -- s.order > t.order, we then flip s and t
      res.order := s.order;
      for i in 0..res.order loop
        res.cff(i) := t.cff(0)*s.cff(i);
        for j in 1..i loop
          exit when j > t.order;
          res.cff(i) := res.cff(i) + t.cff(j)*s.cff(i-j);
        end loop;
      end loop;
    end if;
    return res;
  end "*";

  procedure Mul ( s : in out Series; t : in Series ) is

    res : constant Series := s*t;

  begin
    s := res;
  end Mul;

  function Inverse ( s : Series ) return Series is

    res : Series;

  begin
    res.order := s.order;
    res.cff(0) := 1.0/s.cff(0);
    for i in 1..res.order loop
      res.cff(i) := -s.cff(1)*res.cff(i-1);
      for j in 2..i loop
        res.cff(i) := res.cff(i) - s.cff(j)*res.cff(i-j);
      end loop;
      res.cff(i) := res.cff(i)/s.cff(0);
    end loop;
    return res;
  end Inverse;

  function "/" ( s : Series; c : Complex_Number ) return Series is

    res : Series;

  begin
    res.order := s.order;
    for k in 0..s.order loop
      res.cff(k) := s.cff(k)/c;
    end loop;
    return res;
  end "/";

  function "/" ( c : Complex_Number; s : Series ) return Series is
  begin
    return c*Inverse(s);
  end "/";

  procedure Div ( s : in out Series; c : in Complex_Number ) is
  begin
    for k in 0..s.order loop
      s.cff(k) := s.cff(k)/c;
    end loop;
  end Div;

  function "/" ( s,t : Series ) return Series is
  begin
    return s*Inverse(t);
  end "/";

  procedure Div ( s : in out Series; t : in Series ) is

    invt : constant Series := Inverse(t);

  begin
    Mul(s,invt);
  end Div;

  function "**" ( s : Series; p : integer ) return Series is

    res : Series;

  begin
    if p = 0 then
      res := Create(1);
    elsif p > 0 then
      res := s;
      for k in 2..p loop
        Mul(res,s);
      end loop;
    else -- p < 0
      res := s;
      for k in 2..(-p) loop
        Mul(res,s);
      end loop;
      res := Inverse(res);
    end if;
    return res;
  end "**";

  function "**" ( s : Series; p : natural32 ) return Series is

    res : Series;

  begin
    if p = 0 then
      res := Create(1);
    else
      res := s;
      for k in 2..p loop
        Mul(res,s);
      end loop;
    end if;
    return res;
  end "**";

-- EVALUATORS :

  function Eval ( s : Series; t : double_float ) return Complex_Number is

    res : Complex_Number := s.cff(0);
    pwt : double_float := t;

  begin
    for i in 1..(s.order-1) loop
      res := res + s.cff(i)*pwt;
      pwt := pwt*t;
    end loop;
    res := res + s.cff(s.order)*pwt;
    return res;
  end Eval;

  function Eval ( s : Series; t : Complex_Number ) return Complex_Number is

    res : Complex_Number := s.cff(0);
    pwt : Complex_Number := t;

  begin
    for i in 1..(s.order-1) loop
      res := res + s.cff(i)*pwt;
      pwt := pwt*t;
    end loop;
    res := res + s.cff(s.order)*pwt;
    return res;
  end Eval;

  procedure Filter ( s : in out Series; tol : in double_float ) is
  begin
    for i in 0..s.order loop
      if AbsVal(s.cff(i)) < tol
       then s.cff(i) := Create(0.0);
      end if;
    end loop;
  end Filter;

-- DESTRUCTOR :

  procedure Clear ( s : in out Series ) is
  begin
    s.order := 0;
    for i in s.cff'range loop
      s.cff(i) := Create(0.0);
    end loop;
  end Clear;

end Standard_Dense_Series;
