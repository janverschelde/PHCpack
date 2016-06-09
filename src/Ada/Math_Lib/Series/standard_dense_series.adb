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
    for i in 0..s.order loop
      t.cff(i) := s.cff(i);
    end loop;
    if t.order > s.order then
      for i in s.order+1..t.order loop
        t.cff(i) := zero;
      end loop;
    end if;
  end Copy;

-- ARITHMETICAL OPERATORS :

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
    if t.order >= s.order then  -- ignore terms of order > s.order
      for i in 0..s.order loop
        s.cff(i) := s.cff(i) + t.cff(i);
      end loop;
    else          -- add only the coefficients of index <= t.order
      for i in 0..t.order loop
        s.cff(i) := s.cff(i) + t.cff(i);
      end loop;
    end if;
  end Add;

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
    if t.order >= s.order then -- ignore terms of t of index > s.order
      for i in 0..s.order loop
        s.cff(i) := s.cff(i) - t.cff(i);
      end loop;
    else
      for i in 0..t.order loop
        s.cff(i) := s.cff(i) - t.cff(i);
      end loop;
    end if;
  end Sub;

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

    res : Series;

  begin
    res.order := s.order;
    for i in 0..res.order loop
      res.cff(i) := t.cff(0)*s.cff(i);
      for j in 1..i loop
        exit when j > t.order;
        res.cff(i) := res.cff(i) + t.cff(j)*s.cff(i-j);
      end loop;
    end loop;
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

  function "/" ( s,t : Series ) return Series is

    invt : constant Series := Inverse(t); 

  begin
    return s*invt;
  end "/";

  procedure Div ( s : in out Series; t : in Series ) is

    invt : constant Series := Inverse(t);

  begin
    Mul(s,invt);
  end Div;

-- DESTRUCTOR :

  procedure Clear ( s : in out Series ) is
  begin
    s.order := 0;
    for i in s.cff'range loop
      s.cff(i) := Create(0.0);
    end loop;
  end Clear;

end Standard_Dense_Series;
