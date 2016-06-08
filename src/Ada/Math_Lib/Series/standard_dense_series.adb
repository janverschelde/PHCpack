with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;

package body Standard_Dense_Series is

-- CREATORS :

  function Create ( i : integer ) return Series is

    res : Series(0..0);

  begin
    res(0) := Create(i);
    return res;
  end Create;

  function Create ( f : double_float ) return Series is

    res : Series(0..0);

  begin
    res(0) := Create(f);
    return res;
  end Create;

  function Create ( c : Complex_Number ) return Series is

    res : Series(0..0);

  begin
    res(0) := c;
    return res;
  end Create;

  function Create ( i : integer; order : integer32 ) return Series is

    res : Series(0..order);

  begin
    res(0) := Create(i);
    return res;
  end Create;

  function Create ( f : double_float; order : integer32 ) return Series is

    res : Series(0..order);

  begin
    res(0) := Create(f);
    res(1..order) := Series'(1..order => Create(0.0));
    return res;
  end Create;

  function Create ( c : Complex_Number; order : integer32 ) return Series is

    res : Series(0..order);

  begin
    res(0) := c;
    return res;
  end Create;

  function Create ( s : Series; order : integer32 ) return Series is

    res : Series(0..order);
    zero : constant Complex_Number := Create(0.0);

  begin
    if order <= s'last then
      for i in res'range loop
        res(i) := s(i);
      end loop;
    else
      for i in s'range loop
        res(i) := s(i);
      end loop;
      for i in s'last+1..order loop
        res(i) := zero;
      end loop;
    end if;
    return res;
  end Create;

-- EQUALITY AND COPY :

  function Equal ( s,t : Series ) return boolean is

    zero : constant Complex_Number := Create(0.0);

  begin
    if s'last <= t'last then
      for i in s'range loop
        if not Standard_Complex_Numbers.Equal(s(i),t(i))
         then return false;
        end if;
      end loop;
      for i in s'last+1..t'last loop
        if not Standard_Complex_Numbers.Equal(t(i),zero)
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
    for i in s'range loop
      t(i) := s(i);
    end loop;
    if t'last > s'last then
      for i in s'last+1..t'last loop
        t(i) := zero;
      end loop;
    end if;
  end Copy;

-- ARITHMETICAL OPERATORS :

  function "+" ( s,t : Series ) return Series is

    res : Series(0..s'last);

  begin
    for i in s'range loop
      res(i) := s(i) + t(i);
    end loop;
    return res;
  end "+";

  procedure Add ( s : in out Series; t : in Series ) is
  begin
    for i in s'range loop
      s(i) := s(i) + t(i);
    end loop;
  end Add;

  function "-" ( s : Series ) return Series is

    res : Series(0..s'last);

  begin
    for i in s'range loop
      res(i) := -s(i);
    end loop;
    return res;
  end "-";

  procedure Min ( s : in out Series ) is
  begin
    for i in s'range loop
      s(i) := -s(i);
    end loop;
  end Min;

  function "-" ( s,t : Series ) return Series is

    res : Series(0..s'last);

  begin
    for i in s'range loop
      res(i) := s(i) - t(i);
    end loop;
    return res;
  end "-";

  procedure Min ( s : in out Series; t : in Series ) is
  begin
    for i in s'range loop
      s(i) := s(i) - t(i);
    end loop;
  end Min;

  function "*" ( s,t : Series ) return Series is

    res : Series(0..s'last);

  begin
    for i in res'range loop
      res(i) := s(0)*t(i);
      for j in 1..i loop
        res(i) := res(i) + s(j)*t(i-j);
      end loop;
    end loop;
    return res;
  end "*";

  procedure Mul ( s : in out Series; t : in Series ) is

    res : Series(0..s'last);

  begin
    for i in res'range loop
      res(i) := s(0)*t(i);
      for j in 1..i loop
        res(i) := res(i) + s(j)*t(i-j);
      end loop;
    end loop;
    s := res;
  end Mul;

  function Inverse ( s : Series ) return Series is

    res : Series(0..s'last);

  begin
    res(0) := 1.0/s(0);
    for i in 1..res'last loop
      res(i) := -s(1)*res(i-1);
      for j in 2..i loop
        res(i) := res(i) - s(j)*res(i-j);
      end loop;
      res(i) := res(i)/s(0);
    end loop;
    return res;
  end Inverse;

  function "/" ( s,t : Series ) return Series is

    invt : constant Series(t'range) := Inverse(t); 

  begin
    return s*invt;
  end "/";

  procedure Div ( s : in out Series; t : in Series ) is

    invt : constant Series(t'range) := Inverse(t);

  begin
    Mul(s,invt);
  end Div;

-- DESTRUCTOR :

  procedure Clear ( s : in out Series ) is
  begin
    for i in s'range loop
      s(i) := Create(0.0);
    end loop;
  end Clear;

end Standard_Dense_Series;
