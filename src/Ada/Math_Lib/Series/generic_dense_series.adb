with unchecked_deallocation;

package body Generic_Dense_Series is 

  use Ring; -- for the arithmetical operators

-- CONSTRUCTORS :

  function Create ( i : integer ) return Series is
 
    res : Series(0);

  begin
    res.cff(0) := Ring.Create(i);
    return res;
  end Create;

  function Create ( n : Ring.number ) return Series is
 
    res : Series(0);

  begin
    res.cff(0) := n;
    return res;
  end Create;

  function Create ( i : integer ) return Link_to_Series is

    ser : constant Series(0) := Create(i);
    res : constant Link_to_Series := new Series'(ser);

  begin
    return res;
  end Create;

  function Create ( n : Ring.number ) return Link_to_Series is

    ser : constant Series(0) := Create(n);
    res : constant Link_to_Series := new Series'(ser);

  begin
    return res;
  end Create;

  function Create ( i : integer; deg : integer32 ) return Series is

    res : Series(deg);

  begin
    res.cff(0) := Ring.Create(i);
    for k in 1..deg loop
      res.cff(k) := Ring.zero;
    end loop;
    return res;
  end Create;

  function Create ( n : Ring.number; deg : integer32 ) return Series is

    res : Series(deg);

  begin
    res.cff(0) := n;
    for k in 1..deg loop
      res.cff(k) := Ring.zero;
    end loop;
    return res;
  end Create;

  function Create ( i : integer; deg : integer32 ) return Link_to_Series is

    ser : constant Series(deg) := Create(i,deg);
    res : constant Link_to_Series := new Series'(ser);

  begin
    return res;
  end Create;

  function Create ( n : Ring.number;
                    deg : integer32 ) return Link_to_Series is

    ser : constant Series(deg) := Create(n,deg);
    res : constant Link_to_Series := new Series'(ser);

  begin
    return res;
  end Create;

  function Create ( c : Vectors.Vector ) return Series is

    res : Series(c'last);

  begin
    res.cff(0..c'last) := c;
    return res;
  end Create;

  function Create ( c : Vectors.Vector ) return Link_to_Series is

    ser : constant Series(c'last) := Create(c);
    res : constant Link_to_Series := new Series'(ser);

  begin
    return res;
  end Create;

  function Create ( s : Series; deg : integer32 ) return Series is

    res : Series(deg);

  begin
    if deg <= s.deg then
      for i in 0..res.deg loop
        res.cff(i) := s.cff(i);
      end loop;
    else
      for i in 0..s.deg loop
        res.cff(i) := s.cff(i);
      end loop;
      for i in s.deg+1..deg loop
        res.cff(i) := Ring.zero;
      end loop;
    end if;
    return res;
  end Create;

  function Create ( s : Series; deg : integer32 )
                  return Link_to_Series is

    ser : constant Series(deg) := Create(s,deg);
    res : constant Link_to_Series := new Series'(ser);

  begin
    return res;
  end Create;

  procedure Set_Degree ( s : in out Link_to_Series; deg : in integer32 ) is
  begin
    if s.deg /= deg then
      declare
        res : constant Link_to_Series := Create(s.all,deg);
      begin
        Clear(s);
        s := res;
      end;
    end if;
  end Set_Degree;

-- EQUALITY AND COPY :

  function Equal ( s,t : Series ) return boolean is
  begin
    if s.deg <= t.deg then
      for i in 0..s.deg loop
        if not Ring.Equal(s.cff(i),t.cff(i))
         then return false;
        end if;
      end loop;
      for i in s.deg+1..t.deg loop
        if not Ring.Equal(t.cff(i),Ring.zero)
         then return false;
        end if;
      end loop;
      return true;
    else
      return Generic_Dense_Series.Equal(s=>t,t=>s);
    end if;
  end Equal;

  function Equal ( s,t : Link_to_Series ) return boolean is
  begin
    if s = null then
      if t = null
       then return true;
       else return false;
      end if;
    elsif t = null then
      return false;
    else
      return Equal(s.all,t.all);
    end if;
  end Equal;

  procedure Copy ( s : in Series; t : in out Series ) is
  begin
    for i in 0..s.deg loop
      exit when (i > t.deg);
      t.cff(i) := s.cff(i);
    end loop;
  end Copy;

  procedure Copy ( s : in Link_to_Series; t : in out Link_to_Series ) is
  begin
    Clear(t);
    t := Create(s.cff);
  end Copy;

-- ARITHMETICAL OPERATORS :

  function "+" ( s : Series; c : Ring.number ) return Series is

    res : Series := s;

  begin
    res.cff(0) := s.cff(0) + c;
    return res;
  end "+";

  function "+" ( s : Link_to_Series;
                 c : Ring.number ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s = null then
      res := Create(c);
    else
      res := Create(s.cff);
      res.cff(0) := res.cff(0) + c;
    end if;
    return res;
  end "+";

  function "+" ( c : Ring.number; s : Series ) return Series is

    res : Series := s;
  
  begin
    res.cff(0) := c + s.cff(0);
    return res;
  end "+";

  function "+" ( c : Ring.number;
                 s : Link_to_Series ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s = null then
      res := Create(c);
    else
      res := Create(s.cff);
      res.cff(0) := res.cff(0) + c;
    end if;
    return res;
  end "+";

  procedure Add ( s : in out Series; c : in Ring.number ) is
  begin
    s.cff(0) := s.cff(0) + c;
  end Add;

  procedure Add ( s : in out Link_to_Series;
                  c : in Ring.number ) is
  begin
    if s = null
     then s := Create(c);
     else s.cff(0) := s.cff(0) + c;
    end if;
  end Add;

  function "+" ( s : Series ) return Series is

    res : constant Series(s.deg) := s;

  begin
    return res;
  end "+";

  function "+" ( s : Link_to_Series ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s = null then
      return null;
    else
      res := new Series'(s.all);
      return res;
    end if;
  end "+";

  function "+" ( s,t : Series ) return Series is
  begin
    if s.deg = t.deg then
      declare
        res : Series(s.deg);
      begin
        for i in 0..res.deg loop
          res.cff(i) := s.cff(i) + t.cff(i);
        end loop;
        return res;
      end;
    elsif s.deg < t.deg then -- deg of result is t.deg
      declare
        res : Series(t.deg);
      begin
        for i in 0..s.deg loop
          res.cff(i) := s.cff(i) + t.cff(i);
        end loop;
        for i in s.deg+1..t.deg loop -- copy remaining terms from t
          res.cff(i) := t.cff(i);
        end loop;
        return res;
      end;
    else -- s.deg > t.deg and the deg of result is s.deg
      declare
        res : Series(s.deg);
      begin
        for i in 0..t.deg loop
          res.cff(i) := s.cff(i) + t.cff(i);
        end loop;
        for i in t.deg+1..s.deg loop -- copy remaining terms from s
          res.cff(i) := s.cff(i);
        end loop;
        return res;
      end;
    end if;
  end "+";

  function "+" ( s,t : Link_to_Series ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s = null then
      res := +t; -- return a copy of t
    elsif t = null then
      res := +s; -- return a copy of s
    else
      declare
        spt : constant Series := s.all + t.all;
      begin
        res := new Series'(spt);
      end;
    end if;
    return res;
  end "+";

  procedure Add ( s : in out Series; t : in Series ) is
  begin
    for i in 0..t.deg loop
      exit when (i > s.deg); -- ignore higher degree terms of t
      s.cff(i) := s.cff(i) + t.cff(i);
    end loop;
  end Add;

  procedure Add ( s : in out Link_to_Series; t : in Link_to_Series ) is
  begin
    if t /= null then
      if s = null then
        s := new Series'(t.all);
      else
        if t.deg <= s.deg then
          for i in 0..t.deg loop
            s.cff(i) := s.cff(i) + t.cff(i);
          end loop;
        else
          declare
            spt : Series(t.deg);
          begin
            for i in 0..s.deg loop
              spt.cff(i) := s.cff(i) + t.cff(i);
            end loop;
            for i in s.deg+1..t.deg loop
              spt.cff(i) := t.cff(i);
            end loop;
            Clear(s);
            s := new Series'(spt);
          end;
        end if;
      end if;
    end if;
  end Add;

  function "-" ( s : Series; c : Ring.number ) return Series is

    res : Series(s.deg) := s;

  begin
    res.cff(0) := s.cff(0) - c;
    return res;
  end "-";

  function "-" ( s : Link_to_Series;
                 c : Ring.number ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s = null then
      res := Create(-c);
    else
      res := Create(s.cff);
      res.cff(0) := res.cff(0) - c;
    end if;
    return res;
  end "-";

  function "-" ( c : Ring.number; s : Series ) return Series is

    res : Series(s.deg);

  begin
    res.cff(0) := c - s.cff(0);
    for k in 1..res.deg loop
      res.cff(k) := -s.cff(k);
    end loop;
    return res;
  end "-";

  function "-" ( c : Ring.number;
                 s : Link_to_Series ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s = null then
      res := Create(c);
    else
      res := Create(s.cff);
      res.cff(0) := c - res.cff(0);
      for k in 1..res.deg loop
        res.cff(k) := -res.cff(k);
      end loop;
    end if;
    return res;
  end "-";

  procedure Sub ( s : in out Series; c : in Ring.Number ) is
  begin
    s.cff(0) := s.cff(0) - c;
  end Sub;

  procedure Sub ( s : in out Link_to_Series;
                  c : in Ring.number ) is
  begin
    if s = null
     then s := Create(-c);
     else s.cff(0) := s.cff(0) - c;
    end if;
  end Sub;

  function "-" ( s : Series ) return Series is

    res : Series(s.deg);

  begin
    for i in 0..res.deg loop
      res.cff(i) := -s.cff(i);
    end loop;
    return res;
  end "-";

  function "-" ( s : Link_to_Series ) return Link_to_Series is
  begin
    if s = null
     then return null;
     else return new Series'(-s.all);
    end if;
  end "-";

  procedure Min ( s : in out Series ) is
  begin
    for i in 0..s.deg loop
      s.cff(i) := -s.cff(i);
    end loop;
  end Min;

  procedure Min ( s : in out Link_to_Series ) is
  begin
    if s /= null then
      for i in 0..s.deg loop
        s.cff(i) := -s.cff(i);
      end loop;
    end if;
  end Min;

  function "-" ( s,t : Series ) return Series is
  begin
    if s.deg = t.deg then
      declare
        res : Series(s.deg);
      begin
        for i in 0..t.deg loop
          res.cff(i) := s.cff(i) - t.cff(i);
        end loop;
        return res;
      end;
    elsif s.deg < t.deg then -- the deg of the result is t.deg
      declare
        res : Series(t.deg);
      begin
        for i in 0..s.deg loop
          res.cff(i) := s.cff(i) - t.cff(i);
        end loop;
        for i in s.deg+1..t.deg loop -- copy remaining terms
          res.cff(i) := -t.cff(i);       -- of t with minus sign
        end loop;
        return res;
      end;
    else
      declare
        res : Series(s.deg);
      begin
        for i in 0..s.deg loop
          res.cff(i) := s.cff(i) - t.cff(i);
        end loop;
        for i in t.deg+1..s.deg loop -- copy remaining terms of s
          res.cff(i) := s.cff(i);
        end loop;
        return res;
      end;
    end if;
  end "-";

  function "-" ( s,t : Link_to_Series ) return Link_to_Series is
  begin
    if s = null then
      return -t;
    elsif t = null then
      return +s; -- must return a copy of s, not s
    else
      declare
        smt : constant Series := s.all - t.all;
        res : constant Link_to_Series := new Series'(smt);
      begin
        return res;
      end;
    end if;
  end "-";

  procedure Sub ( s : in out Series; t : in Series ) is
  begin
    for i in 0..t.deg loop
      exit when (i > s.deg); -- ignore higher degree terms of t
      s.cff(i) := s.cff(i) - t.cff(i);
    end loop;
  end Sub;

  procedure Sub ( s : in out Link_to_Series;
                  t : in Link_to_Series ) is
  begin
    if t /= null then
      if s = null then
        s := new Series'(t.all);
        for i in s.cff'range loop  -- s = -t
          s.cff(i) := -s.cff(i);
        end loop;
      else
        if t.deg <= s.deg then
          for i in 0..t.deg loop
            s.cff(i) := s.cff(i) - t.cff(i);
          end loop;
        else
          declare
            smt : Series(t.deg);
          begin
            for i in 0..s.deg loop
              smt.cff(i) := s.cff(i) - t.cff(i);
            end loop;
            for i in s.deg+1..t.deg loop
              smt.cff(i) := -t.cff(i);
            end loop;
            Clear(s);
            s := new Series'(smt);
          end;
        end if;
      end if;
    end if;
  end Sub;

  function "*" ( s : Series; c : Ring.number ) return Series is

    res : Series(s.deg);

  begin
    for k in 0..s.deg loop
      res.cff(k) := s.cff(k)*c;
    end loop;
    return res;
  end "*";

  function "*" ( s : Link_to_Series;
                 c : Ring.number ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s = null
     then res := null;
     else res := new Series'(s.all*c);
    end if;
    return res;
  end "*";

  function "*" ( c : Ring.number; s : Series ) return Series is

    res : Series(s.deg);

  begin
    for k in 0..s.deg loop
      res.cff(k) := c*s.cff(k);
    end loop;
    return res;
  end "*";

  function "*" ( c : Ring.number;
                 s : Link_to_Series ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s = null
     then res := null;
     else res := new Series'(c*s.all);
    end if;
    return res;
  end "*";

  procedure Mul ( s : in out Series; c : in Ring.number ) is
  begin
    for i in 0..s.deg loop
      s.cff(i) := s.cff(i)*c;
    end loop;
  end Mul;

  procedure Mul ( s : in out Link_to_Series; c : in Ring.number ) is
  begin
    if s /= null then
      for i in 0..s.deg loop
        s.cff(i) := s.cff(i)*c;
      end loop;
    end if;
  end Mul;

  function "*" ( s,t : Series ) return Series is
  begin
    if s.deg = t.deg then
      declare
        res : Series(s.deg);
      begin
        for i in 0..res.deg loop
          res.cff(i) := s.cff(0)*t.cff(i);
          for j in 1..i loop
            res.cff(i) := res.cff(i) + s.cff(j)*t.cff(i-j);
          end loop;
        end loop;
        return res;
      end;
    elsif s.deg < t.deg then
      declare
        res : Series(t.deg);
      begin
        for i in 0..res.deg loop
          res.cff(i) := s.cff(0)*t.cff(i);
          for j in 1..i loop
            exit when j > s.deg;
            res.cff(i) := res.cff(i) + s.cff(j)*t.cff(i-j);
          end loop;
        end loop;
        return res;
      end;
    else -- s.deg > t.deg, we then flip s and t
      declare
        res : Series(s.deg);
      begin
        for i in 0..res.deg loop
          res.cff(i) := t.cff(0)*s.cff(i);
          for j in 1..i loop
            exit when j > t.deg;
            res.cff(i) := res.cff(i) + t.cff(j)*s.cff(i-j);
          end loop;
        end loop;
        return res;
      end;
    end if;
  end "*";

  function "*" ( s,t : Link_to_Series ) return Link_to_Series is
  begin
    if s = null or t = null
     then return null;
     else return new Series'(s.all*t.all);
    end if;
  end "*";

  procedure Mul ( s : in out Series; t : in Series ) is

    res : constant Series := s*t;

  begin
    s := res;
  end Mul;

  procedure Mul ( s : in out Link_to_Series; t : in Link_to_Series ) is
  begin
    if s /= null then
      if t = null then
        Clear(s);
      else
        if s.deg >= t.deg then
          for i in reverse 1..s.deg loop
            if i > t.deg
             then s.cff(i) := s.cff(i)*t.cff(0);
             else s.cff(i) := s.cff(i)*t.cff(0) + s.cff(0)*t.cff(i); 
            end if;
            for j in 1..(i-1) loop
              exit when j > t.deg;
              s.cff(i) := s.cff(i) + s.cff(i-j)*t.cff(j);
            end loop;
          end loop;
          s.cff(0) := s.cff(0)*t.cff(0);
        else
          declare
            st : constant Series(t.deg) := s.all*t.all;
          begin
            Clear(s);
            s := new Series'(st);
          end;
        end if;
      end if;
    end if;
  end Mul;

  function Inverse ( s : Series ) return Series is

    res : Series(s.deg);

  begin
    res.cff(0) := one/s.cff(0);
    for i in 1..res.deg loop
      res.cff(i) := -s.cff(1)*res.cff(i-1);
      for j in 2..i loop
        res.cff(i) := res.cff(i) - s.cff(j)*res.cff(i-j);
      end loop;
      res.cff(i) := res.cff(i)/s.cff(0);
    end loop;
    return res;
  end Inverse;

  function Inverse ( s : Link_to_Series ) return Link_to_Series is
  begin
    if s = null
     then return null;
     else return new Series'(Inverse(s.all));
    end if;
  end Inverse;

  function "/" ( s : Series; c : Ring.number ) return Series is

    res : Series(s.deg);

  begin
    for k in 0..s.deg loop
      res.cff(k) := s.cff(k)/c;
    end loop;
    return res;
  end "/";

  function "/" ( s : Link_to_Series;
                 c : Ring.number ) return Link_to_Series is
  begin
    if s = null
     then return null;
     else return new Series'(s.all/c);
    end if;
  end "/";

  function "/" ( c : Ring.number; s : Series ) return Series is
  begin
    return c*Inverse(s);
  end "/";

  function "/" ( c : Ring.number;
                 s : Link_to_Series ) return Link_to_Series is
  begin
    if s = null
     then return null;
     else return new Series'(c/Inverse(s.all));
    end if;
  end "/";

  procedure Div ( s : in out Series; c : in Ring.number ) is
  begin
    for k in 0..s.deg loop
      s.cff(k) := s.cff(k)/c;
    end loop;
  end Div;

  procedure Div ( s : in out Link_to_Series; c : in Ring.number ) is
  begin
    if s /= null
     then Div(s.all,c);
    end if;
  end Div;

  function "/" ( s,t : Series ) return Series is
  begin
    return s*Inverse(t);
  end "/";

  function "/" ( s,t : Link_to_Series ) return Link_to_Series is
  begin
    if s = null or t = null
     then return null;
     else return new Series'(s.all/t.all);
    end if;
  end "/";

  procedure Div ( s : in out Series; t : in Series ) is

    invt : constant Series := Inverse(t);

  begin
    Mul(s,invt);
  end Div;

  procedure Div ( s : in out Link_to_Series; t : in Link_to_Series ) is
  begin
    if s /= null and t /= null then
      declare
        invt : constant Series := Inverse(t.all);
        res : constant Series := s.all*invt;
      begin
        if res.deg = s.deg then
          s.cff := res.cff;
        else
          Clear(s);
          s := new Series'(res);
        end if;
      end;
    end if;
  end Div;

  function "**" ( s : Series; p : integer ) return Series is

    res : Series(s.deg);

  begin
    if p = 0 then
      res := Create(1,s.deg); -- keep the same degree as s
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

  function "**" ( s : Link_to_Series;
                  p : integer ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s = null then
      if p = 0
       then res := Create(1); -- 0**0 = 1
       else res := null;
      end if;
    else
      res := new Series'(s.all**p);
    end if;
    return res;
  end "**";

  procedure Power ( result : in out Link_to_Series;
                    s : in Link_to_Series; p : in integer ) is
  begin
    if p = 0 then
      result.cff(0) := Ring.one;
      for i in 1..result.deg loop
        result.cff(i) := Ring.zero;
      end loop;
    else
      for i in 0..result.deg loop
        result.cff(i) := s.cff(i);
      end loop;
      for i in 2..p loop
        Mul(result,s);
      end loop;
    end if;
  end Power;

-- DESTRUCTORS :

  procedure Clear ( s : in out Series ) is
  begin
    for i in s.cff'range loop
      s.cff(i) := Ring.zero;
    end loop;
  end Clear;

  procedure Clear ( s : in out Link_to_Series ) is

    procedure free is new unchecked_deallocation(Series,Link_to_Series);
    
  begin
    if s /= null
     then free(s);
    end if;
  end Clear;

end Generic_Dense_Series;
