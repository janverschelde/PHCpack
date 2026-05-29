with unchecked_deallocation;
with Double_rpSeries_Operations;

package body Double_Real_Power_Series is

-- CONSTRUCTORS :

  function Make ( x : integer32; tdx : integer32 ) return Series is

    res : Series(tdx);

  begin
    res.cff(0) := create(double_float(x));
    for i in 1..tdx loop
      res.cff(i) := create(0.0);
      res.pwt(i) := 0.0;
    end loop;
    return res;
  end Make;

  function Make ( x : integer64; tdx : integer32 ) return Series is

    res : Series(tdx);

  begin
    res.cff(0) := create(double_float(x));
    for i in 1..tdx loop
      res.cff(i) := create(0.0);
      res.pwt(i) := 0.0;
    end loop;
    return res;
  end Make;

  function Make ( x : double_float; tdx : integer32 ) return Series is

    res : Series(tdx);

  begin
    res.cff(0) := create(x);
    for i in 1..tdx loop
      res.cff(i) := create(0.0);
      res.pwt(i) := 0.0;
    end loop;
    return res;
  end Make;

  function Make ( x : complex_number; tdx : integer32 ) return Series is

    res : Series(tdx);

  begin
    res.cff(0) := x;
    for i in 1..tdx loop
      res.cff(i) := create(0.0);
      res.pwt(i) := 0.0;
    end loop;
    return res;
  end Make;

  function Make ( x : integer32; tdx : integer32 ) return Link_to_Series is

    rep : constant Series(tdx) := Make(x,tdx);
    res : constant Link_to_Series := new Series'(rep);

  begin
    return res;
  end Make;

  function Make ( x : integer64; tdx : integer32 ) return Link_to_Series is

    rep : constant Series(tdx) := Make(x,tdx);
    res : constant Link_to_Series := new Series'(rep);

  begin
    return res;
  end Make;

  function Make ( x : double_float; tdx : integer32 ) return Link_to_Series is

    rep : constant Series(tdx) := Make(x,tdx);
    res : constant Link_to_Series := new Series'(rep);

  begin
    return res;
  end Make;

  function Make ( x : complex_number;
                  tdx : integer32 ) return Link_to_Series is

    rep : constant Series(tdx) := Make(x,tdx);
    res : constant Link_to_Series := new Series'(rep);

  begin
    return res;
  end Make;

  function Make ( cff : Standard_Complex_Vectors.Vector;
                  pwt : Standard_Floating_Vectors.Vector ) return Series is

    res : Series(cff'last);

  begin
    res.cff(0) := cff(0);
    for i in 1..cff'last loop
      res.cff(i) := cff(i);
      res.pwt(i) := pwt(i);
    end loop;
    return res;
  end Make;

  function Make ( cff : Standard_Complex_Vectors.Vector;
                  pwt : Standard_Floating_Vectors.Vector )
                return Link_to_Series is

    rep : constant Series(cff'last) := Make(cff,pwt);
    res : constant Link_to_Series := new Series'(rep);

  begin
    return res;
  end Make;

-- EQUALITY AND COPY :

  function Equal ( a,b : Series;
                   tol : double_float := 1.0E-12 ) return boolean is
  begin
    if a.tdx > b.tdx then
      return Equal(a=>b,b=>a,tol=>tol);
    elsif AbsVal(a.cff(0) - b.cff(0)) > tol then
      return false;
    else
      for i in a.pwt'range loop                -- a.tdx <= b.tdx
        if AbsVal(a.cff(i) - b.cff(i)) > tol
         then return false;
        end if;
        if abs(a.pwt(i) - b.pwt(i)) > tol
         then return false;
        end if;
      end loop;
      for i in a.tdx+1..b.tdx loop 
        if AbsVal(b.cff(i)) > tol    -- are extra coefficients zero?
         then return false;
        end if;
        if abs(b.pwt(i)) > tol       -- are extra powers zero?
         then return false;
        end if;
      end loop;
      return true;
    end if;
  end Equal;

  function Equal ( a,b : Link_to_Series;
                   tol : double_float := 1.0E-12 ) return boolean is
  begin
    if a = null then
      if b = null
       then return true;
       else return false;
      end if;
    elsif b = null then
      return false;
    else
      return Equal(a.all,b.all,tol);
    end if;
  end Equal;

  procedure Copy ( a : in Series; b : in out Series ) is
  begin
    b.cff(0) := a.cff(0);
    for i in a.pwt'range loop
      exit when (i > b.pwt'last);
      b.cff(i) := a.cff(i);
      b.pwt(i) := a.pwt(i);
    end loop;
  end Copy;

  procedure Copy ( a : in Link_to_Series; b : in out Link_to_Series ) is
  begin
    Clear(b);
    b := Make(a.cff,b.pwt);
  end Copy;

-- ARITHMETICAL OPERATORS :

  function "+" ( s : Series; c : complex_number ) return Series is

    res : Series := s;

  begin
    res.cff(0) := s.cff(0) + c;
    return res;
  end "+";

  function "+" ( s : Link_to_Series;
                 c : complex_number ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s = null then
      res := Make(c,0);
    else
      res := Make(s.cff,s.pwt);
      res.cff(0) := res.cff(0) + c;
    end if;
    return res;
  end "+";

  function "+" ( c : complex_number; s : Series ) return Series is

    res : Series := s;

  begin
    res.cff(0) := c + s.cff(0);
    return res;
  end "+";

  function "+" ( c : complex_number;
                 s : Link_to_Series ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s = null then
      res := Make(c,0);
    else
      res := Make(s.cff,s.pwt);
      res.cff(0) := c + res.cff(0);
    end if;
    return res;
  end "+";

  procedure Add ( s : in out Series; c : in complex_number ) is
  begin
    s.cff(0) := s.cff(0) + c;
  end Add;

  procedure Add ( s : in out Link_to_Series;
                  c : in complex_number ) is
  begin
    if s = null
     then s := Make(c,0);
     else s.cff(0) := s.cff(0) + c;
    end if;
  end Add;

  function "+" ( s : Series ) return Series is

    res : constant Series(s.tdx) := s;

  begin
    return res;
  end "+";

  function "+" ( s : Link_to_Series ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s /= null
     then Copy(s,res);
    end if;
    return res;
  end "+";

  function "+" ( a,b : Series ) return Series is

    res : Series(a.tdx + b.tdx);

  begin
    Double_rpSeries_Operations.Add(a.cff,b.cff,a.pwt,b.pwt,res.cff,res.pwt);
    return res;
  end "+";

  function "+" ( a,b : Link_to_Series ) return Link_to_Series is

    rep : constant Series(a.tdx + b.tdx) := a.all + b.all;
    res : constant Link_to_Series := new Series'(rep);

  begin
    return res;
  end "+";

  procedure Add ( a : in out Series; b : in Series ) is

    res : Series(a.tdx + b.tdx);

  begin
    Double_rpSeries_Operations.Add(a.cff,b.cff,a.pwt,b.pwt,res.cff,res.pwt);
    a.cff(0) := res.cff(0);
    for i in 1..a.tdx loop
      a.cff(i) := res.cff(i);
      a.pwt(i) := res.pwt(i);
    end loop;
  end Add;

  procedure Add ( a : in out Link_to_Series; b : in Link_to_Series ) is
  begin
    if a /= null and b /= null then
      declare
        res : Series(a.tdx + b.tdx);
      begin
        Double_rpSeries_Operations.Add
          (a.cff,b.cff,a.pwt,b.pwt,res.cff,res.pwt);
        a.cff(0) := res.cff(0);
        for i in 1..a.tdx loop
          a.cff(i) := res.cff(i);
          a.pwt(i) := res.pwt(i);
        end loop;
      end;
    end if;
  end Add;

  function "-" ( s : Series; c : complex_number ) return Series is

    res : Series(s.tdx) := s;

  begin
    res.cff(0) := s.cff(0) - c;
    return res;
  end "-";

  function "-" ( s : Link_to_Series;
                 c : complex_number ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s /= null then
      Copy(s,res);
      res.cff(0) := s.cff(0) - c;
    end if;
    return res;
  end "-";

  function "-" ( c : complex_number; s : Series ) return Series is

    res : Series(s.tdx) := s;

  begin
    res.cff(0) := c - s.cff(0);
    for i in 1..res.cff'last loop
      res.cff(i) := -res.cff(i);
    end loop;
    return res;
  end "-";

  function "-" ( c : complex_number;
                 s : Link_to_Series ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s = null then
      res := Make(c,0);
    else
      res := Make(s.cff,s.pwt);
      res.cff(0) := c - res.cff(0);
      for i in 1..res.cff'last loop
        res.cff(i) := -res.cff(i);
      end loop;
    end if;
    return res;
  end "-";

  procedure Sub ( s : in out Series; c : in complex_number ) is
  begin
    s.cff(0) := s.cff(0) - c;
  end Sub;

  procedure Sub ( s : in out Link_to_Series;
                  c : in complex_number ) is
  begin
    if s = null
     then s := Make(-c,0);
     else s.cff(0) := s.cff(0) - c;
    end if;
  end Sub;

  function "-" ( s : Series ) return Series is

    res : Series(s.tdx) := s;

  begin
    for i in res.cff'range loop
      res.cff(i) := -res.cff(i);
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
    for i in s.cff'range loop
      s.cff(i) := -s.cff(i);
    end loop;
  end Min;

  procedure Min ( s : in out Link_to_Series ) is
  begin
    if s /= null then
      for i in s.cff'range loop
        s.cff(i) := -s.cff(i);
      end loop;
    end if;
  end Min;

  function "-" ( a,b : Series ) return Series is

    res : Series(a.tdx + b.tdx);

  begin
    Double_rpSeries_Operations.Sub(a.cff,b.cff,a.pwt,b.pwt,res.cff,res.pwt);
    return res;
  end "-";

  function "-" ( a,b : Link_to_Series ) return Link_to_Series is

    rep : constant Series(a.tdx + b.tdx) := a.all - b.all;
    res : constant Link_to_Series := new Series'(rep);

  begin
    return res;
  end "-";

  procedure Sub ( a : in out Series; b : in Series ) is

    res : Series(a.tdx + b.tdx);

  begin
    Double_rpSeries_Operations.Sub(a.cff,b.cff,a.pwt,b.pwt,res.cff,res.pwt);
    a.cff(0) := res.cff(0);
    for i in 1..a.tdx loop
      a.cff(i) := res.cff(i);
      a.pwt(i) := res.pwt(i);
    end loop;
  end Sub;

  procedure Sub ( a : in out Link_to_Series; b : in Link_to_Series ) is


  begin
    if a /= null and b /= null then
      declare
        res : Series(a.tdx + b.tdx);
      begin
        Double_rpSeries_Operations.Sub
          (a.cff,b.cff,a.pwt,b.pwt,res.cff,res.pwt);
        a.cff(0) := res.cff(0);
        for i in 1..a.tdx loop
          a.cff(i) := res.cff(i);
          a.pwt(i) := res.pwt(i);
        end loop;
      end;
    end if;
  end Sub;

  function "*" ( s : Series; c : complex_number ) return Series is

    res : Series(s.tdx) := s;

  begin
    for i in res.cff'range loop
      res.cff(i) := res.cff(i)*c;
    end loop;
    return res;
  end "*";

  function "*" ( s : Link_to_Series;
                 c : complex_number ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s /= null then
      Copy(s,res);
      for i in res.cff'range loop
        res.cff(i) := res.cff(i)*c;
      end loop;
    end if;
    return res;
  end "*";

  function "*" ( c : complex_number; s : Series ) return Series is

    res : Series(s.tdx) := s;

  begin
    for i in res.cff'range loop
      res.cff(i) := c*res.cff(i);
    end loop;
    return res;
  end "*";

  function "*" ( c : complex_number;
                 s : Link_to_Series ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s /= null then
      Copy(s,res);
      for i in res.cff'range loop
        res.cff(i) := c*res.cff(i);
      end loop;
    end if;
    return res;
  end "*";

  procedure Mul ( s : in out Series; c : in complex_number ) is
  begin
    for i in s.cff'range loop
      s.cff(i) := s.cff(i)*c;
    end loop;
  end Mul;

  procedure Mul ( s : in out Link_to_Series; c : in complex_number ) is
  begin
    if s /= null then
      for i in s.cff'range loop
        s.cff(i) := s.cff(i)*c;
      end loop;
    end if;
  end Mul;

  function "*" ( a,b : Series ) return Series is

    size : constant integer32 := (a.tdx+1)*(b.tdx+1) - 1;
    res : Series(size);

  begin
    Double_rpSeries_Operations.Mul(a.cff,b.cff,a.pwt,b.pwt,res.cff,res.pwt);
    return res;
  end "*";

  function "*" ( a,b : Link_to_Series ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if a /= null and b /= null then
      declare
        size : constant integer32 := (a.tdx+1)*(b.tdx+1) - 1;
        rep : Series(size);
      begin
        Double_rpSeries_Operations.Mul
          (a.cff,b.cff,a.pwt,b.pwt,rep.cff,rep.pwt);
        res := new Series'(rep);
      end;
    end if;
    return res;
  end "*";
 
  procedure Mul ( a : in out Series; b : in Series ) is

    size : constant integer32 := (a.tdx+1)*(b.tdx+1) - 1;
    res : constant Series(size) := a*b;

  begin
    a.cff(0) := res.cff(0);
    for i in 1..a.tdx loop
      a.cff(i) := res.cff(i);
      a.pwt(i) := res.pwt(i);
    end loop;
  end Mul;

  procedure Mul ( a : in out Link_to_Series; b : in Link_to_Series ) is
  begin
    if b = null then
      Clear(a);
    elsif a /= null then
      declare
        size : constant integer32 := (a.tdx+1)*(b.tdx+1) - 1;
        res : constant Series(size) := a.all*b.all;
      begin
        a.cff(0) := res.cff(0);
        for i in 1..a.tdx loop
          a.cff(i) := res.cff(i);
          a.pwt(i) := res.pwt(i);
        end loop;
      end;
    end if;
  end Mul;

  function Inverse ( s : Series ) return Series is

    res : Series(s.tdx);

  begin
    Double_rpSeries_Operations.Inv(s.cff,s.pwt,res.cff,res.pwt);
    return res;
  end Inverse;

  function Inverse ( s : Link_to_Series ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s /= null then
      declare
        rep : constant Series(s.tdx) := Inverse(s.all);
      begin
        res := new Series'(rep);
      end;
    end if;
    return res;
  end Inverse;

  function "/" ( s : Series; c : complex_number ) return Series is

    res : Series(s.tdx) := s;

  begin
    for i in res.cff'range loop
      res.cff(i) := res.cff(i)/c;
    end loop;
    return res;
  end "/";

  function "/" ( s : Link_to_Series;
                 c : complex_number ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s /= null then
      declare
        rep : Series(s.tdx) := s.all;
      begin
        for i in rep.cff'range loop
          rep.cff(i) := rep.cff(i)/c;
        end loop;
        res := new Series'(rep);
      end;
    end if;
    return res;
  end "/";

  function "/" ( c : complex_number; s : Series ) return Series is

    res : Series(s.tdx) := Inverse(s);

  begin
    for i in res.cff'range loop
      res.cff(i) := res.cff(i)*c;
    end loop;
    return res;
  end "/";

  function "/" ( c : complex_number;
                 s : Link_to_Series ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if s /= null then
      declare
        rep : Series(s.tdx) := Inverse(s.all);
      begin
        for i in rep.cff'range loop
          rep.cff(i) := rep.cff(i)*c;
        end loop;
        res := new Series'(rep);
      end;
    end if;
    return res;
  end "/";

  procedure Div ( s : in out Series; c : in complex_number ) is
  begin
    for i in s.cff'range loop
      s.cff(i) := s.cff(i)/c;
    end loop;
  end Div;

  procedure Div ( s : in out Link_to_Series; c : in complex_number ) is
  begin
    if s /= null then
      for i in s.cff'range loop
        s.cff(i) := s.cff(i)/c;
      end loop;
    end if;
  end Div;

  function "/" ( a,b : Series ) return Series is

    size : constant integer32 := (a.tdx+1)*(b.tdx+1)-1;
    invb : constant Series(b.tdx) := Inverse(b);
    res : constant Series(size) := a*invb;

  begin
    return res;
  end "/";

  function "/" ( a,b : Link_to_Series ) return Link_to_Series is

    res : Link_to_Series;

  begin
    if a /= null and b /= null then
      declare
        size : constant integer32 := (a.tdx+1)*(b.tdx+1)-1;
        invb : constant Series(b.tdx) := Inverse(b.all);
        rep : constant Series(size) := a.all*invb;
      begin
        res := new Series'(rep);
      end;
    end if;
    return res;
  end "/";

  procedure Div ( a : in out Series; b : in Series ) is

    size : constant integer32 := (a.tdx+1)*(b.tdx+1)-1;
    res : constant Series(size) := a/b;

  begin
    a.cff(0) := res.cff(0);
    for i in 1..a.tdx loop
      a.cff(i) := res.cff(i);
      a.pwt(i) := res.pwt(i);
    end loop;
  end Div;

  procedure Div ( a : in out Link_to_Series; b : in Link_to_Series ) is
  begin
    if a /= null and b /= null then
      declare
        size : constant integer32 := (a.tdx+1)*(b.tdx+1)-1;
        res : constant Series(size) := a.all/b.all;
      begin
        a.cff(0) := res.cff(0);
        for i in 1..a.tdx loop
          a.cff(i) := res.cff(i);
          a.pwt(i) := res.pwt(i);
        end loop;
      end;
    end if;
  end Div;

-- DESTRUCTORS :

  procedure Clear ( s : in out Series ) is
  begin
    s.cff(0) := create(0.0);
    for i in 1..s.tdx loop
      s.cff(i) := create(0.0);
      s.pwt(i) := 0.0;
    end loop;
  end Clear;

  procedure Clear ( s : in out Link_to_Series ) is

    procedure free is new unchecked_deallocation(Series,Link_to_Series);
    
  begin
    if s /= null
     then free(s);
    end if;
  end Clear;

end Double_Real_Power_Series;
