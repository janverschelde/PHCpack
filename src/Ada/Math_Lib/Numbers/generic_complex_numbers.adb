package body Generic_Complex_Numbers is

--  TWO : constant number := one+one;

-- CREATORS :

  function Create ( n : natural32 ) return Complex_Number is
  begin
    return Create(integer32(n));
  end Create;

  function Create ( i : integer ) return Complex_Number is
  begin
    return Create(integer32(i));
  end Create;

  function Create ( i : integer32 ) return Complex_Number is

    res : Complex_Number;

  begin
    if i = 0 then
      res.re := +zero;
    elsif i = 1 then
      res.re := +one;
    else 
      res.re := Create(integer(i));
    end if;
    res.im := +zero;
    return res;
  end Create;

  function Create ( f : number ) return Complex_Number is

    res : Complex_Number;

  begin
    res.re := +f;
    res.im := +zero;
    return res;
  end Create;

  function Create ( re,im : number ) return Complex_Number is

    res : Complex_Number;

  begin
    res.RE := +re;
    res.IM := +im;
    return res;
  end Create;

  function Conjugate ( c : Complex_Number ) return Complex_Number is

    res : Complex_Number;

  begin
    res.RE := +c.RE;
    res.IM := -c.IM;
    return res;
  end Conjugate;

-- COMPARISON and COPYING :

  function Equal ( x,y : Complex_Number ) return boolean is
  begin
    return (Equal(x.RE,y.RE) and Equal(x.IM,y.IM));
  end Equal;

  procedure Copy ( x : in Complex_Number; y : in out Complex_Number ) is
  begin
    Copy(x.RE,y.RE);
    Copy(x.IM,y.IM);
  end Copy;

  function "<" ( x,y : Complex_Number ) return boolean is

    avx : number := AbsVal(x);
    avy : number := AbsVal(y);
    res : constant boolean := (avx < avy);

  begin
    Clear(avx);
    Clear(avy);
    return res;
  end "<";

  function ">" ( x,y : Complex_Number ) return boolean is

    avx : number := AbsVal(x);
    avy : number := AbsVal(y);
    res : constant boolean := (avx > avy);

  begin
    Clear(avx);
    Clear(avy);
    return res;
  end ">";

-- SELECTORS :

  function REAL_PART ( x : Complex_Number ) return number is
  begin
    return +x.RE;
  end REAL_PART;

  function IMAG_PART ( x : Complex_Number ) return number is
  begin
    return +x.IM;
  end IMAG_PART;

  function AbsVal ( x : Complex_Number ) return number is

    res : number := AbsVal(x.RE);
    acc : number := AbsVal(x.IM);

  begin
    Add(res,acc);
    Clear(acc);
    return res;
  end AbsVal;

  function AbsVal ( x : Complex_Number ) return Complex_Number is

    abx : number := AbsVal(x);
    res : constant Complex_Number := Create(abx);

  begin
    Clear(abx);
    return res;
  end AbsVal;

-- ARITHMETIC OPERATIONS AS FUNCTIONS :

  function "+" ( x : Complex_Number; y : number ) return Complex_Number is
  begin
    return (x.RE+y,+x.IM);
  end "+";

  function "-" ( x : Complex_Number; y : number ) return Complex_Number is
  begin
    return (x.RE-y,+x.IM);
  end "-";

  function "*" ( x : Complex_Number; y : number ) return Complex_Number is
  begin
    return (x.RE*y,x.IM*y);
  end "*";

  function "/" ( x : Complex_Number; y : number ) return Complex_Number is
  begin
    return (x.RE/y,x.IM/y);
  end "/";

  function "+" ( x : number; y : Complex_Number ) return Complex_Number is
  begin
    return (x+y.RE,+y.IM);
  end "+";

  function "-" ( x : number; y : Complex_Number ) return Complex_Number is
  begin
    return (x-y.RE,-y.IM);
  end "-";

  function "*" ( x : number; y : Complex_Number ) return Complex_Number is
  begin
    return (x*y.RE,x*y.IM);
  end "*";

  function "/" ( x : number; y : Complex_Number ) return Complex_Number is

    res : Complex_Number;
    acc,nrm : number;

  begin
    nrm := y.RE*y.RE;
    acc := y.IM*y.IM;
    Add(nrm,acc);
    Clear(acc);
    acc := x/nrm;
    res.RE := acc*y.RE;
    res.IM := acc*y.IM; Min(res.IM);
    Clear(nrm); Clear(acc);
    return res;
  end "/";

--  function "/" ( x : number; y : Complex_Number ) return Complex_Number is
--
--    res : Complex_Number;
--    acc,avyim,avyre : Number;
--
--  begin
--    if Equal(y.IM,zero) then
--      res.RE := x/y.RE;
--      res.IM := +zero;
--    elsif Equal(y.RE,zero) then
--      res.RE := +zero;
--      res.IM := x/y.IM; Min(res.IM);
--    else
--      avyim := AbsVal(y.IM); avyre := AbsVal(y.RE);
--      if avyim < avyre then
--        acc := y.IM/y.RE;
--        res.RE := x*acc; Mul(acc,y.IM); Add(acc,y.RE); Div(res.RE,acc);
--        res.IM := x/acc; Min(res.IM);
--        Clear(acc);
--      elsif avyim > avyre then
--        acc := y.RE/y.IM;
--        res.IM := x*acc; Min(res.IM); Mul(acc,y.RE); Add(acc,y.IM);
--        res.RE := x/acc;
--        Clear(acc);
--      elsif Equal(y.IM,y.RE) then
--        acc := TWO*y.IM;
--        res.RE := x/acc;
--        res.IM := -res.RE;
--        Clear(acc);
--      else -- y.IM = -y.RE then
--        acc := TWO*y.IM;
--        res.RE := x/acc; Min(res.RE);
--        res.IM := x/acc; Min(res.IM);
--        Clear(acc);
--      end if;
--      Clear(avyim); Clear(avyre);
--    end if;
--    return res;
--  end "/";

  function "+" ( x : Complex_Number ) return Complex_Number is
  begin
    return (+x.RE,+x.IM);
  end "+";

  function "-" ( x : Complex_Number ) return Complex_Number is
  begin
    return (-x.RE,-x.IM);
  end "-";

  function "+" ( x,y : Complex_Number ) return Complex_Number is
  begin
    return (x.RE+y.RE,x.IM+y.IM);
  end "+";

  function "-" ( x,y : Complex_Number ) return Complex_Number is
  begin
    return (x.RE-y.RE,x.IM-y.IM);
  end "-";

  function "*" ( x,y : Complex_Number ) return Complex_Number is

    res : Complex_Number;
    acc : number;

  begin
    acc := x.IM*y.IM;
    res.RE := x.RE*y.RE;
    Sub(res.RE,acc);
    Clear(acc);
    acc := x.IM*y.RE;
    res.IM := x.RE*y.IM;
    Add(res.IM,acc);
    Clear(acc);
    return res;
  end "*";

--  function "*" ( x,y : Complex_Number ) return Complex_Number is
--
--    res : Complex_Number;
--    acc,avyim,avyre : number;
--
--  begin
--    if Equal(y.IM,zero) then
--      res.RE := x.RE*y.RE;
--      res.IM := x.IM*y.RE;
--    elsif Equal(y.RE,zero) then
--      res.RE := x.IM*y.IM; Min(res.RE); 
--      res.IM := x.RE*y.IM;
--    else
--      avyre := AbsVal(y.RE); avyim := AbsVal(y.IM);
--      if avyre < avyim then
--        acc := y.RE/y.IM;
--        res.RE := x.RE*acc; Sub(res.RE,x.IM); Mul(res.RE,y.IM);
--        res.IM := x.IM*acc; Add(res.IM,x.RE); Mul(res.IM,y.IM);
--        Clear(acc);
--      elsif avyre > avyim then
--        acc := y.IM/y.RE;
--        res.RE := x.IM*acc; Sub(res.RE,x.RE); Min(res.RE); Mul(res.RE,y.RE);
--        res.IM := x.RE*acc; Add(res.IM,x.IM); Mul(res.IM,y.RE);
--        Clear(acc);
--      elsif Equal(y.RE,y.IM) then
--        res.RE := x.RE - x.IM; Mul(res.RE,y.RE);
--        res.IM := x.IM + x.RE; Mul(res.IM,y.RE);
--      else -- y.RE = -y.IM then
--        res.RE := x.RE + x.IM; Mul(res.RE,y.RE);
--        res.IM := x.IM - x.RE; Mul(res.IM,y.RE); 
--      end if;
--      Clear(avyre); Clear(avyim);
--    end if;
--    return res;
--  end "*";

  function "/"  ( x,y : Complex_Number ) return Complex_Number is

    res : Complex_Number;
    acc,nrm : number;

  begin
    nrm := y.RE*y.RE;
    acc := y.IM*y.IM;
    Add(nrm,acc);
    Clear(acc);
    res.RE := x.RE*y.RE;
    acc := x.IM*y.IM;
    Add(res.RE,acc);
    Clear(acc);
    res.IM := x.IM*y.RE;
    acc := x.RE*y.IM;
    Sub(res.IM,acc);
    Clear(acc);
    Div(res.RE,nrm);
    Div(res.IM,nrm);
    Clear(nrm);
    return res;
  end "/";

--  function "/"  ( x,y : Complex_Number ) return Complex_Number is
--
--    res : Complex_Number;
--    acc,avyre,avyim : number;
--
--  begin
--    if Equal(y.IM,zero) then
--      res.RE := x.RE/y.RE;
--      res.IM := x.IM/y.RE;
--    elsif Equal(y.RE,zero) then
--      res.RE := x.IM/y.IM;
--      res.IM := x.RE/y.IM; Min(res.IM);
--    else
--      avyre := AbsVal(y.RE); avyim := AbsVal(y.IM);
--      if avyre < avyim then
--        acc := y.RE/y.IM;
--        res.RE := x.RE*acc; Add(res.RE,x.IM);
--        res.IM := x.IM*acc; Sub(res.IM,x.RE);
--        Mul(acc,y.RE); Add(acc,y.IM);
--        Div(res.RE,acc);
--        Div(res.IM,acc);
--        Clear(acc);
--      elsif avyre > avyim then
--        acc := y.IM/y.RE;
--        res.RE := x.IM*acc; Add(res.RE,x.RE);
--        res.IM := x.RE*acc; Sub(res.IM,x.IM); Min(res.IM);
--        Mul(acc,y.IM); Add(acc,y.RE);
--        Div(res.RE,acc);
--        Div(res.IM,acc);
--        Clear(acc);
--      elsif Equal(y.RE,y.IM) then
--        acc := TWO*y.RE;
--        res.RE := x.RE + x.IM; Div(res.RE,acc);
--        res.IM := x.IM - x.RE; Div(res.IM,acc);
--        Clear(acc);
--      else -- y.RE = -y.IM then
--        acc := TWO*y.RE;
--        res.RE := x.RE - x.IM; Div(res.RE,acc);
--        res.IM := x.IM + x.RE; Div(res.IM,acc);
--        Clear(acc);
--      end if;
--      Clear(avyre); Clear(avyim);
--    end if;
--    return res;
--  end "/";

  function "**" ( x : Complex_Number; m : integer ) return Complex_Number is

    res : Complex_Number;

  begin
    if m = 0 then
      res := Create(natural32(1));
    elsif m > 0 then
      res := +x; --Copy(x,res); <- Copy CAUSED BUS ERROR!!
      for j in 2..m loop
        Mul(res,x);
      end loop;
    else
      res := Create(natural32(1));
      for j in 1..(-m) loop
        Div(res,x);
      end loop;
    end if;
    return res;
  end "**";

-- ARITHMETIC OPERATIONS AS PROCEDURES :

  procedure Add ( x : in out Complex_Number; y : in number ) is
  begin
    Add(x.RE,y);
  end Add;

  procedure Sub ( x : in out Complex_Number; y : in number ) is
  begin
    Sub(x.RE,y);
  end Sub;

  procedure Mul ( x : in out Complex_Number; y : in number ) is
  begin
    Mul(x.RE,y);
    Mul(x.IM,y);
  end Mul;

  procedure Div ( x : in out Complex_Number; y : in number ) is
  begin
    Div(x.RE,y);
    Div(x.IM,y);
  end Div;

  procedure Add ( x : in out Complex_Number; y : in Complex_Number ) is
  begin
    Add(x.RE,y.RE);
    Add(x.IM,y.IM);
  end Add;

  procedure Sub ( x : in out Complex_Number; y : in Complex_Number ) is
  begin
    Sub(x.RE,y.RE);
    Sub(x.IM,y.IM);
  end Sub;

  procedure Min ( x : in out Complex_Number ) is
  begin
    Min(x.RE);
    Min(x.IM);
  end Min;

  procedure Mul ( x : in out Complex_Number; y : in Complex_Number ) is

    res : Complex_Number;
    acc : number;

  begin
    acc := x.IM*y.IM;
    res.RE := x.RE*y.RE;
    Sub(res.RE,acc);
    Clear(acc);
    acc := x.IM*y.RE;
    res.IM := x.RE*y.IM;
    Add(res.IM,acc);
    Clear(acc);
    Clear(x);
    x := res;   
  end Mul;

--  procedure Mul ( x : in out Complex_Number; y : in Complex_Number ) is
--
--    res : Complex_Number;
--    acc,avyim,avyre : number;
--
--  begin
--    if Equal(y.IM,zero) then
--      Mul(x.RE,y.RE);
--      Mul(x.IM,y.RE);
--    elsif Equal(y.RE,zero) then
--      res.RE := x.IM*y.IM; Min(res.RE);
--      res.IM := x.RE*y.IM;
--      Clear(x); x := res;
--    else
--      avyre := AbsVal(y.RE); avyim := AbsVal(y.IM);
--      if avyre < avyim then
--        acc := y.RE/y.IM;
--        res.RE := x.RE*acc; Sub(res.RE,x.IM); Mul(res.RE,y.IM);
--        res.IM := x.IM*acc; Add(res.IM,x.RE); Mul(res.IM,y.IM);
--        Clear(acc);
--      elsif avyre > avyim then
--        acc := y.IM/y.RE;
--        res.RE := x.IM*acc; Sub(res.RE,x.RE); Min(res.RE); Mul(res.RE,y.RE);
--        res.IM := x.RE*acc; Add(res.IM,x.IM); Mul(res.IM,y.RE);
--        Clear(acc);
--      elsif Equal(y.RE,y.IM) then
--        res.RE := x.RE - x.IM; Mul(res.RE,y.RE);
--        res.IM := x.IM + x.RE; Mul(res.IM,y.RE);
--      else -- y.RE = -y.IM then
--        res.RE := x.RE + x.IM; Mul(res.RE,y.RE);
--        res.IM := x.IM - x.RE; Mul(res.IM,y.RE);
--      end if;
--      Clear(avyre); Clear(avyim);
--      Clear(x); x := res;
--    end if;
--  end Mul;

  procedure Div ( x : in out Complex_Number; y : in Complex_Number ) is

    res : Complex_Number;
    acc,nrm : number;

  begin
    nrm := y.RE*y.RE;
    acc := y.IM*y.IM;
    Add(nrm,acc);
    Clear(acc);
    res.RE := x.RE*y.RE;
    acc := x.IM*y.IM;
    Add(res.RE,acc);
    Clear(acc);
    res.IM := x.IM*y.RE;
    acc := x.RE*y.IM;
    Sub(res.IM,acc);
    Clear(acc);
    Div(res.RE,nrm);
    Div(res.IM,nrm);
    Clear(nrm);
    Clear(x);
    x := res;
  end Div;

--  procedure Div ( x : in out Complex_Number; y : in Complex_Number ) is
--
--    res : Complex_Number;
--    acc,avyre,avyim : number;
--
--  begin
--    if Equal(y.IM,zero) then
--      Div(x.RE,y.RE);
--      Div(x.IM,y.RE);
--    elsif Equal(y.IM,zero) then
--      res.RE := x.IM/y.IM;
--      res.IM := x.RE/y.IM; Min(res.IM);
--      Clear(x); x := res;
--    else
--      avyre := AbsVal(y.RE); avyim := AbsVal(y.IM);
--      if avyre < avyim then
--        acc := y.RE/y.IM;
--        res.RE := x.RE*acc; Add(res.RE,x.IM);
--        res.IM := x.IM*acc; Sub(res.IM,x.RE);
--        Mul(acc,y.RE); Add(acc,y.IM);
--        Div(res.RE,acc);
--        Div(res.IM,acc);
--        Clear(acc);
--      elsif avyre > avyim then
--        acc := y.IM/y.RE;
--        res.RE := x.IM*acc; Add(res.RE,x.RE);
--        res.IM := x.RE*acc; Sub(res.IM,x.IM); Min(res.IM);
--        Mul(acc,y.IM); Add(acc,y.RE);
--        Div(res.RE,acc);
--        Div(res.IM,acc);
--        Clear(acc);
--      elsif Equal(y.RE,y.IM) then
--        acc := TWO*y.RE;
--        res.RE := x.RE + x.IM; Div(res.RE,acc);
--        res.IM := x.IM - x.RE; Div(res.IM,acc);
--        Clear(acc);
--      else -- y.RE = -y.IM then
--        acc := TWO*y.RE;
--        res.RE := x.RE - x.IM; Div(res.RE,acc);
--        res.IM := x.IM + x.RE; Div(res.IM,acc);
--        Clear(acc);
--      end if;
--      Clear(avyre); Clear(avyim);
--      Clear(x); x := res;
--    end if;
--  end Div;

-- DESTRUCTOR :

  procedure Clear ( x : in out Complex_Number ) is
  begin
    Clear(x.RE);
    Clear(x.IM);
  end Clear;

end Generic_Complex_Numbers;
