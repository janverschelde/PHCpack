with text_io;                            use text_io;
with Standard_Mathematical_Functions;
with Double_Double_Basics;               use Double_Double_Basics;

package body Double_Double_Numbers is

-- CONSTRUCTORS :

  function Create ( i : integer ) return double_double is

    d : double_double;

  begin
    d.hi := double_float(i);
    d.lo := 0.0;
    return d;
  end Create;

  function Create ( n : natural32 ) return double_double is

    d : double_double;

  begin
    d.hi := double_float(n);
    d.lo := 0.0;
    return d;
  end Create;

  function Create ( n : natural64 ) return double_double is

    d : double_double;

  begin
    d.hi := double_float(n);
    d.lo := 0.0;
    return d;
  end Create;

  function Create ( i : integer32 ) return double_double is

    d : double_double;

  begin
    d.hi := double_float(i);
    d.lo := 0.0;
    return d;
  end Create;

  function Create ( i : integer64 ) return double_double is

    d : double_double;

  begin
    d.hi := double_float(i);
    d.lo := 0.0;
    return d;
  end Create;

  function Create ( f : double_float ) return double_double is

    d : double_double;

  begin
    d.hi := f;
    d.lo := 0.0;
    return d;
  end Create;

  function Create ( hi,lo : double_float ) return double_double is

    d : double_double;

  begin
    d.hi := hi;
    d.lo := lo;
    return d;
  end Create;

-- SELECTORS :

  function lo_part ( d : Double_Double ) return double_float is
  begin
    return d.lo;
  end lo_part;

  function hi_part ( d : Double_Double ) return double_float is
  begin
    return d.hi;
  end hi_part;

-- COMPARISON and COPYING :

  function is_zero ( d : double_double ) return boolean is
  begin
    return ((d.hi = 0.0) and (d.lo = 0.0));
  end is_zero;

  function is_one ( d : double_double ) return boolean is
  begin
    return ((d.hi = 1.0) and (d.lo = 0.0));
  end is_one;

  function is_positive ( d : double_double ) return boolean is
  begin
    return (d.hi > 0.0);
  end is_positive;

  function is_negative ( d : double_double ) return boolean is
  begin
    return (d.hi < 0.0);
  end is_negative;

  function equal ( x,y : double_double ) return boolean is
  begin
    return (x.hi = y.hi and x.lo = y.lo);
  end equal;

  function equal ( x : double_double; y : double_float ) return boolean is
  begin
    return (x.hi = y and x.lo = 0.0);
  end equal;

  function "<" ( x,y : double_double ) return boolean is
  begin
    return (x.hi < y.hi) or (x.hi = y.hi and x.lo < y.lo);
  end "<";

  function "<" ( x : double_double; y : double_float ) return boolean is
  begin
    return (x.hi < y) or (x.hi = y and x.lo < 0.0);
  end "<";

  function "<" ( x : double_float; y : double_double ) return boolean is
  begin
    return (x < y.hi) or (x = y.hi and y.lo > 0.0);
  end "<";

  function "<=" ( x,y : double_double ) return boolean is
  begin
    return (x.hi < y.hi) or (x.hi = y.hi and x.lo <= y.lo);
  end "<=";

  function "<=" ( x : double_double; y : double_float ) return boolean is
  begin
    return (x.hi < y) or (x.hi = y and x.lo <= 0.0);
  end "<=";

  function "<=" ( x : double_float; y : double_double ) return boolean is
  begin
    return (x < y.hi) or (x = y.hi and y.lo >= 0.0);
  end "<=";

  function ">" ( x,y : double_double ) return boolean is
  begin
    return (x.hi > y.hi) or (x.hi = y.hi and x.lo > y.lo);
  end ">";

  function ">" ( x : double_double; y : double_float ) return boolean is
  begin
    return (x.hi > y) or (x.hi = y and x.lo > 0.0);
  end ">";

  function ">" ( x : double_float; y : double_double ) return boolean is
  begin
    return (x > y.hi) or (x = y.hi and y.lo < 0.0);
  end ">";

  function ">=" ( x,y : double_double ) return boolean is
  begin
    return (x.hi > y.hi) or (x.hi = y.hi and x.lo >= y.lo);
  end ">=";

  function ">=" ( x : double_double; y : double_float ) return boolean is
  begin
    return (x.hi > y) or (x.hi = y and x.lo >= 0.0);
  end ">=";

  function ">=" ( x : double_float; y : double_double ) return boolean is
  begin
    return (x > y.hi) or (x = y.hi and y.lo <= 0.0);
  end ">=";

  procedure copy ( x : in double_double; y : in out double_double ) is
  begin
    y.hi := x.hi;
    y.lo := x.lo;
  end copy;

-- Absolute value and type casts :

  function to_int ( x : double_double ) return integer32 is
  begin
    return integer32(x.hi);
  end to_int;

  function to_double ( x : double_double ) return double_float is
  begin
    return x.hi;
  end to_double;

  function "abs" ( x : double_double ) return double_double is

    res : double_double;

  begin
    if x.hi < 0.0
     then res.hi := -x.hi; res.lo := -x.lo;
     else res.hi :=  x.hi; res.lo :=  x.lo;
    end if;
    return res;
  end "abs";

  function AbsVal ( x : double_double ) return double_double is

    res : double_double;

  begin
    if x.hi < 0.0
     then res.hi := -x.hi; res.lo := -x.lo;
     else res.hi :=  x.hi; res.lo :=  x.lo;
    end if;
    return res;
  end AbsVal;

  function floor ( x : double_double ) return double_double is

    res : double_double;
    hi : constant double_float := double_float'floor(x.hi);
    lo : double_float := 0.0;

  begin
    if(hi = x.hi) then   -- high word is already integer, round low word
      lo := double_float'floor(x.lo);
      quick_two_sum(hi,lo,res.hi,res.lo);
    else
      res.hi := hi; res.lo := lo;
    end if;
    return res;
  end floor;

  function nint ( x : double_float ) return double_float is

    f : constant double_float := double_float'floor(x);

  begin
    if x = f then
      return x;
    else
      declare
        ff : constant double_float := double_float'floor(x+0.5);
      begin
        return ff;
      end;
    end if; 
  end nint;

  function nint ( x : double_double ) return double_double is

    f : constant double_double := floor(x);

  begin
    if x = f then
      return x;
    else
      declare
        ff : constant double_double := floor(x+0.5);
      begin
        return ff;
      end;
    end if; 
  end nint;

-- ARITHMETICAL OPERATIONS :

  function "+" ( x,y : double_double ) return double_double is

    res : double_double;
    s1,s2,t1,t2 : double_float;

  begin
    two_sum(x.hi,y.hi,s1,s2);
    two_sum(x.lo,y.lo,t1,t2);
    s2 := s2 + t1;
    quick_two_sum(s1,s2,s1,s2);
    s2 := s2 + t2;
    quick_two_sum(s1,s2,res.hi,res.lo);
    return res;
  end "+";

  function "+" ( x : double_double; y : double_float ) return double_double is

    res : double_double;
    s1,s2 : double_float;

  begin
    two_sum(x.hi,y,s1,s2);
    s2 := s2 + x.lo;
    quick_two_sum(s1,s2,res.hi,res.lo);
    return res;
  end "+";

  function "+" ( x : double_float; y : double_double ) return double_double is

    res : double_double;
    s1,s2 : double_float;

  begin
    two_sum(x,y.hi,s1,s2);
    s2 := s2 + y.lo;
    quick_two_sum(s1,s2,res.hi,res.lo);
    return res;
  end "+";

  function "+" ( x,y : double_float ) return double_double is

    res : double_double;

  begin
    two_sum(x,y,res.hi,res.lo);
    return res;
  end "+";

  function "+" ( x : double_double ) return double_double is

    res : double_double;

  begin
    res.hi := x.hi;
    res.lo := x.lo;
    return res;
  end "+";

  procedure Add ( x : in out double_double; y : in double_double ) is

    s1,s2,t1,t2 : double_float;

  begin
    two_sum(x.hi,y.hi,s1,s2);
    two_sum(x.lo,y.lo,t1,t2);
    s2 := s2 + t1;
    quick_two_sum(s1,s2,s1,s2);
    s2 := s2 + t2;
    quick_two_sum(s1,s2,x.hi,x.lo);
  end Add;

  procedure Add ( x : in out double_double; y : in double_float ) is
 
    s1,s2 : double_float;

  begin
    two_sum(x.hi,y,s1,s2);
    s2 := s2 + x.lo;
    quick_two_sum(s1,s2,x.hi,x.lo);
  end Add;

  function "-" ( x,y : double_double ) return double_double is

    res : double_double;
    s1,s2,t1,t2 : double_float;

  begin
    two_diff(x.hi,y.hi,s1,s2);
    two_diff(x.lo,y.lo,t1,t2);
    s2 := s2 + t1;
    quick_two_sum(s1,s2,s1,s2);
    s2 := s2 + t2;
    quick_two_sum(s1,s2,res.hi,res.lo);
    return res;
  end "-";

  function "-" ( x : double_double; y : double_float ) return double_double is

    res : double_double;
    s1,s2 : double_float;

  begin
    two_diff(x.hi,y,s1,s2);
    s2 := s2 + x.lo;
    quick_two_sum(s1,s2,res.hi,res.lo);
    return res;
  end "-";

  function "-" ( x : double_float; y : double_double ) return double_double is

    res : double_double;
    s1,s2 : double_float;

  begin
    two_diff(x,y.hi,s1,s2);
    s2 := s2 - y.lo;
    quick_two_sum(s1,s2,res.hi,res.lo);
    return res;
  end "-";

  function "-" ( x,y : double_float ) return double_double is

    res : double_double;

  begin
    two_diff(x,y,res.hi,res.lo);
    return res;
  end "-";

  function "-" ( x : double_double ) return double_double is

    res : double_double;

  begin
    res.hi := -x.hi;
    res.lo := -x.lo;
    return res;
  end "-";

  procedure Min ( x : in out double_double ) is
  begin
    x.hi := -x.hi;
    x.lo := -x.lo;
  end Min;

  procedure Sub ( x : in out double_double; y : in double_double ) is

    s1,s2,t1,t2 : double_float;

  begin
    two_diff(x.hi,y.hi,s1,s2);
    two_diff(x.lo,y.lo,t1,t2);
    s2 := s2 + t1;
    quick_two_sum(s1,s2,s1,s2);
    s2 := s2 + t2;
    quick_two_sum(s1,s2,x.hi,x.lo);
  end Sub;

  procedure Sub ( x : in out double_double; y : in double_float ) is

    s1,s2 : double_float;

  begin
    two_diff(x.hi,y,s1,s2);
    s2 := s2 + x.lo;
    quick_two_sum(s1,s2,x.hi,x.lo);
  end Sub;

  function "*" ( x,y : double_double ) return double_double is

    res : double_double;
    p1,p2 : double_float;

  begin
    two_prod(x.hi,y.hi,p1,p2);
    p2 := p2 + (x.hi * y.lo + x.lo * y.hi);
    quick_two_sum(p1,p2,res.hi,res.lo);
    return res;
  end "*";

  function "*" ( x : double_double; y : double_float ) return double_double is

    res : double_double;
    p1,p2 : double_float;

  begin
    two_prod(x.hi,y,p1,p2);
    p2 := p2 + (x.lo * y);
    quick_two_sum(p1,p2,res.hi,res.lo);
    return res;
  end "*";

  function "*" ( x : double_float; y : double_double ) return double_double is

    res : double_double;
    p1,p2 : double_float;

  begin
    two_prod(x,y.hi,p1,p2);
    p2 := p2 + (y.lo * x);
    quick_two_sum(p1,p2,res.hi,res.lo);
    return res;
  end "*";

  function "*" ( x,y : double_float ) return double_double is

    res : double_double;

  begin
    two_prod(x,y,res.hi,res.lo);
    return res;
  end "*";

  procedure Mul ( x : in out double_double; y : in double_double ) is

    p1,p2 : double_float;

  begin
    two_prod(x.hi,y.hi,p1,p2);
    p2 := p2 + y.lo * x.hi;
    p2 := p2 + y.hi * x.lo;
    quick_two_sum(p1,p2,x.hi,x.lo);
  end Mul;

  procedure Mul ( x : in out double_double; y : in double_float ) is

    p1,p2 : double_float;

  begin
    two_prod(x.hi,y,p1,p2);
    p2 := p2 + x.lo * y;
    quick_two_sum(p1,p2,x.hi,x.lo);
  end Mul;

  function Mul_pwr2 ( x : double_double; y : double_float )
                    return double_double is

    res : double_double;

  begin
    res.hi := x.hi*y;
    res.lo := x.lo*y;
    return res;
  end Mul_pwr2;

  procedure Mul_pwr2 ( x : in out double_double; y : in double_float ) is
  begin
    x.hi := x.hi*y;
    x.lo := x.lo*y;
  end Mul_pwr2;

  function "/" ( x,y : double_double ) return double_double is

    res,acc : double_double;
    q1,q2,q3 : double_float;

  begin
    q1 := x.hi/y.hi;        -- approximate quotient
    acc := q1*y;
    res := x - acc;         -- res = x - q1 * b
    q2 := res.hi/y.hi;
    acc := q2*y;           
    Sub(res,acc);           -- res -= (q2 * b)
    q3 := res.hi/y.hi;
    quick_two_sum(q1,q2,res.hi,res.lo);
    Add(res,q3);            -- res = dd_real(q1,q2) + q3
    return res;
  end "/";

  function "/" ( x : double_double; y : double_float ) return double_double is

    res : double_double;
    q1,q2,p1,p2,s,e : double_float;

  begin
    q1 := x.hi/y;                          -- approximate quotient
    two_prod(q1,y,p1,p2);                  -- compute x - q1 * y 
    two_diff(x.hi,p1,s,e);
    e := e + x.lo;
    e := e - p2;
    q2 := (s + e)/y;                       -- get next approximation
    quick_two_sum(q1,q2,res.hi,res.lo);    -- renormalize
    return res;
  end "/";

  function "/" ( x : double_float; y : double_double ) return double_double is

    xx : double_double;

  begin
    xx.hi := x; 
    xx.lo := 0.0;
    return xx/y;
  end "/";

  function "/" ( x,y : double_float ) return double_double is

    res : double_double;
    q1,q2,p1,p2,s,e : double_float;

  begin
    q1 := x/y;
    two_prod(q1,y,p1,p2);    -- compute x - q1 * y
    two_diff(x,p1,s,e);
    e := e - p2;
    q2 := (s + e)/y;         -- get next approximation
    quick_two_sum(q1,q2,res.hi,res.lo);
    return res;
  end "/";

  procedure Div ( x : in out double_double; y : in double_double ) is

    acc : double_double;
    q1,q2,q3 : double_float;

  begin
    q1 := x.hi/y.hi;        -- approximate quotient
    acc := q1*y;
    Sub(x,acc);             -- x -= q1 * b
    q2 := x.hi/y.hi;
    acc := q2*y;
    Sub(x,acc);             -- x -= (q2 * b)
    q3 := x.hi/y.hi;
    quick_two_sum(q1,q2,x.hi,x.lo);
    Add(x,q3);              -- x = dd_real(q1,q2) + q3
  end Div;

  procedure Div ( x : in out double_double; y : in double_float ) is

    q1,q2,p1,p2,s,e : double_float;

  begin
    q1 := x.hi/y;                          -- approximate quotient
    two_prod(q1,y,p1,p2);                  -- compute x - q1 * y 
    two_diff(x.hi,p1,s,e);
    e := e + x.lo;
    e := e - p2;
    q2 := (s + e)/y;                       -- get next approximation
    quick_two_sum(q1,q2,x.hi,x.lo);        -- renormalize
  end Div;

  function sqr ( x : double_float ) return double_double is

    res : double_double;

  begin
    two_sqr(x,res.hi,res.lo);
    return res;
  end sqr;
 
  function sqr ( x : double_double ) return double_double is

    res : double_double;
    p1,p2 : double_float;
 
  begin
    two_sqr(x.hi,p1,p2);
    p2 := p2 + 2.0 * x.hi * x.lo;
    p2 := p2 + x.lo * x.lo;
    quick_two_sum(p1,p2,res.hi,res.lo);
    return res;
  end sqr;

  function "**" ( x : double_double; n : integer ) return double_double is

    res,acc : double_double;
    absn : natural;

  begin
    if n = 0 then
      res.hi := 1.0; res.lo := 0.0;
    else
      if n > 0
       then absn := n;
       else absn := -n;
      end if;
      res.hi := x.hi; res.lo := x.lo;
      acc.hi := 1.0;  acc.lo := 0.0;
      if absn > 1 then          -- use binary exponentiation
        while absn > 0 loop
          if absn mod 2 = 1
           then Mul(acc,res);
          end if;
          absn := absn/2;
          if absn > 0
           then res := sqr(res);
          end if;
        end loop;
      else
        acc.hi := res.hi; acc.lo := res.lo;
      end if;
      if n < 0
       then res := 1.0/acc;          -- compute reciprocal
       else res.hi := acc.hi; res.lo := acc.lo;
      end if;
    end if;
    return res;
  end "**";

  function "**" ( x : double_double; n : integer32 ) return double_double is
  begin
    return x**integer(n);
  end "**";

  function "**" ( x : double_double; n : integer64 ) return double_double is

    res,acc : double_double;
    absn : natural64;

  begin
    if n = 0 then
      res.hi := 1.0; res.lo := 0.0;
    else
      if n > 0
       then absn := natural64(n);
       else absn := natural64(-n);
      end if;
      res.hi := x.hi; res.lo := x.lo;
      acc.hi := 1.0;  acc.lo := 0.0;
      if absn > 1 then          -- use binary exponentiation
        while absn > 0 loop
          if absn mod 2 = 1
           then Mul(acc,res);
          end if;
          absn := absn/2;
          if absn > 0
           then res := sqr(res);
          end if;
        end loop;
      else
        acc.hi := res.hi; acc.lo := res.lo;
      end if;
      if n < 0
       then res := 1.0/acc;          -- compute reciprocal
       else res.hi := acc.hi; res.lo := acc.lo;
      end if;
    end if;
    return res;
  end "**";

  function ldexp ( x : double_double; n : integer ) return double_double is

    res : double_double;
    function C_ldexp ( x : double_float; n : integer ) return double_float;
    pragma Import(C,C_ldexp,External_Name => "ldexp");

  begin
    res.hi := C_ldexp(x.hi,n);
    res.lo := C_ldexp(x.lo,n);
    return res;
  end ldexp;

  function "**" ( x,y : double_double ) return double_double is
  begin
    return exp(y*log(x));
  end "**";

  function "**" ( x : double_double; y : double_float ) return double_double is

    dd_y : constant double_double := create(y);

  begin
    return x**dd_y;
  end "**";

  function exp ( x : double_double ) return double_double is

  -- Strategy:  We first reduce the size of x by noting that
  --        exp(kr + m * log(2)) = 2^m * exp(r)^k
  --   where m and k are integers.  By choosing m appropriately
  --   we can make |kr| <= log(2) / 2 = 0.347.  
  --   Then exp(r) is evaluated using the familiar Taylor series.  
  --   Reducing the argument substantially speeds up the convergence.

    res : double_double;
    k : constant double_float := 512.0;
    inv_k : constant double_float := 1.0/k;
    dd_eps : constant double_float := 4.93038065763132e-32;     -- 2^-104
    tol : constant double_float := inv_k*dd_eps;
    exp1hi : constant double_float := 2.718281828459045091e+00; -- exp(1)
    exp1lo : constant double_float := 1.445646891729250158e-16;
    exp1 : constant double_double := Create(exp1hi,exp1lo);
    log2hi : constant double_float := 6.931471805599452862e-01; -- log(2)
    log2lo : constant double_float := 2.319046813846299558e-17;
    log2 : constant double_double := Create(log2hi,log2lo);
    m : constant double_float := double_float'floor(x.hi/log2hi + 0.5);
    i_fac : array(0..14) of double_double;
      -- inverse factorials for Taylor expansion
    p,s,t : double_double;
    cnt : integer;

  begin
    if x.hi <= -709.0 then
      res.hi := 0.0; res.lo := 0.0;
    elsif x.hi >= 709.0 then
      res.hi := -1.0; res.lo := 0.0;
    elsif is_zero(x) then
      res.hi := 1.0; res.lo := 0.0;
    elsif is_one(x) then
      res := exp1;
    else
      i_fac(0)  := Create(1.66666666666666657E-01, 9.25185853854297066E-18);
      i_fac(1)  := Create(4.16666666666666644E-02, 2.31296463463574266E-18);
      i_fac(2)  := Create(8.33333333333333322E-03, 1.15648231731787138E-19);
      i_fac(3)  := Create(1.38888888888888894E-03,-5.30054395437357706E-20);
      i_fac(4)  := Create(1.98412698412698413E-04, 1.72095582934207053E-22);
      i_fac(5)  := Create(2.48015873015873016E-05, 2.15119478667758816E-23);
      i_fac(6)  := Create(2.75573192239858925E-06,-1.85839327404647208E-22);
      i_fac(7)  := Create(2.75573192239858883E-07, 2.37677146222502973E-23);
      i_fac(8)  := Create(2.50521083854417202E-08,-1.44881407093591197E-24);
      i_fac(9)  := Create(2.08767569878681002E-09,-1.20734505911325997E-25);
      i_fac(10) := Create(1.60590438368216133E-10, 1.25852945887520981E-26);
      i_fac(11) := Create(1.14707455977297245E-11, 2.06555127528307454E-28);
      i_fac(12) := Create(7.64716373181981641E-13, 7.03872877733453001E-30);
      i_fac(13) := Create(4.77947733238738525E-14, 4.39920548583408126E-31);
      i_fac(14) := Create(2.81145725434552060E-15, 1.65088427308614326E-31);
      res := mul_pwr2(x - m*log2,inv_k);
      p := sqr(res);
      s := res + mul_pwr2(p,0.5);
      p := p*res;
      t := p*i_fac(0);
      cnt := 0;
      loop
        s := s + t;
        p := p*res;
        t := i_fac(cnt+1)*p;
        exit when abs(t.hi) <= tol;
        cnt := cnt + 1;
        exit when (cnt >= 5);
      end loop;
      s := s + t;
      for i in 1..9 loop
        s := mul_pwr2(s,2.0) + sqr(s);
      end loop;
      res := s + 1.0;
      cnt := integer(m);
      res := ldexp(res,cnt);
    end if;
    return res;
  end exp;

  function log ( x : double_double ) return double_double is

  -- Strategy.  The Taylor series for log converges much more
  --   slowly than that of exp, due to the lack of the factorial
  --   term in the denominator.  Hence this routine instead tries
  --   to determine the root of the function f(x) = exp(x) - a
  --   using Newton iteration.  The iteration is given by
  --       x' = x - f(x)/f'(x)
  --          = x - (1 - a * exp(-x))
  --          = x + a * exp(-x) - 1.
  --   Only one iteration is needed, since Newton's iteration
  --   approximately doubles the number of digits per iteration.

   -- function C_ln ( x : double_float ) return double_float;
   -- pragma Import(C,C_ln,External_Name => "log");

    res,acc : double_double;

  begin
    if is_one(x) then
      res.hi := 0.0; res.lo := 0.0;
    elsif x.hi <= 0.0 then
      put_line("dd_log: argument is not positive");
      res.hi := -1.0; res.lo := 0.0; 
    else
      res.hi := Standard_Mathematical_Functions.LN(x.hi);
     -- res.hi := C_ln(x.hi);  -- initial approximation
      res.lo := 0.0;       
      acc := x*exp(-res);    -- y = x + a*exp(-x) - 1
      res := res + acc;
      res := res - 1.0;
    end if;
    return res;
  end log;

  function log10 ( x : double_double ) return double_double is

    res : double_double;
    log10hi : constant double_float := 2.302585092994045901e+00; -- high word
    log10lo : constant double_float := -2.170756223382249351e-16; -- low word
    logten : constant double_double := Create(log10hi,log10lo);
 
  begin
    res := log(x)/logten;
    return res;
  end log10;

-- DESTRUCTOR :

  procedure clear ( d : in out double_double ) is
  begin
    d.hi := 0.0; d.lo := 0.0;
  end clear;

end Double_Double_Numbers;
