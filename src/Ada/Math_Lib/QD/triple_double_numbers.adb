with text_io;                            use text_io;
with Standard_Mathematical_Functions;
with Double_Double_Basics;
with Double_Double_Numbers;
with Quad_Double_Renormalizations;
with Fast_Double_Renormalizations;       use Fast_Double_Renormalizations;

package body Triple_Double_Numbers is

-- CONSTRUCTORS :

  function create ( i : integer ) return triple_double is

    res : triple_double;

  begin
    res.hi := double_float(i);
    res.mi := 0.0;
    res.lo := 0.0;
    return res;
  end create;

  function create ( n : natural32 ) return triple_double is

    res : triple_double;

  begin
    res.hi := double_float(n);
    res.mi := 0.0;
    res.lo := 0.0;
    return res;
  end create;

  function create ( n : natural64 ) return triple_double is

    res : triple_double;

  begin
    res.hi := double_float(n);
    res.mi := 0.0;
    res.lo := 0.0;
    return res;
  end create;

  function create ( i : integer32 ) return triple_double is

    res : triple_double;

  begin
    res.hi := double_float(i);
    res.mi := 0.0;
    res.lo := 0.0;
    return res;
  end create;

  function create ( i : integer64 ) return triple_double is

    res : triple_double;

  begin
    res.hi := double_float(i);
    res.mi := 0.0;
    res.lo := 0.0;
    return res;
  end create;

  function create ( x : double_float ) return triple_double is

    res : triple_double;

  begin
    res.hi := x;
    res.mi := 0.0;
    res.lo := 0.0;
    return res;
  end create;

  function create ( hi,mi,lo : double_float ) return triple_double is

    res : triple_double;

  begin
    res.hi := hi;
    res.mi := mi;
    res.lo := lo;
    return res;
  end create;

  function "abs" ( x : triple_double ) return triple_double is

    res : triple_double;

  begin
    if x.hi < 0.0
     then res.hi := -x.hi; res.mi := -x.mi; res.lo := -x.lo;
     else res.hi :=  x.hi; res.mi :=  x.mi; res.lo :=  x.lo;
    end if;
    return res;
  end "abs";

  function AbsVal ( x : triple_double ) return triple_double is

    res : triple_double;

  begin
    if x.hi < 0.0
     then res.hi := -x.hi; res.mi := -x.mi; res.lo := -x.lo;
     else res.hi :=  x.hi; res.mi :=  x.mi; res.lo :=  x.lo;
    end if;
    return res;
  end AbsVal;

  function floor ( x : triple_double ) return triple_double is

    res : triple_double;
    x0,x1,x2,x3 : double_float;

  begin
    x1 := 0.0; x2 := 0.0; x3 := 0.0;
    x0 := double_float'floor(x.hi);
    if x0 = x.hi then
      x1 := double_float'floor(x.mi);
      if x1 = x.mi then
        x2 := double_float'floor(x.lo);
      end if;
      Quad_Double_Renormalizations.renorm4(x0,x1,x2,x3);
    end if;
    res := Create(x0,x1,x2);
    return res;
  end floor;

  function nint ( x : triple_double ) return triple_double is

    res : triple_double;
    x0,x1,x2,x3 : double_float;

  begin
    x1 := 0.0; x2 := 0.0; x3 := 0.0;
    x0 := Double_Double_Numbers.nint(x.hi);
    if x0 = x.hi then      -- first double is already an integer
      x1 := Double_Double_Numbers.nint(x.mi);
      if x1 = x.mi then    -- second double is already an integer
        x2 := Double_Double_Numbers.nint(x.lo);
      else
        if abs(x1 - x.mi) = 0.5 and x.lo < 0.0
         then x1 := x1 - 1.0;
        end if;
      end if;
    else                     -- first double is not an integer
      if abs(x0 - x.hi) = 0.5 and x.mi < 0.0
       then x0 := x0 - 1.0;
      end if;
    end if;
    Quad_Double_Renormalizations.renorm4(x0,x1,x2,x3);
    res := Create(x0,x1,x2);
    return res;
  end nint;

-- SELECTORS :

  function hi_part ( x : triple_double ) return double_float is
  begin
    return x.hi;
  end hi_part;

  function mi_part ( x : triple_double ) return double_float is
  begin
    return x.mi;
  end mi_part;

  function lo_part ( x : triple_double ) return double_float is
  begin
    return x.lo;
  end lo_part;

-- TYPE CASTS :

  function to_int ( x : triple_double ) return integer32 is
  begin
    return integer32(x.hi);
  end to_int;

  function to_double ( x : triple_double ) return double_float is
  begin
    return x.hi;
  end to_double;

-- COMPARISON and COPYING :

  function is_zero ( x : triple_double ) return boolean is
  begin
    return ((x.hi = 0.0) and (x.mi = 0.0) and (x.lo = 0.0));
  end is_zero;

  function is_one ( x : triple_double ) return boolean is
  begin
    return ((x.hi = 1.0) and (x.mi = 0.0) and (x.lo = 0.0));
  end is_one;

  function is_positive ( x : triple_double ) return boolean is
  begin
    return (x.hi > 0.0);
  end is_positive;

  function is_negative ( x : triple_double ) return boolean is
  begin
    return (x.hi < 0.0);
  end is_negative;

  function equal ( x,y : triple_double ) return boolean is
  begin
    return ((x.hi = y.hi) and (x.mi = y.mi) and (x.lo = y.lo));
  end equal;

  function equal ( x : triple_double; y : double_float ) return boolean is
  begin
    return ((x.hi = y) and (x.mi = 0.0) and (x.lo = 0.0));
  end equal;

  function "<" ( x,y : triple_double ) return boolean is
  begin
    return ((x.hi < y.hi)
         or (x.hi = y.hi and x.mi < y.mi)
         or (x.mi = y.mi and x.mi = y.mi and x.lo < y.lo));
  end "<";

  function "<" ( x : triple_double; y : double_float ) return boolean is
  begin
    return ((x.hi < y) or (x.hi = y and x.mi < 0.0)
                       or (x.hi = y and x.mi = 0.0 and x.lo < 0.0));
  end "<";

  function "<" ( x : double_float; y : triple_double ) return boolean is
  begin
    return (y > x);
  end "<";

  function "<=" ( x,y : triple_double ) return boolean is
  begin
    return x < y or equal(x,y);
  end "<=";

  function "<=" ( x : triple_double; y : double_float ) return boolean is
  begin
    return ((x.hi < y) or (x.hi = y and x.mi < 0.0)
                       or (x.hi = y and x.mi = 0.0 and x.lo <= 0.0));
  end "<=";

  function "<=" ( x : double_float; y : triple_double ) return boolean is
  begin
    return (y >= x);
  end "<=";

  function ">" ( x,y : triple_double ) return boolean is
  begin
    return ((x.hi > y.hi)
         or (x.hi = y.hi and x.mi > y.mi)
         or (x.hi = y.hi and x.mi = y.mi and x.lo > y.lo));
  end ">";

  function ">" ( x : triple_double; y : double_float ) return boolean is
  begin
    return ((x.hi > y) or (x.hi = y and x.mi > 0.0));
  end ">";

  function ">" ( x : double_float; y : triple_double ) return boolean is
  begin
    return (y < x);
  end ">";

  function ">=" ( x,y : triple_double ) return boolean is
  begin
    return x > y or equal(x,y);
  end ">=";

  function ">=" ( x : triple_double; y : double_float ) return boolean is
  begin
    return ((x.hi > y) or (x.hi = y and x.mi > 0.0)
                       or (x.hi = y and x.mi = 0.0 and x.lo >= 0.0));
  end ">=";

  function ">=" ( x : double_float; y : triple_double ) return boolean is
  begin
    return (y <= x);
  end ">=";

  procedure copy ( x : in triple_double; y : in out triple_double ) is
  begin
    y.hi := x.hi; y.mi := x.mi; y.lo := x.lo;
  end copy;

-- ARITHMETICAL OPERATIONS :

  function "+" ( x,y : triple_double ) return triple_double is

    res : triple_double;
    f0,f1,f2,f3,e : double_float;

  begin
    f3 := 0.0;
    Double_Double_Basics.two_sum(x.lo,y.lo,f2,e);
    f3 := f3 + e;
    Double_Double_Basics.two_sum(x.mi,y.mi,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    f3 := f3 + e;
    Double_Double_Basics.two_sum(x.hi,y.hi,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    f3 := f3 + e;
    fast_renorm(f0,f1,f2,f3,res.hi,res.mi,res.lo);
    return res;
  end "+";

  function "+" ( x : triple_double; y : double_float ) return triple_double is

    res : triple_double;

  begin
    renorm_add1(x.hi,x.mi,x.lo,y,res.hi,res.mi,res.lo);
    return res;
  end "+";

  function "+" ( x : double_float; y : triple_double ) return triple_double is

    res : constant triple_double := y + x;

  begin
    return res;
  end "+";

  function "+" ( x : triple_double ) return triple_double is

    res : triple_double;

  begin
    res.hi := x.hi; res.mi := x.mi; res.lo := x.lo;
    return res;
  end "+";

  function "-" ( x : triple_double ) return triple_double is

    res : triple_double;

  begin
    res.hi := -x.hi;
    res.mi := -x.mi;
    res.lo := -x.lo;
    return res;
  end "-";

  function "-" ( x,y : triple_double ) return triple_double is

    mny : constant triple_double := -y;
    res : constant triple_double := x + mny;

  begin
    return res;
  end "-";

  function "-" ( x : triple_double; y : double_float ) return triple_double is

    mny : constant double_float := -y;
    res : constant triple_double := x + mny;

  begin
    return res;
  end "-";

  function "-" ( x : double_float; y : triple_double ) return triple_double is

    mny : constant triple_double := -y;
    res : constant triple_double := x + mny;

  begin
    return res;
  end "-";

  function "*" ( x,y : triple_double ) return triple_double is

    res : triple_double;
    f0,f1,f2,f3,p,e : double_float;

  begin
    f3 := x.mi*y.lo;
    f3 := f3 + x.lo*y.mi;
    Double_Double_Basics.two_prod(x.hi,y.lo,f2,e);
    f3 := f3 + e;
    Double_Double_Basics.two_prod(x.mi,y.mi,p,e);
    f3 := f3 + e;
    Double_Double_Basics.two_sum(f2,p,f2,e);
    f3 := f3 + e;
    Double_Double_Basics.two_prod(x.lo,y.hi,p,e);
    f3 := f3 + e;
    Double_Double_Basics.two_sum(f2,p,f2,e);
    f3 := f3 + e;
    Double_Double_Basics.two_prod(x.hi,y.mi,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    f3 := f3 + e;
    Double_Double_Basics.two_prod(x.mi,y.hi,p,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    f3 := f3 + e;
    Double_Double_Basics.two_sum(f1,p,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    f3 := f3 + e;
    Double_Double_Basics.two_prod(x.hi,y.hi,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    f3 := f3 + e;
    fast_renorm(f0,f1,f2,f3,res.hi,res.mi,res.lo);
    return res;
  end "*";

  function "*" ( x : triple_double; y : double_float ) return triple_double is

    res : triple_double;
    f0,f1,f2,f3,e : double_float;

  begin
    f3 := 0.0;
    Double_Double_Basics.two_prod(x.lo,y,f2,e);
    f3 := f3 + e;
    Double_Double_Basics.two_prod(x.mi,y,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    f3 := f3 + e;
    Double_Double_Basics.two_prod(x.hi,y,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    f3 := f3 + e;
    fast_renorm(f0,f1,f2,f3,res.hi,res.mi,res.lo);
    return res;
  end "*";

  function "*" ( x : double_float; y : triple_double ) return triple_double is
  begin
    return y*x;
  end "*";

  function Mul_pwr2 ( x : triple_double; y : double_float )
                    return triple_double is

    res : triple_double;

  begin
    res.hi := x.hi*y;
    res.mi := x.mi*y;
    res.lo := x.lo*y;
    return res;
  end Mul_pwr2;

  procedure Mul_pwr2 ( x : in out triple_double; y : in double_float ) is
  begin
    x.hi := x.hi*y;
    x.mi := x.mi*y;
    x.lo := x.lo*y;
  end Mul_pwr2;

  function "/" ( x,y : triple_double ) return triple_double is

    res,acc : triple_double;
    q0,q1,q2,q3 : double_float;

  begin
    q0 := x.hi/y.hi;   acc := q0*y; res := x - acc;
    q1 := res.hi/y.hi; acc := q1*y; res := res - acc;
    q2 := res.hi/y.hi; acc := q2*y; res := res - acc;
    q3 := res.hi/y.hi;
    fast_renorm(q0,q1,q2,q3,res.hi,res.mi,res.lo);
    return res;
  end "/";

  function "/" ( x : double_float; y : triple_double ) return triple_double is

    tdx : constant triple_double := create(x);
    res : constant triple_double := tdx/y;

  begin
    return res;
  end "/";

  function "/" ( x : triple_double; y : double_float ) return triple_double is

    tdy : constant triple_double := create(y);
    res : constant triple_double := x/tdy;

  begin
    return res;
  end "/";

  function sqr ( x : triple_double ) return triple_double is
  begin
    return x*x;
  end sqr;

  function "**" ( x : triple_double; n : integer ) return triple_double is

    res,acc : triple_double;
    absn : natural;

  begin
    if n = 0 then
      res.hi := 1.0; res.mi := 0.0; res.lo := 0.0;
    else
      if n > 0
       then absn := n;
       else absn := -n;
      end if;
      res.hi := x.hi; res.mi := x.mi; res.lo := x.lo;
      acc.hi := 1.0;  acc.mi := 0.0;  acc.lo := 0.0;
      if absn > 1 then          -- use binary exponentiation
        while absn > 0 loop
          if absn mod 2 = 1
           then acc := acc*res;
          end if;
          absn := absn/2;
          if absn > 0
           then res := res*res;
          end if;
        end loop;
      else
        acc.hi := res.hi; acc.mi := res.mi; acc.lo := res.lo;
      end if;
      if n < 0
       then res := 1.0/acc;          -- compute reciprocal
       else res.hi := acc.hi; res.mi := acc.mi; res.lo := acc.lo;
      end if;
    end if;
    return res;
  end "**";

  function ldexp ( x : triple_double; n : integer ) return triple_double is

    res : triple_double;
    function C_ldexp ( x : double_float; n : integer ) return double_float;
    pragma Import(C,C_ldexp,External_Name => "ldexp");

  begin
    res.hi := C_ldexp(x.hi,n);
    res.mi := C_ldexp(x.mi,n);
    res.lo := C_ldexp(x.lo,n);
    return res;
  end ldexp;

  function "**" ( x,y : triple_double ) return triple_double is
  begin
    return exp(y*log(x));
  end "**";

  function "**" ( x : triple_double; y : double_float ) return triple_double is

    td_y : constant triple_double := create(y);

  begin
    return x**td_y;
  end "**";

  function exp ( x : triple_double ) return triple_double is

  -- Strategy:  We first reduce the size of x by noting that
  --    exp(kr + m * log(2)) = 2^m * exp(r)^k
  -- where m and k are integers.  By choosing m appropriately
  -- we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
  -- evaluated using the familiar Taylor series.  Reducing the
  -- argument substantially speeds up the convergence.

    function C_ldexp ( x : double_float; n : integer ) return double_float;
    pragma Import(C,C_ldexp,External_Name => "ldexp");

    res : triple_double;
    k : constant double_float := C_ldexp(1.0,16);
    inv_k : constant double_float := 1.0/k;
    e_0 : constant double_float :=  2.718281828459045091e+00; -- exp(1)[0]
    e_1 : constant double_float :=  1.445646891729250158e-16; -- exp(1)[1]
    e_2 : constant double_float := -2.127717108038176765e-33; -- exp(1)[2]
    exp1 : constant triple_double := Create(e_0,e_1,e_2);
    log2_0 : constant double_float :=  6.931471805599452862e-01; -- log(2)[0]
    log2_1 : constant double_float :=  2.319046813846299558e-17; -- log(2)[1]
    log2_2 : constant double_float :=  5.707708438416212066e-34; -- log(2)[2]
    log2 : constant triple_double := Create(log2_0,log2_1,log2_2);
    td_eps : constant double_float := 5.473822126268817e-48; -- 2.0**(-157)
    tol : constant double_float := inv_k*td_eps;
    m : constant double_float := double_float'floor(x.hi/log2_0 + 0.5);
    i_fac : array(0..14) of triple_double;
      -- inverse factorials for Taylor expansion
    p,s,t : triple_double;
    cnt : integer;

  begin
    if x.hi <= -709.0 then
      res := Create(0.0);
    elsif x.hi >= 709.0 then
      res := Create(-1.0);
    elsif is_zero(x) then
      res := Create(1.0);
    elsif is_one(x) then
      res := exp1;
    else
      i_fac(0)  := Create( 1.66666666666666657e-01, 9.25185853854297066e-18,
                           5.13581318503262866e-34);
      i_fac(1)  := Create( 4.16666666666666644e-02, 2.31296463463574266e-18,
                           1.28395329625815716e-34);
      i_fac(2)  := Create( 8.33333333333333322e-03, 1.15648231731787138e-19,
                           1.60494162032269652e-36);
      i_fac(3)  := Create( 1.38888888888888894e-03,-5.30054395437357706e-20,
                          -1.73868675534958776e-36);
      i_fac(4)  := Create( 1.98412698412698413e-04, 1.72095582934207053e-22,
                           1.49269123913941271e-40);
      i_fac(5)  := Create( 2.48015873015873016e-05, 2.15119478667758816e-23,
                           1.86586404892426588e-41);
      i_fac(6)  := Create( 2.75573192239858925e-06,-1.85839327404647208e-22,
                           8.49175460488199287e-39);
      i_fac(7)  := Create( 2.75573192239858883e-07, 2.37677146222502973e-23,
                          -3.26318890334088294e-40);
      i_fac(8)  := Create( 2.50521083854417202e-08,-1.44881407093591197e-24,
                           2.04267351467144546e-41);
      i_fac(9)  := Create( 2.08767569878681002e-09,-1.20734505911325997e-25,
                           1.70222792889287100e-42);
      i_fac(10) := Create( 1.60590438368216133e-10, 1.25852945887520981e-26,
                          -5.31334602762985031e-43);
      i_fac(11) := Create( 1.14707455977297245e-11, 2.06555127528307454e-28,
                           6.88907923246664603e-45);
      i_fac(12) := Create( 7.64716373181981641e-13, 7.03872877733453001e-30,
                          -7.82753927716258345e-48);
      i_fac(13) := Create( 4.77947733238738525e-14, 4.39920548583408126e-31,
                          -4.89221204822661465e-49);
      i_fac(14) := Create( 2.81145725434552060e-15, 1.65088427308614326e-31,
                          -2.87777179307447918e-50);
      res := mul_pwr2(x - m*log2,inv_k);
      p := res*res;
      s := res + mul_pwr2(p,0.5);
      cnt := 0;
      loop
        p := p*res;
        t := i_fac(cnt)*p;
        s := s + t;
        exit when abs(t.hi) <= tol;
        cnt := cnt + 1;
        exit when (cnt >= 9);
      end loop;
      for i in 1..16 loop -- 16 times s = mul_pwr2(s,2.0) + sqr(s);
        p := Mul_pwr2(s,2.0);
        t := s*s;
        s := p + t;
      end loop;
      res := s + 1.0;
      cnt := integer(m);
      res := ldexp(res,cnt);
    end if;
    return res;
  end exp;

  function log ( x : triple_double ) return triple_double is

  -- Strategy.  The Taylor series for log converges much more
  --   slowly than that of exp, due to the lack of the factorial
  --   term in the denominator.  Hence this routine instead tries
  --   to determine the root of the function f(x) = exp(x) - a
  --   using Newton iteration.  The iteration is given by
  --       x' = x - f(x)/f'(x)
  --          = x - (1 - a * exp(-x))
  --          = x + a * exp(-x) - 1.
  --   Three iterations are needed, since Newton's iteration
  --   approximately doubles the number of digits per iteration.

    res,acc : triple_double;

  begin
    if is_one(x) then
      res := Create(0.0);
    elsif x.hi <= 0.0 then
      put_line("td_log: argument is not positive");
      res := Create(-1.0);
    else
      res := Create(Standard_Mathematical_Functions.LN(x.hi));
      for i in 1..3 loop       -- r = r + x*exp(-r) - 1
        acc := x*exp(-res);    -- acc = x*exp(-res)      
        res := res + acc;      -- res = res + x*exp(-res)
        res := res - 1.0;
      end loop;
    end if;
    return res;
  end log;

  function log10 ( x : triple_double ) return triple_double is

    res : triple_double;
    log10_0 : constant double_float :=  2.302585092994045901e+00;
    log10_1 : constant double_float := -2.170756223382249351e-16;
    log10_2 : constant double_float := -9.984262454465776570e-33;
    logten : constant triple_double := Create(log10_0,log10_1,log10_2);

  begin
    res := log(x)/logten;
    return res;
  end log10;

-- ARITHMETICAL OPERATIONS AS PROCEDURES :

  procedure Add ( x : in out triple_double; y : in triple_double ) is
  begin
    x := x + y;
  end Add;

  procedure Sub ( x : in out triple_double; y : in triple_double ) is
  begin
    x := x - y;
  end Sub;

  procedure Min ( x : in out triple_double ) is
  begin
    x := -x;
  end Min;

  procedure Mul ( x : in out triple_double; y : in triple_double ) is
  begin
    x := x*y;
  end Mul;

  procedure Div ( x : in out triple_double; y : in triple_double ) is
  begin
    x := x/y;
  end Div;

-- DESTRUCTOR :

  procedure clear ( x : in out triple_double ) is
  begin
    x.hi := 0.0; x.mi := 0.0; x.lo := 0.0;
  end Clear;

end Triple_Double_Numbers;
