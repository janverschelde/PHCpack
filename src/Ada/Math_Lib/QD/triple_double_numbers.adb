with Double_Double_Basics;
with Fast_Double_Renormalizations;       use Fast_Double_Renormalizations;

package body Triple_Double_Numbers is

-- CONSTRUCTORS :

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

end Triple_Double_Numbers;
