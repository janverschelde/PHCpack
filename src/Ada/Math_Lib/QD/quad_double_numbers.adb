with text_io;                            use text_io;
with Standard_Mathematical_Functions;
with Double_Double_Basics;               use Double_Double_Basics;
with Quad_Double_Renormalizations;       use Quad_Double_Renormalizations;

package body Quad_Double_Numbers is

-- CONSTRUCTORS :

  function Create ( i : integer ) return quad_double is
  begin
    return Create(integer32(i));
  end Create;

  function Create ( n : natural32 ) return quad_double is
  begin
    return Create(integer32(n));
  end Create;

  function Create ( n : natural64 ) return quad_double is
  begin
    return Create(integer64(n));
  end Create;

  function Create ( i : integer32 ) return quad_double is

    q : quad_double;

  begin
    q.hihi := double_float(i);
    q.lohi := 0.0; q.hilo := 0.0; q.lolo := 0.0;
    return q;
  end Create;

  function Create ( i : integer64 ) return quad_double is

    q : quad_double;

  begin
    q.hihi := double_float(i);
    q.lohi := 0.0; q.hilo := 0.0; q.lolo := 0.0;
    return q;
  end Create;

  function Create ( f : double_float ) return quad_double is

    q : quad_double;

  begin
    q.hihi := f;
    q.lohi := 0.0; q.hilo := 0.0; q.lolo := 0.0;
    return q;
  end Create;

  function Create ( d : double_double ) return quad_double is

    q : quad_double;

  begin
    q.hihi := hi_part(d);
    q.lohi := lo_part(d);
    q.hilo := 0.0; q.lolo := 0.0;
    return q;
  end Create;

  function Create ( hi,lo : double_double ) return quad_double is

    q : quad_double;

  begin
    q.hihi := hi_part(hi);
    q.lohi := lo_part(hi);
    q.hilo := hi_part(lo);
    q.lolo := lo_part(lo);
    return q;
  end Create;

  function Create ( hihi,lohi,hilo,lolo : double_float ) return quad_double is

    q : quad_double;

  begin
    q.hihi := hihi;
    q.lohi := lohi;
    q.hilo := hilo;
    q.lolo := lolo;
    return q;
  end Create;

-- SELECTORS :

  function hi_part ( q : quad_double ) return double_double is

    d : constant double_double := Create(q.hihi,q.lohi);

  begin
    return d;
  end hi_part;

  function hihi_part ( q : quad_double ) return double_float is
  begin
    return q.hihi;
  end hihi_part;

  function lohi_part ( q : quad_double ) return double_float is
  begin
    return q.lohi;
  end lohi_part;

  function lo_part ( q : quad_double ) return double_double is

    d : constant double_double := Create(q.hilo,q.lolo);

  begin
    return d;
  end lo_part;

  function hilo_part ( q : quad_double ) return double_float is
  begin
    return q.hilo;
  end hilo_part;

  function lolo_part ( q : quad_double ) return double_float is
  begin
    return q.lolo;
  end lolo_part;

-- COMPARISON and COPYING :

  function is_zero ( q : quad_double ) return boolean is
  begin
    return ((q.hihi = 0.0) and (q.lohi = 0.0)
        and (q.hilo = 0.0) and (q.lolo = 0.0));
  end is_zero;

  function is_one ( q : quad_double ) return boolean is
  begin
    return ((q.hihi = 1.0) and (q.lohi = 0.0)
        and (q.hilo = 0.0) and (q.lolo = 0.0));
  end is_one;

  function is_positive ( q : quad_double ) return boolean is
  begin
    return (q.hihi > 0.0);
  end is_positive;

  function is_negative ( q : quad_double ) return boolean is
  begin
    return (q.hihi < 0.0);
  end is_negative;

  function equal ( x,y : quad_double ) return boolean is
  begin
    return ((x.hihi = y.hihi) and (x.lohi = y.lohi)
        and (x.hilo = y.hilo) and (x.lolo = y.lolo));
  end equal;

  function equal ( x : quad_double; y : double_double ) return boolean is
  begin
    return ((x.hihi = hi_part(y)) and (x.lohi = lo_part(y))
        and (x.hilo = 0.0) and (x.lolo = 0.0));
  end equal;

  function equal ( x : quad_double; y : double_float ) return boolean is
  begin
    return ((x.hihi = y) and (x.lohi = 0.0)
        and (x.hilo = 0.0) and (x.lolo = 0.0));
  end equal;

  function "<" ( x,y : quad_double ) return boolean is
  begin
    return ((x.hihi < y.hihi)
         or (x.hihi = y.hihi and x.lohi < y.lohi)
         or (x.hihi = y.hihi and x.lohi = y.lohi and x.hilo < y.hilo)
         or (x.hihi = y.hihi and x.lohi = y.lohi and x.hilo = y.hilo and
             x.lolo < y.lolo));
  end "<";

  function "<" ( x : quad_double; y : double_double ) return boolean is
  begin
    return ((x.hihi < hi_part(y))
         or (x.hihi = hi_part(y) and x.lohi < lo_part(y))
         or (x.hihi = hi_part(y) and x.lohi = lo_part(y) and x.hilo < 0.0)
         or (x.hihi = hi_part(y) and x.lohi = lo_part(y) and x.hilo = 0.0 and
             x.lolo < 0.0));
  end "<";

  function "<" ( x : double_double; y : quad_double ) return boolean is
  begin
    return (y > x);
  end "<";

  function "<" ( x : quad_double; y : double_float ) return boolean is
  begin
    return ((x.hihi < y)
         or (x.hihi = y and x.lohi < 0.0)
         or (x.hihi = y and x.lohi = 0.0 and x.hilo < 0.0)
         or (x.hihi = y and x.lohi = 0.0 and x.hilo = 0.0 and x.lolo < 0.0));
  end "<";

  function "<" ( x : double_float; y : quad_double ) return boolean is
  begin
    return (y > x);
  end "<";

  function "<=" ( x,y : quad_double ) return boolean is
  begin
    return x < y or equal(x,y);
  end "<=";

  function "<=" ( x : quad_double; y : double_double ) return boolean is
  begin
    return ((x.hihi < hi_part(y))
         or (x.hihi = hi_part(y) and x.lohi < lo_part(y))
         or (x.hihi = lo_part(y) and x.lohi < 0.0)
         or (x.hihi = lo_part(y) and x.lohi = 0.0 and x.hilo < 0.0)
         or (x.hihi = lo_part(y) and x.lohi = 0.0 and x.hilo = 0.0 and
             x.lolo <= 0.0));
  end "<=";

  function "<=" ( x : double_double; y : quad_double ) return boolean is
  begin
    return (y >= x);
  end "<=";

  function "<=" ( x : quad_double; y : double_float ) return boolean is
  begin
    return ((x.hihi < y)
         or (x.hihi = y and x.lohi < 0.0)
         or (x.hihi = y and x.lohi = 0.0 and x.hilo < 0.0)
         or (x.hihi = y and x.lohi = 0.0 and x.hilo = 0.0 and x.lolo <= 0.0));
  end "<=";

  function "<=" ( x : double_float; y : quad_double ) return boolean is
  begin
    return (y >= x);
  end "<=";

  function ">" ( x,y : quad_double ) return boolean is
  begin
    return ((x.hihi > y.hihi)
         or (x.hihi = y.hihi and x.lohi > y.lohi)
         or (x.hihi = y.hihi and x.lohi = y.lohi and x.hilo > y.hilo)
         or (x.hihi = y.hihi and x.lohi = y.lohi and x.hilo = y.hilo and
             x.lolo > y.lolo));
  end ">";

  function ">" ( x : quad_double; y : double_double ) return boolean is
  begin
    return ((x.hihi > hi_part(y))
         or (x.hihi = hi_part(y) and x.lohi > lo_part(y))
         or (x.hihi = hi_part(y) and x.lohi = lo_part(y) and x.hilo > 0.0)
         or (x.hihi = hi_part(y) and x.lohi = lo_part(y) and x.hilo = 0.0 and
             x.lolo > 0.0));
  end ">";

  function ">" ( x : double_double; y : quad_double ) return boolean is
  begin
    return (y < x);
  end ">";

  function ">" ( x : quad_double; y : double_float ) return boolean is
  begin
    return ((x.hihi > y)
         or (x.hihi = y and x.lohi > 0.0)
         or (x.hihi = y and x.lohi = 0.0 and x.hilo > 0.0)
         or (x.hihi = y and x.lohi = 0.0 and x.hilo = 0.0 and x.lolo > 0.0));
  end ">";

  function ">" ( x : double_float; y : quad_double ) return boolean is
  begin
    return (y < x);
  end ">";

  function ">=" ( x,y : quad_double ) return boolean is
  begin
    return x > y or equal(x,y);
  end ">=";

  function ">=" ( x : quad_double; y : double_double ) return boolean is
  begin
    return ((x.hihi > hi_part(y))
         or (x.hihi = hi_part(y) and x.lohi > lo_part(y))
         or (x.hihi = hi_part(y) and x.lohi = lo_part(y) and x.lohi > 0.0)
         or (x.hihi = hi_part(y) and x.lohi = lo_part(y) and x.lohi = 0.0 and
             x.lolo > 0.0));
  end ">=";

  function ">=" ( x : double_double; y : quad_double ) return boolean is
  begin
    return (y <= x);
  end ">=";

  function ">=" ( x : quad_double; y : double_float ) return boolean is
  begin
    return ((x.hihi > y)
         or (x.hihi = y and x.lohi > 0.0)
         or (x.hihi = y and x.lohi = 0.0 and x.hilo > 0.0)
         or (x.hihi = y and x.lohi = 0.0 and x.hilo = 0.0 and x.lolo >= 0.0));
  end ">=";

  function ">=" ( x : double_float; y : quad_double ) return boolean is
  begin
    return (y <= x);
  end ">=";

  procedure copy ( x : in quad_double; y : in out quad_double ) is
  begin
    y.hihi := x.hihi; y.lohi := x.lohi;
    y.hilo := x.hilo; y.lolo := x.lolo;
  end copy;

-- Absolute value and type casts :

  function to_int ( x : quad_double ) return integer32 is
  begin
    return integer32(x.hihi);
  end to_int;

  function to_double ( x : quad_double ) return double_float is
  begin
    return x.hihi;
  end to_double;

  function to_double_double ( x : quad_double ) return double_double is

    d : constant double_double := hi_part(x);

  begin
    return d;
  end to_double_double;

  function to_triple_double ( x : quad_double ) return triple_double is

    t : constant triple_double := create(x.hihi,x.lohi,x.hilo);

  begin
    return t;
  end to_triple_double;

  function "abs" ( x : quad_double ) return quad_double is

    q : quad_double;

  begin
    if x.hihi < 0.0 then
      q.hihi := -x.hihi; q.lohi := -x.lohi;
      q.hilo := -x.hilo; q.lolo := -x.lolo;
    else
      q.hihi := x.hihi; q.lohi := x.lohi;
      q.hilo := x.hilo; q.lolo := x.lolo;
    end if;
    return q;
  end "abs";

  function AbsVal ( x : quad_double ) return quad_double is

    q : quad_double;

  begin
    if x.hihi < 0.0 then
      q.hihi := -x.hihi; q.lohi := -x.lohi;
      q.hilo := -x.hilo; q.lolo := -x.lolo;
    else
      q.hihi := x.hihi; q.lohi := x.lohi;
      q.hilo := x.hilo; q.lolo := x.lolo;
    end if;
    return q;
  end AbsVal;

  function floor ( x : quad_double ) return quad_double is

    res : quad_double;
    x0,x1,x2,x3 : double_float;

  begin
    x1 := 0.0; x2 := 0.0; x3 := 0.0;
    x0 := double_float'floor(x.hihi);
    if x0 = x.hihi then
      x1 := double_float'floor(x.lohi);
      if x1 = x.lohi then
        x2 := double_float'floor(x.hilo);
        if x2 = x.hilo
         then x3 := double_float'floor(x.lolo);
        end if;
      end if;
      renorm4(x0,x1,x2,x3);
    end if;
    res := Create(x0,x1,x2,x3);
    return res;
  end floor;

  function nint ( x : quad_double ) return quad_double is

    res : quad_double;
    x0,x1,x2,x3 : double_float;

  begin
    x1 := 0.0; x2 := 0.0; x3 := 0.0;
    x0 := Double_Double_Numbers.nint(x.hihi);
    if x0 = x.hihi then      -- first double is already an integer
      x1 := Double_Double_Numbers.nint(x.lohi);
      if x1 = x.lohi then    -- second double is already an integer
        x2 := Double_Double_Numbers.nint(x.hilo);
        if x2 = x.hilo then  -- third double is already an integer
          x3 := Double_Double_Numbers.nint(x.lolo);
        else
          if abs(x2 - x.hilo) = 0.5 and x.lolo < 0.0
           then x2 := x2 - 1.0;
          end if;
        end if;
      else
        if abs(x1 - x.lohi) = 0.5 and x.hilo < 0.0
         then x1 := x1 - 1.0;
        end if;
      end if;
    else                     -- first double is not an integer
      if abs(x0 - x.hihi) = 0.5 and x.lohi < 0.0
       then x0 := x0 - 1.0;
      end if;
    end if;
    renorm4(x0,x1,x2,x3);
    res := Create(x0,x1,x2,x3);
    return res;
  end nint;

-- ARITHMETICAL OPERATIONS :

  function get ( q : quad_double; i : natural ) return double_float is

  -- DESCRIPTION :
  --   Returns the component of q corresponding to the index i
  --   of the array representation of the quad double.

  begin
    case i is
      when 0 => return q.hihi;
      when 1 => return q.lohi;
      when 2 => return q.hilo;
      when others => return q.lolo;
    end case;
  end get;

  procedure Assign ( q : in out quad_double; i : natural;
                     v : in double_float ) is

  -- DESCRIPTION :
  --   Assigns the value v to the component of q defined by the 
  --   index i of the array representation of the quad double.

  begin
    case i is
      when 0 => q.hihi := v;
      when 1 => q.lohi := v;
      when 2 => q.hilo := v;
      when others => q.lolo := v;
    end case;
  end Assign;

  procedure Update ( q : in out quad_double; i : natural;
                     v : in double_float ) is

  -- DESCRIPTION :
  --   Updates the component of q defined by the index i of the array
  --   representation of the quad double with the value v.

  begin
    case i is
      when 0 => q.hihi := q.hihi + v;
      when 1 => q.lohi := q.lohi + v;
      when 2 => q.hilo := q.hilo +v;
      when others => q.lolo := q.lolo + v;
    end case;
  end Update;

  function "+" ( x,y : quad_double ) return quad_double is

    res : quad_double := Create(0.0);
    i,j,k : integer;
    s,t : double_float;
    u,v : double_float;  -- double-length accumulator

  begin
    i := 0; j := 0; k := 0;
    if abs(get(x,i)) > abs(get(y,j))
     then u := get(x,i); i := i+1;
     else u := get(y,j); j := j+1;
    end if;
    if abs(get(x,i)) > abs(get(y,j))
     then v := get(x,i); i := i+1;  
     else v := get(y,j); j := j+1;
    end if;
    quick_two_sum(u,v,u,v);
    while k < 4 loop
      if (i >= 4 and j >= 4) then
        Assign(res,k,u);
        if k < 3
         then k := k+1; Assign(res,k,v);
        end if;
        exit;
      end if;
      if i >= 4 then
        t := get(y,j); j := j+1;
      elsif j >= 4  then
        t := get(x,i); i := i+1;
      elsif abs(get(x,i)) > abs(get(y,j)) then
        t := get(x,i); i := i+1;
      else 
        t := get(y,j); j := j+1;
      end if;
      s := 0.0;
      quick_three_accum(u,v,s,t);
      if s /= 0.0
       then Assign(res,k,s); k := k+1;
      end if;
    end loop;
    for k in i..3 loop                    -- add the rest
      Update(res,3,get(x,k));
    end loop;
    for k in j..3 loop
      Update(res,3,get(y,k));
    end loop;
    renorm4(res.hihi,res.lohi,res.hilo,res.lolo);
    return res;
  end "+";

  function "+" ( x : quad_double; y : double_double ) return quad_double is

    res : quad_double;
    s0,s1,s2,s3,t0,t1 : double_float;

  begin
    two_sum(x.hihi,hi_part(y),s0,t0);
    two_sum(x.lohi,lo_part(y),s1,t1);
    two_sum(s1,t0,s1,t0);
    s2 := x.hilo;
    three_sum(s2,t0,t1);
    two_sum(t0,x.lolo,s3,t0);
    t0 := t0 + t1;
    renorm5(s0,s1,s2,s3,t0);
    res.hihi := s0; res.lohi := s1;
    res.hilo := s2; res.lolo := s3;
    return res;
  end "+";

  function "+" ( x : double_double; y : quad_double ) return quad_double is

    res : constant quad_double := y + x;

  begin
    return res;
  end "+";

  function "+" ( x : quad_double; y : double_float ) return quad_double is

    res : quad_double;
    e : double_float;

  begin
    two_sum(x.hihi,y,res.hihi,e);
    two_sum(x.lohi,e,res.lohi,e);
    two_sum(x.hilo,e,res.hilo,e);
    two_sum(x.lolo,e,res.lolo,e);
    renorm5(res.hihi,res.lohi,res.hilo,res.lolo,e);
    return res;
  end "+";

  function "+" ( x : double_float; y : quad_double ) return quad_double is

    res : constant quad_double := y + x;

  begin
    return res;
  end "+";

  function "+" ( x,y : double_double ) return quad_double is

    res : quad_double;
    qdx : constant quad_double := Create(x);

  begin
    res := qdx + y;
    return res;
  end "+";

  function "+" ( x : double_double; y : double_float ) return quad_double is

    res : quad_double;
    qdx : constant quad_double := Create(x);

  begin
    res := qdx + y;
    return res;
  end "+";

  function "+" ( x : double_float; y : double_double ) return quad_double is

    res : constant quad_double := y + x;

  begin
    return res;
  end "+";

  function "+" ( x,y : double_float ) return quad_double is

    res : quad_double;
    qdx : constant quad_double := Create(x);

  begin
    res := qdx + y;
    return res;
  end "+";

  function "+" ( x : quad_double ) return quad_double is

    res : quad_double;

  begin
    res.hihi := x.hihi; res.lohi := x.lohi;
    res.hilo := x.hilo; res.lolo := x.lolo;
    return res;
  end "+";

  procedure Add ( x : in out quad_double; y : in quad_double ) is
  begin
    x := x + y;
  end Add;

  procedure Add ( x : in out quad_double; y : in double_double ) is
  begin
    x := x + y;
  end Add;

  procedure Add ( x : in out quad_double; y : in double_float ) is
  begin
    x := x + y;
  end Add;

  function "-" ( x,y : quad_double ) return quad_double is

    res : quad_double;
    min_y : constant quad_double := -y;

  begin
    res := x + min_y;
    return res;
  end "-";

  function "-" ( x : quad_double; y : double_double ) return quad_double is

    res : quad_double;
    min_y : constant double_double := -y;

  begin
    res := x + min_y;
    return res;
  end "-";

  function "-" ( x : double_double; y : quad_double ) return quad_double is

    res : quad_double;
    min_y : constant quad_double := -y;

  begin
    res := x + min_y;
    return res;
  end "-";

  function "-" ( x : quad_double; y : double_float ) return quad_double is

    res : quad_double;
    min_y : constant double_float := -y;

  begin
    res := x + min_y;
    return res;
  end "-";

  function "-" ( x : double_float; y : quad_double ) return quad_double is

    res : quad_double;
    min_y : constant quad_double := -y;

  begin
    res := x + min_y;
    return res;
  end "-";

  function "-" ( x,y : double_double ) return quad_double is

    res : quad_double;
    min_y : constant double_double := -y;

  begin
    res := x + min_y;
    return res;
  end "-";

  function "-" ( x : double_double; y : double_float ) return quad_double is

    res : quad_double;
    min_y : constant double_float := -y;

  begin
    res := x + min_y; 
    return res;
  end "-";

  function "-" ( x : double_float; y : double_double ) return quad_double is

    res : quad_double;
    min_y : constant double_double := -y;

  begin
    res := x + min_y;
    return res;
  end "-";

  function "-" ( x,y : double_float ) return quad_double is

    res : quad_double;
    min_y : constant double_float := -y;

  begin
    res := x + min_y;
    return res;
  end "-";

  function "-" ( x : quad_double ) return quad_double is

    res : quad_double;

  begin
    res.hihi := -x.hihi; res.lohi := -x.lohi;
    res.hilo := -x.hilo; res.lolo := -x.lolo;
    return res;
  end "-";

  procedure Min ( x : in out quad_double ) is
  begin
    x.hihi := -x.hihi; x.lohi := -x.lohi;
    x.hilo := -x.hilo; x.lolo := -x.lolo;
  end Min;

  procedure Sub ( x : in out quad_double; y : in quad_double ) is
  begin
    x := x - y;
  end Sub;

  procedure Sub ( x : in out quad_double; y : in double_double ) is
  begin
    x := x - y;
  end Sub;

  procedure Sub ( x : in out quad_double; y : in double_float ) is
  begin
    x := x - y;
  end Sub;

  function "*" ( x,y : quad_double ) return quad_double is

  -- a0 * b0                    0
  --      a0 * b1               1
  --      a1 * b0               2
  --           a0 * b2          3
  --           a1 * b1          4
  --           a2 * b0          5
  --                a0 * b3     6
  --                a1 * b2     7
  --                a2 * b1     8
  --                a3 * b0     9
  -- note: accurate version satisfying IEEE error bound

    res : quad_double;
    p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,q5 : double_float;
    p6,p7,p8,p9,q6,q7,q8,q9,r0,r1,t0,t1,s0,s1,s2 : double_float;

  begin
    two_prod(x.hihi,y.hihi,p0,q0);
    two_prod(x.hihi,y.lohi,p1,q1);
    two_prod(x.lohi,y.hihi,p2,q2);
    two_prod(x.hihi,y.hilo,p3,q3);
    two_prod(x.lohi,y.lohi,p4,q4);
    two_prod(x.hilo,y.hihi,p5,q5);
    three_sum(p1,p2,q0);              -- start accumulation
    three_sum(p2,q1,q2);              -- six-three sum 
    three_sum(p3,p4,p5);              --  of p2, q1, q2, p3, p4, p5
    two_sum(p2,p3,s0,t0);             -- compute (s0, s1, s2)
    two_sum(q1,p4,s1,t1);             --  = (p2, q1, q2) + (p3, p4, p5)
    s2 := q2 + p5;
    two_sum(s1,t0,s1,t0);
    s2 := s2 + (t0 + t1);
    two_prod(x.hihi,y.lolo,p6,q6);    -- O(eps^3) order terms
    two_prod(x.lohi,y.hilo,p7,q7);
    two_prod(x.hilo,y.lohi,p8,q8);
    two_prod(x.lolo,y.hihi,p9,q9);
    two_sum(q0,q3,q0,q3);             -- nine-two-sum of q0, s1, q3,
    two_sum(q4,q5,q4,q5);             --  q4, q5, p6, p7, p8, p9
    two_sum(p6,p7,p6,p7);
    two_sum(p8,p9,p8,p9);
    two_sum(q0,q4,t0,t1);             -- compute (t0, t1) 
    t1 := t1 + (q3 + q5);             --  = (q0, q3) + (q4, q5)
    two_sum(p6,p8,r0,r1);             -- compute (r0, r1) 
    r1 := r1 + (p7 + p9);             --  = (p6, p7) + (p8, p9) 
    two_sum(t0,r0,q3,q4);             -- compute (q3, q4) 
    q4 := q4 + (t1 + r1);             --  = (t0, t1) + (r0, r1) 
    two_sum(q3,s1,t0,t1);             -- compute (t0, t1) 
    t1 := t1 + q4;                    --  = (q3, q4) + s1 
                                      -- O(eps^4) terms -- nine-one-sum 
    t1 := t1 + x.lohi * y.lolo + x.hilo * y.hilo + x.lolo * y.lohi
        + q6 + q7 + q8 + q9 + s2;
    renorm5(p0,p1,s0,t0,t1);
    res.hihi := p0; res.lohi := p1;
    res.hilo := s0; res.lolo := t0;
    return res;
  end "*";

  function "*" ( x : quad_double; y : double_double ) return quad_double is

    res : quad_double;
    p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,s0,s1,s2,t0,t1 : double_float;

  begin
    two_prod(x.hihi,hi_part(y),p0,q0);
    two_prod(x.hihi,lo_part(y),p1,q1);
    two_prod(x.lohi,hi_part(y),p2,q2);
    two_prod(x.lohi,lo_part(y),p3,q3);
    two_prod(x.hilo,hi_part(y),p4,q4);
    three_sum(p1,p2,q0);
    three_sum(p2,p3,p4);          -- five-three-sum
    two_sum(q1,q2,q1,q2);
    two_sum(p2,q1,s0,t0);
    two_sum(p3,q2,s1,t1);
    two_sum(s1,t0,s1,t0);
    s2 := t0 + t1 + p4;
    p2 := s0;
    p3 := x.hilo * hi_part(y) + x.lolo * lo_part(y) + q3 + q4;
    three_sum2(p3,q0,s1);
    p4 := q0 + s2;
    renorm5(p0,p1,p2,p3,p4);
    res.hihi := p0; res.lohi := p1;
    res.hilo := p2; res.lolo := p3;
    return res;
  end "*";

  function "*" ( x : double_double; y : quad_double ) return quad_double is

    res : constant quad_double := y * x;

  begin
    return res;
  end "*";

  function "*" ( x : quad_double; y : double_float ) return quad_double is

    res : quad_double;
    p0,p1,p2,p3,q0,q1,q2,s4 : double_float;

  begin
    two_prod(x.hihi,y,p0,q0);
    two_prod(x.lohi,y,p1,q1);
    two_prod(x.hilo,y,p2,q2);
    p3 := x.lolo * y;
    res.hihi := p0;
    two_sum(q0,p1,res.lohi,res.hilo);
    three_sum(res.hilo,q1,p2);
    three_sum2(q1,q2,p3);
    res.lolo := q1;
    s4 := q2 + p2;
    renorm5(res.hihi,res.lohi,res.hilo,res.lolo,s4);
    return res;
  end "*";

  function "*" ( x : double_float; y : quad_double ) return quad_double is

    res : constant quad_double := y * x;

  begin
    return res;
  end "*";

  function "*" ( x,y : double_double ) return quad_double is

    res : quad_double;
    qdx : constant quad_double := Create(x);

  begin
    res := qdx * y;
    return res;
  end "*";

  function "*" ( x : double_double; y : double_float ) return quad_double is

    res : quad_double;
    qdx : constant quad_double := Create(x);

  begin
    res := qdx * y;
    return res;
  end "*";

  function "*" ( x : double_float; y : double_double ) return quad_double is

    res : quad_double;
    qdx : constant quad_double := Create(x);

  begin
    res := qdx * y;
    return res;
  end "*";

  function "*" ( x,y : double_float ) return quad_double is

    res : quad_double;
    qdx : constant quad_double := Create(x);

  begin
    res := qdx * y;
    return res;
  end "*";

  procedure Mul ( x : in out quad_double; y : in quad_double ) is
  begin
    x := x*y;
  end Mul;

  procedure Mul ( x : in out quad_double; y : in double_double ) is
  begin
    x := x*y;
  end Mul;

  procedure Mul ( x : in out quad_double; y : in double_float ) is
  begin
    x := x*y;
  end Mul;

  function Mul_pwr2 ( x : quad_double; y : double_float )
                    return quad_double is

    res : quad_double;

  begin
    res.hihi := x.hihi*y; res.lohi := x.lohi*y;
    res.hilo := x.hilo*y; res.lolo := x.lolo*y;
    return res;
  end Mul_pwr2;

  procedure Mul_pwr2 ( x : in out quad_double; y : in double_float ) is
  begin
    x.hihi := x.hihi*y; x.lohi := x.lohi*y;
    x.hilo := x.hilo*y; x.lolo := x.lolo*y;
  end Mul_pwr2;

  function "/" ( x,y : quad_double ) return quad_double is

    res,acc : quad_double;
    q0,q1,q2,q3,q4 : double_float;

  begin
    q0 := x.hihi/y.hihi;
    acc := q0*y;
    res := x - acc;          -- r = a - (b * q0);
    q1 := res.hihi/y.hihi;
    acc := q1*y;
    Sub(res,acc);            -- r -= (b * q1);
    q2 := res.hihi/y.hihi;
    acc := q2*y;
    Sub(res,acc);            -- r -= (b * q2);
    q3 := res.hihi/y.hihi;
    acc := q3*y;
    Sub(res,acc);            -- r -= (b * q3);
    q4 := res.hihi/y.hihi;
    renorm5(q0,q1,q2,q3,q4);
    res.hihi := q0; res.lohi := q1; res.hilo := q2; res.lolo := q3;
    return res;
  end "/";

  function "/" ( x : quad_double; y : double_double ) return quad_double is

    res : quad_double;            -- note: this is the accurate division
    q0,q1,q2,q3,q4 : double_float;
    qd_y : constant quad_double := Create(y);
    acc,r : quad_double;

  begin
    q0 := x.hihi/hi_part(y);
    acc := q0*qd_y;             -- r = x - q0 * qd_y
    r := x - acc;
    q1 := r.hihi/hi_part(y);
    acc := q1*qd_y;
    Sub(r,acc);                 -- r -= (q1 * qd_b);
    q2 := r.hihi/hi_part(y);
    acc := q2*qd_y;
    Sub(r,acc);                 -- r -= (q2 * qd_b);
    q3 := r.hihi/hi_part(y);
    acc := q3*qd_y;
    Sub(r,acc);                 -- r -= (q3 * qd_b);
    q4 := r.hihi/hi_part(y);
    renorm5(q0,q1,q2,q3,q4);
    res.hihi := q0; res.lohi := q1; res.hilo := q2; res.lolo := q3;
    return res;
  end "/";

  function "/" ( x : double_double; y : quad_double ) return quad_double is

    res : quad_double;
    qdx : constant quad_double := Create(x);

  begin
    res := qdx/y;
    return res;
  end "/";

  function "/" ( x : quad_double; y : double_float ) return quad_double is

  -- Strategy: compute approximate quotient using high order doubles,
  -- and then correct it 3 times using the remainder (like long division).

    res : quad_double;
    t0,t1,q0,q1,q2,q3 : double_float;
    dd_t : double_double;

  begin
    q0 := x.hihi/y;                  -- approximate quotient
    two_prod(q0,y,t0,t1);            -- compute the remainder a - q0*b
    dd_t := Create(t0,t1);           -- make double_double for sub
    res := x - dd_t;                 -- c = a - dd_real(t0, t1);
    q1 := res.hihi/y;                -- compute the first correction
    two_prod(q1,y,t0,t1);
    dd_t := Create(t0,t1);
    res := res - dd_t;               -- c -= dd_real(t0, t1); 
    q2 := res.hihi/y;                -- second correction to the quotient 
    two_prod(q2,y,t0,t1);
    dd_t := Create(t0,t1);
    res := res - dd_t;               -- c -= dd_real(t0, t1);
    q3 := res.hihi/y;                -- final correction to the quotient
    renorm4(q0,q1,q2,q3);
    res.hihi := q0; res.lohi := q1; res.hilo := q2; res.lolo := q3;
    return res;
  end "/";

  function "/" ( x : double_float; y : quad_double ) return quad_double is

    res : quad_double;
    qdx : constant quad_double := Create(x);

  begin
    res := qdx/y;
    return res;
  end "/";

  function "/" ( x,y : double_double ) return quad_double is

    res : quad_double;
    qdx : constant quad_double := Create(x);

  begin
    res := qdx/y;
    return res;
  end "/";

  function "/" ( x : double_double; y : double_float ) return quad_double is

    res : quad_double;
    qdx : constant quad_double := Create(x);

  begin
    res := qdx/y;
    return res;
  end "/";

  function "/" ( x : double_float; y : double_double ) return quad_double is

    res : quad_double;
    qdx : constant quad_double := Create(x);

  begin
    res := qdx/y;
    return res;
  end "/";

  function "/" ( x,y : double_float ) return quad_double is

    res : quad_double;
    qdx : constant quad_double := Create(x);

  begin
    res := qdx/y;
    return res;
  end "/";

  procedure Div ( x : in out quad_double; y : in quad_double ) is
  begin
    x := x/y;
  end Div;

  procedure Div ( x : in out quad_double; y : in double_double ) is
  begin
    x := x/y;
  end Div;

  procedure Div ( x : in out quad_double; y : in double_float ) is
  begin
    x := x/y;
  end Div;

  function sqr ( x : double_float ) return quad_double is

    res : quad_double;
    qdx : constant quad_double := Create(x);

  begin
    res := sqr(qdx);
    return res;
  end sqr;

  function sqr ( x : double_double ) return quad_double is

    res : quad_double;
    qdx : constant quad_double := Create(x);

  begin
    res := sqr(qdx);
    return res;
  end sqr;

  function sqr ( x : quad_double ) return quad_double is

    res : quad_double;

  -- quad-double ^ 2  = (x0 + x1 + x2 + x3) ^ 2
  -- = x0 ^ 2 + 2 x0 * x1 + (2 x0 * x2 + x1 ^ 2) + (2 x0 * x3 + 2 x1 * x2)

    p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,s0,s1,t0,t1 : double_float;

  begin
    two_sqr(x.hihi,p0,q0);
    two_prod(2.0*x.hihi,x.lohi,p1,q1);
    two_prod(2.0*x.hihi,x.hilo,p2,q2);
    two_sqr(x.lohi,p3,q3);
    two_sum(q0,p1,p1,q0);
    two_sum(q0,q1,q0,q1);
    two_sum(p2,p3,p2,p3);
    two_sum(q0,p2,s0,t0);
    two_sum(q1,p3,s1,t1);
    two_sum(s1,t0,s1,t0);
    t0 := t0 + t1;
    quick_two_sum(s1,t0,s1,t0);
    quick_two_sum(s0,s1,p2,t1);
    quick_two_sum(t1,t0,p3,q0);
    p4 := 2.0*x.hihi*x.lolo;
    p5 := 2.0*x.lohi*x.hilo;
    two_sum(p4,p5,p4,p5);
    two_sum(q2,q3,q2,q3);
    two_sum(p4,q2,t0,t1);
    t1 := t1 + p5 + q3;
    two_sum(p3,t0,p3,p4);
    p4 := p4 + q0 + t1;
    renorm5(p0,p1,p2,p3,p4);
    res.hihi := p0; res.lohi := p1; res.hilo := p2; res.lolo := p3;
    return res;
  end sqr;

  function "**" ( x : quad_double; n : integer ) return quad_double is

    res,acc : quad_double;
    absn : natural;

  begin
    if n = 0 then
      res := Create(1.0);
    else
      if n > 0
       then absn := n;
       else absn := -n;
      end if;
      copy(x,res);
      acc := Create(1.0);
      if absn > 1 then         -- use binary exponentiation
        while absn > 0 loop     
          if absn mod 2 = 1    -- if odd multiply by res
           then Mul(acc,res);  -- eventually absn = 1, so this executes
          end if;
          absn := absn/2;
          if absn > 0
           then res := sqr(res);
          end if;
        end loop;
      else
        copy(res,acc);
      end if;
      if n < 0                 -- compute reciprocal
       then res := 1.0/acc;
       else copy(acc,res);
      end if;
    end if;
    return res;
  end "**";

  function "**" ( x : quad_double; n : integer32 ) return quad_double is
  begin
    return x**integer(n);
  end "**";

  function "**" ( x : quad_double; n : integer64 ) return quad_double is

    res,acc : quad_double;
    absn : natural64;

  begin
    if n = 0 then
      res := Create(1.0);
    else
      if n > 0
       then absn := natural64(n);
       else absn := natural64(-n);
      end if;
      copy(x,res);
      acc := Create(1.0);
      if absn > 1 then         -- use binary exponentiation
        while absn > 0 loop     
          if absn mod 2 = 1    -- if odd multiply by res
           then Mul(acc,res);  -- eventually absn = 1, so this executes
          end if;
          absn := absn/2;
          if absn > 0
           then res := sqr(res);
          end if;
        end loop;
      else
        copy(res,acc);
      end if;
      if n < 0                 -- compute reciprocal
       then res := 1.0/acc;
       else copy(acc,res);
      end if;
    end if;
    return res;
  end "**";

  function ldexp ( x : quad_double; n : integer ) return quad_double is

    res : quad_double;
    function C_ldexp ( x : double_float; n : integer ) return double_float;
    pragma Import(C,C_ldexp,External_Name => "ldexp");

  begin
    res.hihi := C_ldexp(x.hihi,n);
    res.lohi := C_ldexp(x.lohi,n);
    res.hilo := C_ldexp(x.hilo,n);
    res.lolo := C_ldexp(x.lolo,n);
    return res;
  end ldexp;

  function "**" ( x,y : quad_double ) return quad_double is
  begin
    return exp(y*log(x));
  end "**";

  function "**" ( x : quad_double; y : double_double ) return quad_double is

    qd_y : constant quad_double := create(y);

  begin
    return x**qd_y;
  end "**";

  function "**" ( x : quad_double; y : double_float ) return quad_double is

    qd_y : constant quad_double := create(y);

  begin
    return x**qd_y;
  end "**";
 
  function exp ( x : quad_double ) return quad_double is

  -- Strategy:  We first reduce the size of x by noting that
  --    exp(kr + m * log(2)) = 2^m * exp(r)^k
  -- where m and k are integers.  By choosing m appropriately
  -- we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
  -- evaluated using the familiar Taylor series.  Reducing the
  -- argument substantially speeds up the convergence.

    function C_ldexp ( x : double_float; n : integer ) return double_float;
    pragma Import(C,C_ldexp,External_Name => "ldexp");

    res : quad_double;
    k : constant double_float := C_ldexp(1.0,16);
    inv_k : constant double_float := 1.0/k;
    e_0 : constant double_float :=  2.718281828459045091e+00; -- exp(1)[0]
    e_1 : constant double_float :=  1.445646891729250158e-16; -- exp(1)[1]
    e_2 : constant double_float := -2.127717108038176765e-33; -- exp(1)[2]
    e_3 : constant double_float :=  1.515630159841218954e-49; -- exp(1)[3]
    exp1 : constant quad_double := Create(e_0,e_1,e_2,e_3);
    log2_0 : constant double_float :=  6.931471805599452862e-01; -- log(2)[0]
    log2_1 : constant double_float :=  2.319046813846299558e-17; -- log(2)[1]
    log2_2 : constant double_float :=  5.707708438416212066e-34; -- log(2)[2]
    log2_3 : constant double_float := -3.582432210601811423e-50; -- log(2)[3]
    log2 : constant quad_double := Create(log2_0,log2_1,log2_2,log2_3);
    qd_eps : constant double_float := 6.077163357286271e-64; -- 2.0**(-54*4-2)
    --  := 1.21543267145725e-63;      -- 2^-209 -- but it should be 2^(-210)
    tol : constant double_float := inv_k*qd_eps;
    m : constant double_float := double_float'floor(x.hihi/log2_0 + 0.5);
    i_fac : array(0..14) of quad_double;
      -- inverse factorials for Taylor expansion
    p,s,t : quad_double;
    cnt : integer;

  begin
    if x.hihi <= -709.0 then
      res := Create(0.0);
    elsif x.hihi >= 709.0 then
      res := Create(-1.0);
    elsif is_zero(x) then
      res := Create(1.0);
    elsif is_one(x) then
      res := exp1;
    else
      i_fac(0)  := Create( 1.66666666666666657e-01, 9.25185853854297066e-18,
                           5.13581318503262866e-34, 2.85094902409834186e-50);
      i_fac(1)  := Create( 4.16666666666666644e-02, 2.31296463463574266e-18,
                           1.28395329625815716e-34, 7.12737256024585466e-51);
      i_fac(2)  := Create( 8.33333333333333322e-03, 1.15648231731787138e-19,
                           1.60494162032269652e-36, 2.22730392507682967e-53);
      i_fac(3)  := Create( 1.38888888888888894e-03,-5.30054395437357706e-20,
                          -1.73868675534958776e-36,-1.63335621172300840e-52);
      i_fac(4)  := Create( 1.98412698412698413e-04, 1.72095582934207053e-22,
                           1.49269123913941271e-40, 1.29470326746002471e-58);
      i_fac(5)  := Create( 2.48015873015873016e-05, 2.15119478667758816e-23,
                           1.86586404892426588e-41, 1.61837908432503088e-59);
      i_fac(6)  := Create( 2.75573192239858925e-06,-1.85839327404647208e-22,
                           8.49175460488199287e-39,-5.72661640789429621e-55);
      i_fac(7)  := Create( 2.75573192239858883e-07, 2.37677146222502973e-23,
                          -3.26318890334088294e-40, 1.61435111860404415e-56);
      i_fac(8)  := Create( 2.50521083854417202e-08,-1.44881407093591197e-24,
                           2.04267351467144546e-41,-8.49632672007163175e-58);
      i_fac(9)  := Create( 2.08767569878681002e-09,-1.20734505911325997e-25,
                           1.70222792889287100e-42, 1.41609532150396700e-58);
      i_fac(10) := Create( 1.60590438368216133e-10, 1.25852945887520981e-26,
                          -5.31334602762985031e-43, 3.54021472597605528e-59);
      i_fac(11) := Create( 1.14707455977297245e-11, 2.06555127528307454e-28,
                           6.88907923246664603e-45, 5.72920002655109095e-61);
      i_fac(12) := Create( 7.64716373181981641e-13, 7.03872877733453001e-30,
                          -7.82753927716258345e-48, 1.92138649443790242e-64);
      i_fac(13) := Create( 4.77947733238738525e-14, 4.39920548583408126e-31,
                          -4.89221204822661465e-49, 1.20086655902368901e-65);
      i_fac(14) := Create( 2.81145725434552060e-15, 1.65088427308614326e-31,
                          -2.87777179307447918e-50, 4.27110689256293549e-67);
      res := mul_pwr2(x - m*log2,inv_k);
      p := sqr(res);
      s := res + mul_pwr2(p,0.5);
      cnt := 0;
      loop
        Mul(p,res);
        t := i_fac(cnt)*p;
        Add(s,t);
        exit when abs(t.hihi) <= tol;
        cnt := cnt + 1;
        exit when (cnt >= 9);
      end loop;
      for i in 1..16 loop -- 16 times s = mul_pwr2(s,2.0) + sqr(s);
        p := Mul_pwr2(s,2.0);
        t := sqr(s);
        s := p + t;
      end loop;
      res := s + 1.0;
      cnt := integer(m);
      res := ldexp(res,cnt);
    end if;
    return res;
  end exp;

  function log ( x : quad_double ) return quad_double is

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

    res,acc : quad_double;

  begin
    if is_one(x) then
      res := Create(0.0);
    elsif x.hihi <= 0.0 then
      put_line("qd_log: argument is not positive");
      res := Create(-1.0);
    else
      res := Create(Standard_Mathematical_Functions.LN(x.hihi));
      for i in 1..3 loop       -- r = r + x*exp(-r) - 1
        acc := x*exp(-res);    -- acc = x*exp(-res)      
        Add(res,acc);          -- res = res + x*exp(-res)
        Sub(res,1.0);
      end loop;
    end if;
    return res;
  end log;

  function log10 ( x : quad_double ) return quad_double is

    res : quad_double;
    log10_0 : constant double_float :=  2.302585092994045901e+00;
    log10_1 : constant double_float := -2.170756223382249351e-16;
    log10_2 : constant double_float := -9.984262454465776570e-33;
    log10_3 : constant double_float := -4.023357454450206379e-49;
    logten : constant quad_double := Create(log10_0,log10_1,log10_2,log10_3);

  begin
    res := log(x)/logten;
    return res;
  end log10;

-- DESTRUCTOR :

  procedure clear ( q : in out quad_double ) is
  begin
    q.hihi := 0.0; q.lohi := 0.0;
    q.hilo := 0.0; q.lolo := 0.0;
  end clear;

end Quad_Double_Numbers;
