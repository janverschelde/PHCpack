with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;

with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package body Multprec_Natural64_Coefficients is

  function Powers_of_Ten ( k : natural32 ) return Array_of_Naturals is

  -- DESCRIPTION :
  --   Returns the first k powers of ten in an array of range 1..k.

    res : Array_of_Naturals(1..k);

  begin
    res(1) := 10;
    for i in 2..k loop
      res(i) := res(i-1)*10;
    end loop;
    return res;
  end Powers_of_Ten;

-- BASIC CONSTANTS :

  fact : constant natural64 := 10;
  expo : constant natural32 := 16;    -- instead of 8, doubled the size
 -- fact : constant natural := 16;  -- a bit faster, but less accurate
 -- expo : constant natural := 6;   -- also causes i/o problems for floats
  sqrt_base : constant natural64 := fact**(natural(expo/2));
  the_base : constant natural64 := sqrt_base*sqrt_base;
  sub_base : constant natural64 := the_base/fact;
  powten : constant Array_of_Naturals(1..expo) := Powers_of_Ten(expo);

-- SELECTORS :

  function Radix return natural64 is
  begin
    return fact;
  end Radix;

  function Base return natural64 is
  begin
    return the_base;
  end Base;

  function Exponent return natural32 is
  begin
    return expo;
  end Exponent;

  function Size_of_Coefficient ( n : natural64 ) return natural32 is
  begin
    for i in powten'range loop
      if n < powten(i)
       then return i;
      end if;
    end loop;
    return expo+1;
  end Size_of_Coefficient;

  function Digits_to_Normal ( n : Array_of_Naturals ) return natural32 is

    res : natural32 := 0;
    ind : natural32 := n'last+1;

  begin
    for i in reverse n'range loop
      if n(i) = 0
       then res := res + expo;
       else ind := i; exit;
      end if;
    end loop;
    if ind > n'last
     then return 0;
     else return expo - Size_of_Coefficient(n(ind));
    end if;
  end Digits_to_Normal;

-- COMPARISONS :

  function Equal ( n1,n2 : Array_of_Naturals ) return boolean is
  begin
    if n1'first /= n2'first then
      return false;
    elsif n1'last /= n2'last then
      return false;
    else
      for i in n1'range loop
        if n1(i) /= n2(i)
         then return false;
        end if;
      end loop;
      return true;
    end if;
  end Equal;

  function "<" ( n1,n2 : Array_of_Naturals ) return boolean is

    indlast : natural32;

  begin
    if n1'last < n2'last then                             -- compare sizes
      for i in reverse n1'last+1..n2'last loop
        if n2(i) > 0
         then return true;
        end if;
      end loop;
      indlast := n1'last;
    elsif n1'last > n2'last then
      for i in reverse n2'last+1..n1'last loop
        if n1(i) > 0
         then return false;
        end if;
      end loop;
      indlast := n2'last;
    else
      indlast := n1'last;
    end if;
    for i in reverse 0..indlast loop     -- look for first different leading
      if n1(i) < n2(i) then
        return true;
      elsif n1(i) > n2(i) then
        return false;
      end if;
    end loop;
    return false;                               -- equal coefficient vectors
  end "<";

  function ">" ( n1,n2 : Array_of_Naturals ) return boolean is

    indlast : natural32;

  begin
    if n1'last < n2'last then                               -- compare sizes
      for i in reverse n1'last+1..n2'last loop
        if n2(i) > 0
         then return false;
        end if;
      end loop;
      indlast := n1'last;
    elsif n1'last > n2'last then
      for i in reverse n2'last+1..n1'last loop
        if n1(i) > 0
         then return true;
        end if;
      end loop;
      indlast := n2'last;
    else
      indlast := n1'last;
    end if;
    for i in reverse 0..indlast loop     -- look for first different leading
      if n1(i) > n2(i) then
        return true;
      elsif n1(i) < n2(i) then
        return false;
      end if;
    end loop;
    return false;                               -- equal coefficient vectors
  end ">";

-- SHIFTS AND NORMALIZATIONS :

  procedure Digits_to_Left ( n : in out Array_of_Naturals; k : in natural32 ) is

    pt,pm,highest,lowest,carry : natural64;
    km,kd : natural32;

  begin
    if k < expo then
      km := k;
    else
      kd := k/expo;
      for i in reverse kd..n'last loop
        n(i) := n(i-kd);
      end loop;
      for i in 0..(kd-1) loop
        n(i) := 0;
      end loop;
      km := k mod expo;
    end if;
    if km /= 0 then
      pt := 10**(natural((expo-km)));
      pm := 10**natural(km);
      carry := 0;
      for i in n'range loop
        highest := n(i)/pt;
        lowest := n(i) mod pt;
        if lowest = 0
         then n(i) := carry;
         else n(i) := lowest*pm + carry;
        end if;
        carry := highest;
      end loop;
    end if;
  end Digits_to_Left;

  procedure Shift_Left ( n : in out Array_of_Naturals; ns : out natural32 ) is

    indlast,difind : natural32;

  begin
    indlast := n'last+1;
    for i in reverse n'range loop
      if n(i) /= 0
       then indlast := i; exit;
      end if;
    end loop;                             
    if indlast > n'last then
      ns := 0;                                       -- n = 0
    else
      if indlast /= n'last then                      -- n(indlast) /= 0
        difind := n'last - indlast;
        for i in reverse difind..n'last loop
          n(i) := n(i-difind);
        end loop;                                    -- n(n'last) /= 0
        for i in 0..difind-1 loop
          n(i) := 0;
        end loop;
        ns := difind*expo;
      else
        ns := 0;
      end if;
      while n(n'last) < sub_base loop
        Mul_Fact(n,fact);
        ns := ns+1;
      end loop;
    end if;
  end Shift_Left;

  procedure Shift_Right ( n : in out Array_of_Naturals; ns : out natural32 ) is

    indfirst,difind : natural32;
    cff : natural64;

  begin
    indfirst := n'last+1;
    for i in n'range loop
      if n(i) /= 0
       then indfirst := i; exit;
      end if;
    end loop;
    if indfirst > n'last then
      ns := 0;                                            -- n = 0
    else
      if indfirst /= 0 then                               -- n(indfirst) /= 0
        difind := n'last-indfirst;
        for i in 0..difind loop
          n(i) := n(i+indfirst);
        end loop;
        for i in difind+1..n'last loop
          n(i) := 0;
        end loop;                                 
        ns := indfirst*expo;
      else
        ns := 0;
      end if;                                            -- n(0) /= 0
      cff := n(0) mod fact;
      while (cff = 0) loop
        Small_Div(n,fact);
        ns := ns+1;
        cff := n(0) mod fact;
      end loop;
    end if;
  end Shift_Right;

-- OPERATIONS :

  procedure Add ( n1 : in out Array_of_Naturals; n2 : in natural64 ) is

    sum,carry : natural64;

  begin
    carry := n2;
    for i in n1'range loop
      sum := n1(i) + carry;
      n1(i) := sum mod the_base;
      carry := sum/the_base;
      exit when (carry = 0);
    end loop;
  end Add;

  procedure Add ( n1 : in out Array_of_Naturals;
                  n2 : in Array_of_Naturals ) is

    sum,carry : natural64;

  begin
    carry := 0;
    for i in n1'range loop
      if i <= n2'last then
        sum := n1(i) + n2(i) + carry;
        n1(i) := sum mod the_base;
        carry := sum/the_base;
      else
        sum := n1(i) + carry;
        n1(i) := sum mod the_base;
        carry := sum/the_base;
      end if;
      exit when ((i > n2'last) and (carry = 0));
    end loop;
  end Add;

  procedure Sub ( n1 : in out Array_of_Naturals;
                  n2 : in Array_of_Naturals ) is

    diff,carry : integer64;

  begin
    carry := 0;
    for i in n1'range loop
      if i <= n2'last then
        diff := integer64(n1(i)) - integer64(n2(i)) - carry;
        if diff >= 0 then
          n1(i) := natural64(diff);
          carry := 0;
        else
          diff := integer64(the_base) + diff;
          n1(i) := natural64(diff mod integer64(the_base));
          carry := 1;
        end if;
      elsif carry /= 0 then
        diff := integer64(n1(i)) - carry;
        if diff >= 0 then
          n1(i) := natural64(diff);
          carry := 0;
        else
          diff := integer64(the_base) + diff;
          n1(i) := natural64(diff mod integer64(the_base));
          carry := 1;
        end if;
      end if;
      exit when ((i > n2'last) and (carry = 0));
    end loop;
  end Sub;

  procedure Mul_Fact ( n : in out Array_of_Naturals; f : in natural64 ) is

    prod,carry : natural64;

  begin
    carry := 0;
    for i in n'range loop
      prod := n(i)*f + carry;
      n(i) := prod mod the_base;
      carry := prod/the_base;
    end loop;
  end Mul_Fact;

  procedure Acc_Add ( n : in out Array_of_Naturals; m1,m0 : in natural64;
                      k : in natural32; c : in natural64 ) is

    sum,carry : natural64;

  begin
    sum := n(k) + m0;
    n(k) := sum mod the_base;
    carry := sum/the_base;
    sum := n(k+1) + m1 + carry;
    n(k+1) := sum mod the_base;
    carry := sum/the_base;
    if carry /= 0
     then n(k+2) := n(k+2) + carry;
    end if;
    if c /= 0
     then n(k+2) := n(k+2) + c;
    end if;
  end Acc_Add;

  procedure Mul ( n1 : in out Array_of_Naturals; n2 : in natural64 ) is

    res : Array_of_Naturals(0..n1'last) := (0..n1'last => 0);
    acc0,acc1,n10,n11,n20,n21,prod,sum0,sum1,cr1 : natural64;

  begin
    n20 := n2 mod sqrt_base;
    n21 := n2/sqrt_base;
    for i in 0..n1'last-2 loop
      n10 := n1(i) mod sqrt_base;
      n11 := n1(i)/sqrt_base;
      acc1 := n21*n11;
      acc0 := n20*n10;
      Acc_Add(res,acc1,acc0,i,0);
      prod := n21*n10;
      if prod < sqrt_base then
        acc0 := prod*sqrt_base;
        acc1 := 0;
      else
        acc0 := (prod mod sqrt_base)*sqrt_base;
        acc1 := prod/sqrt_base;
      end if;
      cr1 := 0;
      prod := n20*n11;
      if prod < sqrt_base then
        sum0 := acc0 + prod*sqrt_base;
        if sum0 < the_base then
          acc0 := sum0;
        else
          acc0 := sum0 mod the_base;
          sum1 := acc1 + sum0/the_base;
          if sum1 < the_base
           then acc1 := sum1;
           else acc1 := sum1 mod the_base;
                cr1 := sum1/the_base;
          end if;
        end if;
      else
        sum0 := acc0 + (prod mod sqrt_base)*sqrt_base;
        if sum0 < the_base then
          acc0 := sum0;
        else
          acc0 := sum0 mod the_base;
          sum1 := acc1 + sum0/the_base;
          if sum1 < the_base
           then acc1 := sum1;
           else acc1 := sum1 mod the_base;
                cr1 := sum1/the_base;
          end if;
        end if;
        sum1 := acc1 + prod/sqrt_base;
        if sum1 < the_base
         then acc1 := sum1;
         else acc1 := sum1 mod the_base;
              cr1 := cr1 + sum1/the_base;
        end if;
      end if;
      Acc_Add(res,acc1,acc0,i,cr1);
    end loop;
    n1 := res;
  end Mul;

  function Mul ( n1,n2 : Array_of_Naturals ) return Array_of_Naturals is

    res : Array_of_Naturals(0..n1'last+n2'last+2);
    prod : Array_of_Naturals(0..n1'last+2);
    sum,carry : natural64;
    sizeres : natural32;

  begin
    prod(n1'range) := n1;
    prod(prod'last-1) := 0;
    prod(prod'last) := 0;
    Mul(prod,n2(0));
    res(prod'range) := prod;
    for i in prod'last+1..res'last loop
      res(i) := 0;
    end loop;
    for i in 1..n2'last loop
      prod(n1'range) := n1;
      prod(prod'last-1) := 0;
      prod(prod'last) := 0;
      Mul(prod,n2(i));
      carry := 0;
      for j in prod'range loop
        sum := res(i+j) + prod(j) + carry;
        res(i+j) := sum mod the_base;
        carry := sum/the_base;
      end loop;
    end loop;
    sizeres := 0;
    for i in reverse res'range loop
      if res(i) /= 0
       then sizeres := i;
      end if;
      exit when (sizeres /= 0);
    end loop;
    if sizeres < n1'last
     then sizeres := n1'last;
    end if;
    return res(0..sizeres);
  end Mul;

  procedure Small_Div ( n1 : in out Array_of_Naturals; n2 : in natural64 ) is

    carry : natural64 := 0;
    sum : natural64;

  begin
    for i in reverse 1..n1'last loop
      sum := n1(i) + carry;
      n1(i) := sum/n2;
      carry := sum mod n2;
      carry := carry*the_base;
    end loop;
    sum := n1(0) + carry;
    n1(0) := sum/n2;
  end Small_Div;

  procedure Small_Div ( n1 : in Array_of_Naturals; n2 : in natural64;
                        q : out Array_of_Naturals; r : out natural64 ) is

    carry : natural64 := 0;
    sum : natural64;

  begin
    for i in reverse 1..n1'last loop
      sum := n1(i) + carry;
      q(i) := sum/n2;
      carry := sum mod n2;
      carry := carry*the_base;
    end loop;
    sum := n1(0) + carry;
    q(0) := sum/n2;
    r := sum mod n2;
  exception 
    when others => put("n2 = "); put(integer64(n2));
                   put("  n1'last = "); put(n1'last,1);
                   put("  carry = "); put(integer64(carry)); new_line;
      for i in reverse n1'range loop
        put("  n1("); put(i,1); put(") = "); put(integer64(n1(i)));
        if n1(i) <= sqrt_base
         then put(" <= "); put(integer64(sqrt_base)); put_line("  okay");
         else put(" > "); put(integer64(sqrt_base)); put_line("  bug!");
        end if;
      end loop;
      new_line; raise;
  end Small_Div;

  procedure Small_Div ( n1 : in out Array_of_Naturals; n2 : in natural64;
                        r : out natural64 ) is

    carry : natural64 := 0;
    sum : natural64;

  begin
    for i in reverse 1..n1'last loop
      sum := n1(i) + carry;
      n1(i) := sum/n2;
      carry := sum mod n2;
      carry := carry*the_base;
    end loop;
    sum := n1(0) + carry;
    n1(0) := sum/n2;
    r := sum mod n2;
  end Small_Div;

  procedure Acc_Div ( n1 : in out Array_of_Naturals; leadn1 : in natural32;
                      n2,qb,rb,rbhead,rbtail : in natural64;
                      q1,q0 : out natural64 ) is

  -- DESCRIPTION :
  --   This accumulated divisor handles division of n1 by a standard natural
  --   n2 (apply for n2 > radix), where n1 has at least two leading nonzeros.
  --   If n2 <= radix, then Small_Div is more appropriate and efficient.
  --   Repeated application of this Acc_Div implements the Big_Div below.

  -- REQUIRED : n2 /= 0, n1(leadn1) /= 0 and leadn1 > 0.

  -- ON ENTRY :
  --   n1       array with n1(leadn1) /= 0 and n1(i) = 0, i > leadn1;
  --   leadn1   leadn1 > 0, index to leading nonzero in n1;
  --   n2       divisor, different from zero, and preferably n2 > radix;
  --   qb       quotient of the base divided by n2, i.e.: qb = Base/n2;
  --   rb       remainder of the base after division by n2: rb = Base mod n2;
  --   rbhead   rb / sqrt_base;
  --   rbtail   rb mod sqrt_base.

  -- ON RETURN :
  --   n1       n1(leadn1) and n1(leadn1-1) contain remainder,
  --            note that repeated application is needed to ensure the
  --            remainder will become less than n2;
  --   q1       leading coefficient of quotient;
  --   q0       trailing coefficient of quotient.

    q11,q10,r1,r0,sum,acc,prod0,prod1,r1head,r1tail : natural64;

    procedure Accumulate_Products1 is
    begin
      if prod1 < sqrt_base then
        sum := prod0 + prod1*sqrt_base;
        if sum < the_base
         then prod0 := sum;
              prod1 := 0;
         else prod0 := sum mod the_base;
              prod1 := sum/the_base;
        end if;
      else
        sum := prod0 + (prod1 mod sqrt_base)*sqrt_base;
        if sum < the_base
         then prod0 := sum;
              prod1 := prod1/sqrt_base;
         else prod0 := sum mod the_base;
              prod1 := sum/the_base + prod1/sqrt_base;
        end if;
      end if;
    end Accumulate_Products1;

    procedure Accumulate_Products2 is
    begin
      if acc < sqrt_base then
        sum := prod0 + acc*sqrt_base;
        if sum < the_base
         then prod0 := sum;
         else prod0 := sum mod the_base;
              prod1 := prod1 + sum/the_base;
        end if;
      else
        sum := prod0 + (acc mod sqrt_base)*sqrt_base;
        if sum < the_base
         then prod0 := sum;
              prod1 := prod1 + acc/sqrt_base;
         else prod0 := sum mod the_base;
              prod1 := prod1 + acc/sqrt_base + sum/the_base;
        end if;
      end if;
    end Accumulate_Products2;

  begin
    q11 := n1(leadn1)/n2;
    r1 := n1(leadn1) mod n2;
    q10 := n1(leadn1-1)/n2;
    r0 := n1(leadn1-1) mod n2;
    sum := q10 + r1*qb;           -- r1 < n2 and qb = Base/n2 => r1*qb < Base
    if sum < the_base
     then q0 := sum;
          q1 := q11;
     else q0 := sum mod the_base;
          q1 := q11 + sum/the_base;
    end if;
   -- sum := r0 + r1*rb;       -- r1*rb could be as large as (base-1)(base-2)
    if r1 < sqrt_base then
      if rb < sqrt_base
       then prod0 := r1*rb;
            prod1 := 0;
       else prod0 := r1*rbtail;
            prod1 := r1*rbhead;
            Accumulate_Products1;
      end if;
    elsif rb < sqrt_base then
      prod0 := rb*(r1 mod sqrt_base);
      prod1 := rb*(r1/sqrt_base);
      Accumulate_Products1;
    else
      r1head := r1/sqrt_base; r1tail := r1 mod sqrt_base;
      prod0 := r1tail*rbtail;
      prod1 := r1head*rbhead;
      acc := rbtail*r1head;
      Accumulate_Products2;
      acc := rbhead*r1tail;
      Accumulate_Products2;
    end if;
    sum := r0 + prod0;
    if sum < the_base
     then n1(leadn1-1) := sum;
          n1(leadn1) := prod1;
     else n1(leadn1-1) := sum mod the_base;
          n1(leadn1) := prod1 + sum/the_base;
    end if;
  end Acc_Div;

  procedure Big_Div ( n1 : in Array_of_Naturals; n2 : in natural64;
                      q : out Array_of_Naturals; r : out natural64 ) is

    leadn1 : natural32;
    qb,rb,rbhead,rbtail,q1,q0,sum,carry : natural64;
    rest : Array_of_Naturals(n1'range);
    quot : Array_of_Naturals(q'range);

  begin
    for i in reverse n1'range loop
      if n1(i) /= 0
       then leadn1 := i; exit;
      end if;
    end loop;
    if leadn1 = 0 then
      r := n1(0) mod n2;
      q(0) := n1(0)/n2;
      q(1..q'last) := (1..q'last => 0);
    else
      qb := the_base/n2;
      rb := the_base mod n2;
      rbhead := rb/sqrt_base;
      rbtail := rb mod sqrt_base;
      rest := n1;
      quot := (quot'range => 0);
      while leadn1 > 0 loop
        Acc_Div(rest,leadn1,n2,qb,rb,rbhead,rbtail,q1,q0);
        sum := quot(leadn1-1) + q0;
        quot(leadn1-1) := sum mod the_base;
        carry := sum/the_base;
        sum := quot(leadn1) + q1 + carry;
        quot(leadn1) := sum mod the_base;
        carry := sum/the_base;
        if carry /= 0 then
          for i in leadn1+1..quot'last loop
            sum := quot(i) + carry;
            quot(i) := sum mod the_base;
            carry := sum/the_base;
            exit when (carry = 0);
          end loop;
        end if;
        if rest(leadn1) = 0 then
          leadn1 := leadn1 - 1;
          if ((leadn1 > 0) and then (rest(leadn1) = 0))
           then leadn1 := leadn1 - 1;
          end if;
        end if;
      end loop;
      r := rest(0) mod n2;
      sum := quot(0) + rest(0)/n2;
      q(0) := sum mod the_base;
      carry := sum/the_base;
      for i in 1..quot'last loop
        sum := quot(i) + carry;
        q(i) := sum mod the_base;
        carry := sum/the_base;
      end loop;
    end if;
  end Big_Div;

  procedure Big_Div ( n1 : in out Array_of_Naturals; n2 : in natural64;
                      r : out natural64 ) is

    quot : Array_of_Naturals(n1'range);

  begin
    Big_Div(n1,n2,quot,r);
    for i in quot'range loop
      n1(i) := quot(i);
    end loop;
  end Big_Div;

  procedure Aligned_Acc_Div
               ( n1 : in out Array_of_Naturals; leadn1 : in natural32;
                 n2 : in Array_of_Naturals; q : out natural64 ) is

  -- DESCRIPTION :
  --   This accumulated division is auxiliary in the division routine of
  --   two long coefficient vectors.  This version assumes that n1 and n2
  --   are ligned up, in the sense that n1(leadn1) >= n2(n2'last).

  -- REQUIRED :
  --   n1(leadn1) /= 0, leadn1 >= n2'last, n1(leadn1) >= n2(n2'last),
  --   and n2(n2'last) > Base/Radix.

  -- ON ENTRY :
  --   n1        number to divide: n1(leadn1) >= n2(n2'last);
  --   leadn1    n1(leadn1) /= 0 and leadn1 >= n2'last;
  --   n2        divisor, normalized so that n2(n2'last) > Base/Radix.

  -- ON RETURN :
  --   n1        contains the remainder in its leading entries;
  --   q         quotient of n1 divided by n2.

    quot : natural64 := n1(leadn1)/n2(n2'last); -- first approx for quotient
    rest : Array_of_Naturals(n1'range);
    prod : Array_of_Naturals(n2'range);
    inc : natural32;
    acc,carry : natural64;
    diff : integer64;

  begin
   -- put_line("Entered Aligned Acc Div");
    loop                                        -- compute n2*quot
      carry := 0;
      for i in n2'range loop
        acc := n2(i)*quot + carry;
        prod(i) := acc mod the_base;
        carry := acc/the_base;
      end loop;
      exit when ((carry = 0) and then prod(prod'last) <= n1(leadn1));
      quot := quot - 1;
    end loop;
   -- put("     n1 : ");
   -- for i in n1'range loop
   --   put(n1(i));
   -- end loop;
   -- new_line;
   -- put("n2*quot : ");
   -- for i in prod'range loop
   --   put(prod(i));
   -- end loop;
   -- new_line;
    inc := leadn1-n2'last;
    rest := n1;
    for i in prod'range loop                    -- subtract n2*quot from n1
      diff := integer64(n1(i+inc)) - integer64(prod(i)) - integer64(carry);
      if diff >= 0
       then rest(i+inc) := natural64(diff);
            carry := 0;
       else rest(i+inc) := natural64(integer64(the_base) + diff);
            carry := 1;
      end if; 
    end loop;
   -- put("diff carry : "); put(carry,1);
    if carry > 0 then
      --put_line("   OVERFLOW in Aligned_Div");
      quot := quot - 1;                       -- decrease quot
      carry := 0;                             -- compute n2*quot - n1
      for i in n2'range loop
        diff := integer64(prod(i)) - integer64(n1(i+inc)) - integer64(carry);
        if diff >= 0
         then rest(i+inc) := natural64(diff); carry := 0;
         else rest(i+inc) := natural64(integer64(the_base) + diff);
              carry := 1;
        end if; 
      end loop;
      -- put("diff carry : "); put(carry,1); new_line;
      -- put("n2*quot - n1 : ");
      -- for i in rest'range loop
      --   put(rest(i));
      -- end loop;
      -- new_line;
      carry := 0;                             -- rest = n2 - (n2*quot-n1)
      for i in n2'range loop                  -- add n2 to the rest
        diff := integer64(n2(i)) - integer64(rest(i)) - integer64(carry);
        if diff >= 0
         then n1(i) := natural64(diff); carry := 0;
         else n1(i) := natural64(integer64(the_base) + diff); carry := 1;
        end if;
      end loop;
      for i in n2'last+1..n1'last loop
        n1(i) := rest(i);
      end loop;
    else --new_line;
      n1 := rest;
    end if;
   -- put("quot : "); put(quot,1); new_line;
   -- put("n1 : ");
   -- for i in n1'range loop
   --   put(n1(i));
   -- end loop;
   -- new_line;
    q := quot;
  end Aligned_Acc_Div;

  procedure Uneven_Acc_Div
               ( n1 : in out Array_of_Naturals; leadn1 : in natural32;
                 n2 : in Array_of_Naturals; q : out natural64;
                 shifts : out natural32 ) is

  -- DESCRIPTION :
  --   This accumulated division is auxiliary in the division routine of
  --   two long coefficient vectors.  This version assumes that n1 and n2
  --   are ligned up, in the sense that n1(leadn1) >= n2(n2'last).

  -- REQUIRED :
  --   n1(leadn1) /= 0, leadn1 > n2'last, n1(leadn1) < n2(n2'last),
  --   and n2(n2'last) > Base/Radix.

  -- ON ENTRY :
  --   n1        number to divide: n1(leadn1) < n2(n2'last);
  --   leadn1    n1(leadn1) /= 0 and leadn1 > n2'last;
  --   n2        divisor, normalized so that n2(n2'last) > Base/Radix.

  -- ON RETURN :
  --   n1        contains the remainder in its leading entries;
  --   q         quotient of n1 divided by n2;
  --   shifts    the number of shifts performed.

    inc,shiftctr : natural32;
    head,tail,quot,acc,carry,lostprod : natural64;
    prod : Array_of_Naturals(n2'range);
    rest : Array_of_Naturals(n1'range);
    diff : integer64;

  begin
   -- put_line("Entered Uneven Acc Div");
    head := n1(leadn1);                   -- head will gain digits from tail
    tail := n1(leadn1-1);                 --        to become >= n2(n2'last)
    shiftctr := 0;                        -- counts the number of shifts
    loop
      shiftctr := shiftctr + 1;
      inc := natural32(tail/the_base);    -- increment is first digit of tail
      head := head*fact + natural64(inc); -- shift and add increment to head
      exit when (head >= n2(n2'last));
      tail := tail - natural64(inc)*the_base;  -- remove increment from tail
      tail := tail*fact;                       -- shift tail up
    end loop;
    quot := head/n2(n2'last);
    prod(prod'last) := 0;
    loop                                  -- compute n2*quot
      carry := 0;
      for i in n2'range loop
        acc := n2(i)*quot + carry;
        if i < n2'last                    -- head > the_base !
         then prod(i) := acc mod the_base;
              carry := acc/the_base;
         else prod(i) := acc;
        end if;
      end loop;
      exit when (prod(n2'last) <= head);
      quot := quot - 1;
    end loop;
   -- put("quot : "); put(quot,1); new_line;
    lostprod := 0;                         -- adjust with shifts by dividing
    for i in 1..shiftctr loop 
      lostprod := (prod(0) mod fact)*sub_base + lostprod/fact;
      Small_Div(prod,fact);
      if quot < fact
       then quot := quot*sub_base;
       else quot := quot/fact;
      end if;
    end loop;
   -- put("quot after adjustments : "); put(quot,1); new_line;
   -- put("     n1 : ");
   -- for i in n1'range loop
   --   put(n1(i));
   -- end loop;
   -- new_line;
   -- put("     n2 : ");
   -- for i in n2'range loop
   --   put(n2(i));
   -- end loop;
   -- new_line;
   -- put("n2*quot : ");
   -- put(lostprod);
   -- for i in prod'range loop
   --   put(prod(i));
   -- end loop;
   -- new_line;
    inc := leadn1-n2'last;                       -- subtract n2*quot from n1
   -- put("inc : "); put(inc,1); new_line;
    rest := n1;
    diff := integer64(n1(inc-1)) - integer64(lostprod);
    if diff >= 0
     then rest(inc-1) := natural64(diff); carry := 0;
     else rest(inc-1) := natural64(diff + integer64(the_base)); carry := 1;
    end if;
    for i in n2'range loop
      diff := integer64(n1(i+inc)) - integer64(prod(i)) - integer64(carry);
      if diff >= 0
       then rest(i+inc) := natural64(diff); carry := 0;
       else rest(i+inc) := natural64(integer64(the_base) + diff); carry := 1;
      end if; 
    end loop;
   -- put("diff carry : "); put(carry,1);
    if carry > 0 then
      --put_line("   OVERFLOW in Uneven_Acc_Div");
      quot := quot - 1;           -- decrease quot and compute n2*quot - n1
      diff := integer64(lostprod) - integer64(n1(inc-1)); 
      if diff >= 0
       then rest(inc-1) := natural64(diff); carry := 0;
       else rest(inc-1) := natural64(diff + integer64(the_base)); carry := 1;
      end if;
      for i in n2'range loop
        diff := integer64(prod(i)) - integer64(n1(i+inc)) - integer64(carry);
        if diff >= 0
         then rest(i+inc) := natural64(diff); carry := 0;
         else rest(i+inc) := natural64(integer64(the_base) + diff); carry := 1;
        end if; 
      end loop;
      -- put("diff carry : "); put(carry,1); new_line;
      -- put("n2*quot - n1 : ");
      -- for i in rest'range loop
      --   put(rest(i));
      -- end loop;
      -- new_line;
      carry := 0;                             -- rest = n2 - (n2*quot-n1)
      for i in n2'range loop                  -- add n2 to the rest
        diff := integer64(n2(i)) - integer64(rest(i)) - integer64(carry);
        if diff >= 0
         then n1(i) := natural64(diff); carry := 0;
         else n1(i) := natural64(integer64(the_base) + diff); carry := 1;
        end if;
      end loop;
      for i in n2'last+1..n1'last loop
        n1(i) := rest(i);
      end loop;
    else --new_line;
      n1 := rest;
    end if;
   -- put("quot : "); put(quot,1); new_line;
   -- put("n1      : "); 
   -- for i in n1'range loop
   --   put(n1(i));
   -- end loop;
   -- new_line;
   -- put("leadn1 : "); put(leadn1,1); new_line;
    q := quot;
    shifts := shiftctr;
  end Uneven_Acc_Div;

  procedure Dominating ( n1,n2 : in Array_of_Naturals;
                         leadn1,leadn2 : in natural32;
                         dom,eq : out boolean ) is

  -- DESCRIPTION :
  --   Returns true if at least n1 >= n2, compared for indices in the
  --   common ranges of 0..leadn1 and 0..leadn2. 

    indn1,indn2 : natural32;

  begin
    dom := false;  -- is this the correct default ?
    eq := false;
    indn1 := leadn1;
    indn2 := leadn2;
    loop
      if n1(indn1) > n2(indn2) then
        dom := true; return;
      elsif n1(indn1) < n2(indn2) then
        dom := false; return;
      end if;
      exit when ((indn1 = 0) or (indn2 = 0));
      indn1 := indn1 - 1;
      indn2 := indn2 - 1;
    end loop;
    eq := true;   -- both are equal
  end Dominating;

  procedure Div ( n1,n2 : in Array_of_Naturals;
                  q,r : out Array_of_Naturals ) is

    dn1 : Array_of_Naturals(0..n1'last+1);
    dn2 : Array_of_Naturals(n2'range);
    dn1size,multnorm : natural32;
    quot : Array_of_Naturals(q'range) := (q'range => 0);

  begin
   -- put_line("Entering Div");
    dn1(n1'range) := n1;                             -- dn1 is shifted n1
    dn1(dn1'last) := 0;
    dn2 := n2;
    multnorm := 0;
    while dn2(dn2'last) < sub_base loop              -- dn2 is normalized n2
      Mul_Fact(dn1,fact);
      Mul_Fact(dn2,fact);
      multnorm := multnorm+1;
    end loop;
   -- multnorm := Digits_to_Normal(dn2);
   -- if multnorm /= 0 then
   --   Digits_to_Left(dn1,multnorm);
   --   Digits_to_Left(dn2,multnorm);
   -- end if;
    if dn1(dn1'last) = 0                             -- determine work size
     then dn1size := n1'last;                        -- of the rest of n1
     else dn1size := dn1'last;
    end if;
   -- put_line("In the Div on Array_of_Naturals : ");
   -- put("dn1 : ");
   -- for i in dn1'range loop
   --   put(dn1(i));
   -- end loop;
   -- new_line;
   -- put("dn2 : ");
   -- for i in dn2'range loop
   --   put(dn2(i));
   -- end loop;
   -- new_line;
   -- put("multnorm : "); put(multnorm,1); new_line;
    declare
      rest : Array_of_Naturals(0..dn1size);
      leadn1,shifts,quotind : natural32;
      quotfac,sum,carry : natural64;
      dom,eq : boolean;
    begin
      rest := dn1(rest'range);
      shifts := 0;
      for i in reverse rest'range loop               -- determine leadn1
        if rest(i) /= 0
         then leadn1 := i; exit;
        end if;
      end loop;
      while leadn1 >= dn2'last loop
       -- put("leadn1 : "); put(leadn1,1);
       -- put("  rest("); put(leadn1,1); put(") : ");
       -- put(rest(leadn1),1); new_line;
       -- put("dn2'last : "); put(dn2'last,1);
       -- put("  dn2("); put(dn2'last,1); put(") : ");
       -- put(dn2(dn2'last),1); new_line;
        quotind := 0;
        if rest(leadn1) = dn2(dn2'last) then
          --put("Checking dominating : ");
          Dominating(rest,dn2,leadn1-1,dn2'last-1,dom,eq);
          --if dom
          -- then put_line("okay.");
          -- else put_line("not okay.");
          --end if;
          if not dom then
            --put("rest : ");
            --for i in rest'range loop
            --  put(rest(i));
            --end loop;
            --new_line;
            --for i in dn2'range loop
            --  put(dn2(i));
            --end loop;
            --new_line;
            if eq then
              Sub(rest,dn2);
              quotind := leadn1 - dn2'last;
              carry := 1;
              for i in quotind..quot'last loop
                sum := quot(i) + carry;
                quot(i) := sum mod the_base;
                carry := sum/the_base;
                exit when (carry = 0);
              end loop;
              quotind := quot'last+1;  -- stop the iteration only if eq !!!
            end if;
            --put("rest(leadn1)  : "); put(rest(leadn1),1); new_line;
            --put("dn2(dn2'last) : "); put(dn2(dn2'last),1); new_line;
          end if;
        else
          dom := true;
        end if;
        exit when (quotind > quot'last);
        if rest(leadn1) >= dn2(dn2'last) and dom then
          Aligned_Acc_Div(rest,leadn1,dn2,quotfac);
          quotind := leadn1 - dn2'last;
        elsif leadn1 > dn2'last then
          Uneven_Acc_Div(rest,leadn1,dn2,quotfac,shifts);
          quotind := leadn1 - dn2'last - 1;
        else
          quotind := quot'last+1;            -- stop the iteration
        end if;
       -- put("quotind : "); put(quotind,1); new_line;
        exit when (quotind > quot'last);
        carry := quotfac;
        for i in quotind..quot'last loop             -- update the quotient
          sum := quot(i) + carry;
          quot(i) := sum mod the_base;
          carry := sum/the_base;
          exit when (carry = 0);
        end loop;
       -- put("accquot : ");
       -- for i in quot'range loop
       --   put(quot(i));
       -- end loop;
       -- new_line;
       -- put("n2*quot : ");
       -- declare
       --   chck : constant Array_of_Naturals := Mul(quot,n2);
       -- begin
       --   for i in chck'range loop
       --     put(chck(i));
       --   end loop;
       --   new_line;
       -- end;
        while (leadn1 >= dn2'last) and (rest(leadn1) = 0) loop
          leadn1 := leadn1-1;
        end loop;
      end loop;
      if multnorm > 0 then                               -- shift the rest 
        for i in 1..multnorm loop
          carry := 0;
          for i in reverse 1..leadn1 loop
            sum := rest(i) + carry;
            rest(i) := sum/fact;
            carry := sum mod fact;
            carry := carry*the_base;
          end loop;
          sum := rest(0) + carry;
          rest(0) := sum/fact;
        end loop;
      end if;
      r := rest(dn2'range);
    end;
    q := quot;
   -- put_line("Leaving Div");
  end Div;

  procedure Div ( n1 : in out Array_of_Naturals; n2 : in Array_of_Naturals;
                  r : out Array_of_Naturals ) is

    q : Array_of_Naturals(n1'range);

  begin
    Div(n1,n2,q,r);
    for i in q'range loop
      n1(i) := q(i);
    end loop;
  end Div;

end Multprec_Natural64_Coefficients;
