--with text_io;  use text_io;
--with Multprec_Natural_Coefficients_io; use Multprec_Natural_Coefficients_io;
--with Standard_Natural_Numbers_io; use Standard_Natural_Numbers_io;
--with Standard_Integer_Numbers_io; use Standard_Integer_Numbers_io;

with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package body Multprec_Natural_Coefficients is

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

  fact : constant natural32 := 10;
  expo : constant natural32 := 8;
 -- fact : constant natural := 16;  -- a bit faster, but less accurate
 -- expo : constant natural := 6;   -- also causes i/o problems for floats
  sqrt_base : constant natural32 := fact**natural(expo/2);
  the_base : constant natural32 := sqrt_base*sqrt_base;
  base64 : constant natural64 := natural64(the_base);
  sub_base : constant natural32 := the_base/fact;
  powten : constant Array_of_Naturals(1..expo) := Powers_of_Ten(expo);

-- SELECTORS :

  function Radix return natural32 is
  begin
    return fact;
  end Radix;

  function Base return natural32 is
  begin
    return the_base;
  end Base;

  function Exponent return natural32 is
  begin
    return expo;
  end Exponent;

  function Size_of_Coefficient ( n : natural32 ) return natural32 is
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

    km,kd,pt,pm,highest,lowest,carry : natural32;

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
      pt := 10**natural(expo-km);
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

    indfirst,difind,cff : natural32;

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

  procedure Add ( n1 : in out Array_of_Naturals; n2 : in natural32 ) is

    sum,carry : natural32;

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

    sum,carry : natural32;

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

    diff,carry : integer32;

  begin
    carry := 0;
    for i in n1'range loop
      if i <= n2'last then
        diff := integer32(n1(i)) - integer32(n2(i)) - carry;
        if diff >= 0 then
          n1(i) := natural32(diff);
          carry := 0;
        else
          diff := integer32(the_base) + diff;
          n1(i) := natural32(diff) mod the_base;
          carry := 1;
        end if;
      elsif carry /= 0 then
        diff := integer32(n1(i)) - carry;
        if diff >= 0 then
          n1(i) := natural32(diff);
          carry := 0;
        else
          diff := integer32(the_base) + diff;
          n1(i) := natural32(diff) mod the_base;
          carry := 1;
        end if;
      end if;
      exit when ((i > n2'last) and (carry = 0));
    end loop;
  end Sub;

  procedure Mul_Fact ( n : in out Array_of_Naturals; f : in natural32 ) is

    prod,carry : natural32;

  begin
    carry := 0;
    for i in n'range loop
      prod := n(i)*f + carry;
      n(i) := prod mod the_base;
      carry := prod/the_base;
    end loop;
  end Mul_Fact;

  procedure Short_Mul ( n1 : in out Array_of_Naturals; n2 : in natural32 ) is

    n2m : constant natural64 := natural64(n2);
    p : natural64 := natural64(n1(0))*n2m;
    carry : natural64 := p/base64;

  begin
    n1(0) := natural32(p mod base64);
    for i in 1..n1'last loop
      p := natural64(n1(i))*n2m + carry;
      n1(i) := natural32(p mod base64);
      carry := p/base64;
    end loop;
  end Short_Mul;

  function Mul ( n1,n2 : Array_of_Naturals ) return Array_of_Naturals is

    res : Array_of_Naturals(0..n1'last+n2'last+2);
    prod : Array_of_Naturals(0..n1'last+2);
    sum,carry,sizeres : natural32;

  begin
    prod(n1'range) := n1;
    prod(prod'last-1) := 0;
    prod(prod'last) := 0;
    Short_Mul(prod,n2(0)); 
    res(prod'range) := prod;
    for i in prod'last+1..res'last loop
      res(i) := 0;
    end loop;
    for i in 1..n2'last loop
      prod(n1'range) := n1;
      prod(prod'last-1) := 0;
      prod(prod'last) := 0;
      Short_Mul(prod,n2(i)); 
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

  procedure Small_Div ( n1 : in out Array_of_Naturals; n2 : in natural32 ) is

    carry : natural32 := 0;
    sum : natural32;

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

  procedure Small_Div ( n1 : in Array_of_Naturals; n2 : in natural32;
                        q : out Array_of_Naturals; r : out natural32 ) is

    carry : natural32 := 0;
    sum : natural32;

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
  end Small_Div;

  procedure Small_Div ( n1 : in out Array_of_Naturals; n2 : in natural32;
                        r : out natural32 ) is

    carry : natural32 := 0;
    sum : natural32;

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

  procedure Short_Div ( n1 : in Array_of_Naturals; n2 : in natural32;
                        q : out Array_of_Naturals; r : out natural32 ) is

    n2d : constant natural64 := natural64(n2);
    acc : natural64 := natural64(n1(n1'last));
    quot : natural64 := acc/n2d;
    rest : natural64 := acc mod n2d;

  begin
    q(q'last) := natural32(quot);
    for i in reverse n1'first..n1'last-1 loop
      acc := rest*base64 + natural64(n1(i));
      quot := acc/n2d;
      rest := acc mod n2d;
      if quot < base64 then
        q(i) := natural32(quot);
      else
        q(i) := natural32(quot mod base64);
        q(i+1) := q(i+1) + natural32(quot/base64);
      end if;
    end loop;
    r := natural32(rest);
  end Short_Div;

  procedure Short_Div ( n1 : in out Array_of_Naturals; n2 : in natural32;
                        r : out natural32 ) is

    n2d : constant natural64 := natural64(n2);
    acc : natural64 := natural64(n1(n1'last));
    quot : natural64 := acc/n2d;
    rest : natural64 := acc mod n2d;

  begin
    n1(n1'last) := natural32(quot);
    for i in reverse n1'first..n1'last-1 loop
      acc := rest*base64 + natural64(n1(i));
      quot := acc/n2d;
      rest := acc mod n2d;
      if quot < base64 then
        n1(i) := natural32(quot);
      else
        n1(i) := natural32(quot mod base64);
        n1(i+1) := n1(i+1) + natural32(quot/base64);
      end if;
    end loop;
    r := natural32(rest);
  end Short_Div;

  function To_Natural64 ( n : Array_of_Naturals ) return natural64 is

  -- DESCRIPTION :
  --   Returns the representation of n as a 64-bit natural number,
  --   using only n(0) and n(1), so the answer will be wrong for
  --   larger values of n.

    res : constant natural64
        := natural64(n(1))*base64 + natural64(n(0));

  begin
    return res;
  end To_Natural64;

  function To_Natural64 ( n : Array_of_Naturals; lead : natural32 )
                        return natural64 is

  -- DESCRIPTION :
  --   Returns the representation of n as a 64-bit natural number,
  --   using only n(lead-1) and n(lead).

  -- REQUIRED : lead and lead-1 are in n'range.

    res : constant natural64
        := natural64(n(lead))*base64 + natural64(n(lead-1));

  begin
    return res;
  end To_Natural64;

  --function Double ( n : Array_of_Naturals; lead : natural32 )
  --                return double_float is

  -- DESCRIPTION :
  --   Returns the double float approximation of the leading part
  --   of the number in n.

  -- REQUIRED : 
  --   The highest nonzero position of n is in lead and lead > 1.

  --  res : double_float;

  --begin
  --  res := double_float(n(lead));
  --  res := double_float(n(lead-1)) + double_float(the_base)*res;
  --  res := double_float(n(lead-2)) + double_float(the_base)*res;
  --  return res;
  --end Double;

  function Alignment2 ( w : Array_of_Naturals; lead : natural32;
                        p : Array_of_Naturals ) return integer32 is

  -- DESCRIPTION :
  --   Determines the alignment of p with respect to w
  --   before subtracting p from w, auxiliary for Sub_Div2.

  -- REQUIRED : w(lead) /= 0 and p'range = 0..2.

  begin
    if p(2) /= 0 then
      return integer32(lead) - 2;
    elsif w(lead) < p(1) then
      return integer32(lead) - 2;
    else
      return integer32(lead) - 1;
    end if;
  end Alignment2;

  procedure Sub_Div2 ( w : in out Array_of_Naturals;
                       p : in Array_of_Naturals; least : in natural32;
                       fail : out boolean ) is

  -- DESCRIPTION :
  --   Subtracts from w the product of quotient times divisor in p.
  --   This is a subroutine of Short_Div2.
  --   If w - p turns out negative, then fail is true 
  --   and the estimate for the quotient was too high.

  -- REQUIRED :
  --   w(lead) is highest nonzero number in w and p'range = 0..2.

  -- ON ENTRY :
  --   w         accumulated remainder;
  --   p         product of quotient with divisor;
  --   least     result of Alignment2(w,lead,p).

  -- ON RETURN :
  --   w         result of w - p, provided fail is false;
  --   fail      true if w - p < 0, false otherwise.

    backup : constant Array_of_Naturals(w'range) := w;
    diff : integer32;
    carry : natural32 := 0;

  begin
   -- put("before sub w = "); put(w); new_line;
   -- put("           p = "); put(p); new_line;
    fail := false;
    for i in 0..natural32(1) loop
      diff := integer32(w(least+i)) - integer32(p(i)) - integer32(carry);
      if diff >= 0 then
        w(least+i) := natural32(diff);
        carry := 0;
      else
        w(least+i) := natural32(integer32(the_base) + diff);
        carry := 1;
      end if;
    end loop;
    if least+2 <= w'last then
      diff := integer32(w(least+2)) - integer32(p(2)) - integer32(carry);
      if diff >= 0 then
        w(least+2) := natural32(diff);
      else
        fail := true;
      end if;
    elsif p(2) /= 0 or carry /= 0 then
      fail := true;
    end if;
    if fail
     then w := backup; -- put_line("restoring to backup ...");
    end if;
   -- put(" after sub w = "); put(w); new_line;
  end Sub_Div2;

  procedure Add_Quot2 ( q : in out Array_of_Naturals; qlast : in out natural32;
                        quot : in natural32; undershot : in integer32 ) is

  -- DESCRIPTION :
  --   Updates the accumulated quotient in q, for the range 0..qlast,
  --   with the new quotient quot.  If undershot is positive, then this
  --   addition is just a shift, otherwise (if zero) a true addition.
  --   This is an auxiliary routine to Short_Div2.

  begin
   -- put("in Add_Quot2, undershot = "); put(undershot,1); new_line;
   -- put("before q = "); put(q); new_line;
    if undershot > 0 then -- just shift
      for i in reverse 0..qlast loop
        if i+natural32(undershot) <= q'last        -- patch
         then q(i+natural32(undershot)) := q(i);
        end if;
      end loop;
      q(0) := quot;
    else
      if undershot = 0
       then Add(q,quot);
       else Add(q(natural32(-undershot)..q'last),quot);
      end if;
      if qlast < q'last and then q(qlast+1) /= 0
       then qlast := qlast + 1;
      end if;
    end if;
   -- put(" after q = "); put(q); new_line;
  end Add_Quot2;

  function Est_Quot2 ( n1 : Array_of_Naturals; lead : natural32;
                       -- n2 : Array_of_Naturals;
                       n2d : natural64 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the estimate for the quotient as subroutine of Short_Div2.

  -- ON ENTRY :
  --   n1       numerator, number to be divided with n1(lead) /= 0;
  --   lead     largest index in n1 so n1(lead) /= 0;
  --   n2       denominator, only n2(1) and n2(0) play a role;
  --   n2d      representation of n2 as a 64-bit natural number.   

    res : natural32;
    n1d : constant natural64 := To_Natural64(n1,lead);
    f1,f2,r : double_float;

  begin
   -- put("in Est_Quot2, lead = "); put(lead,1); new_line;
   -- put("n1 = "); put(n1(0..lead)); new_line;
   -- put("n2 = "); put(n2); new_line;
   -- put("n1d = "); Standard_Natural_Numbers_io.put(n1d); new_line;
   -- put("n2d = "); Standard_Natural_Numbers_io.put(n2d); new_line;
    if n1d >= n2d then
      res := natural32(n1d/n2d);
     -- put("computing estimate the natural way : "); put(res,1); new_line;
    elsif lead <= 1 then
      res := 0;
    else
      f1 := double_float(n1(lead));
      f1 := double_float(n1(lead-1)) + double_float(the_base)*f1;
      f1 := double_float(n1(lead-2)) + double_float(the_base)*f1;
      f2 := double_float(n2d);
      r := f1/f2;
      res := natural32(r);
     -- put("using floats to estimate quotient : "); put(res,1); new_line;
      if r < double_float(res) then
        if res > 0
         then res := res - 1;
        end if;
      end if;
    end if;
   -- if res > the_base
   --  then put_line("Estimate larger than base, BUG!!!");
         -- put("n1 = "); put(n1(0..lead)); new_line;
         -- put("n2 = "); put(n2); new_line;
         -- put("n1d = "); Standard_Natural_Numbers_io.put(n1d); new_line;
         -- put("n2d = "); Standard_Natural_Numbers_io.put(n2d); new_line;
   -- end if;
    return res;
  end Est_Quot2;

  procedure Make_Prod2 ( p : in out Array_of_Naturals;
                         q : in natural64; n2 : in Array_of_Naturals ) is

  -- DESCRIPTION :
  --   Computes in p the product of q with n2.

  -- REQUIRED : p'range = 0..2, q < Base.

    n2d0 : constant natural64 := natural64(n2(0));
    n2d1 : constant natural64 := natural64(n2(1));
    prod,carry : natural64;

  begin
    prod := n2d0*q;
    p(0) := natural32(prod mod base64);
    carry := prod/base64;
    prod := n2d1*q + carry;
    p(1) := natural32(prod mod base64);
    carry := prod/base64;
    p(2) := natural32(carry);
  end Make_Prod2;

  procedure Short_Div2 ( n1,n2 : in Array_of_Naturals;
                         q,r : out Array_of_Naturals ) is

    leadn1 : natural32 := n1'last;
    qlast : natural32 := 1;
    w : Array_of_Naturals(n1'range) := n1;
    n2d : constant natural64 := To_Natural64(n2);
    quot : natural32;
    undershot : integer32;
    prevleast : natural32 := leadn1+1;
    least : integer32;
    prod : Array_of_Naturals(0..2);
    overshot : boolean;

  begin
   -- if n2(1) = 0
   --  then put_line("in Short_Div2 with n2(1) = 0, BUG!!!");
   -- end if;
    q := (q'range => 0);
    while w(leadn1) = 0 and leadn1 > 1 loop
      leadn1 := leadn1 - 1;
    end loop;
    quot := Est_Quot2(n1,leadn1,n2d); -- Est_Quot2(n1,leadn1,n2,n2d);
   -- put("estimate for quotient : "); put(quot,1);
   -- put("  leadn1 : "); put(leadn1,1); new_line;
    undershot := 0;
    for i in n1'range loop
     -- put("********* STEP "); put(i,1); put_line(" *********");
     -- put("q = "); put(q); new_line;
      loop
        Make_Prod2(prod,natural64(quot),n2);
        least := Alignment2(w,leadn1,prod);
        Sub_Div2(W,prod,natural32(least),overshot);
       -- put("prod = "); put(prod); new_line;
       -- put("W = "); put(W); new_line;
        exit when not overshot;
       -- put_line("overshot quotient!");
        quot := quot - 1;
       -- if quot = 0 and least > 0 then
       --   prevleast := least;
       --   least := least - 1;
       --   quot := the_base - 1;
       --   Sub
       -- end if;
        exit when quot = 0;
      end loop;
      exit when quot = 0;
      if prevleast <= leadn1
       then undershot := integer32(prevleast) - least;
      end if;
      if undershot > 0
       then qlast := qlast + natural32(undershot);
      end if;
      Add_Quot2(q,qlast,quot,undershot);
      while w(leadn1) = 0 and leadn1 > 0 loop
        leadn1 := leadn1 - 1;
      end loop;
      exit when (w(leadn1) = 0);
      exit when (w(0..leadn1) < n2);
      prevleast := natural32(least);
      quot := Est_Quot2(w,leadn1,n2d); -- Est_Quot2(w,leadn1,n2,n2d);
     -- put("estimate for quotient : "); put(quot,1);
     -- put("  leadn1 : "); put(leadn1,1); new_line;
      exit when (quot = 0);
    end loop;
   -- put_line("at the exit w = "); put(W(0..leadn1)); new_line;
   -- put("q = "); put(q); new_line;
   -- put("  qlast = "); put(qlast,1);
   -- put("  least = "); put(least,1); new_line;
    if least > 0 then
      for i in 0..(qlast-1) loop
        exit when (i+natural32(least) > q'last);
        q(i+natural32(least)) := q(i);
      end loop;
      for i in 0..least-1 loop
        q(natural32(i)) := 0;
      end loop;
    end if;
    r := W(r'range);
  end Short_Div2;

  procedure Short_Div2 ( n1 : in out Array_of_Naturals;
                         n2 : in Array_of_Naturals;
                         r : out Array_of_Naturals ) is

    q : Array_of_Naturals(n1'range);

  begin
    Short_Div2(n1,n2,q,r);
    n1 := q;
  end Short_Div2;

  procedure Aligned_Acc_Div
               ( n1 : in out Array_of_Naturals; leadn1 : in natural32;
                 n2 : in Array_of_Naturals; q : out natural32 ) is

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

    quot : natural32 := n1(leadn1)/n2(n2'last);  -- first approx for quotient
    rest : Array_of_Naturals(n1'range);
    prod : Array_of_Naturals(n2'range);
    acc,carry,inc : natural32;
    diff : integer32;

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
   -- put("   quot : "); put(quot,1); new_line;
   -- put("     n1 : "); put(n1); new_line;
   -- put("     n2 : "); put(n2); new_line;
   -- put("n2*quot : "); put(prod); new_line;
    inc := leadn1-n2'last;
    rest := n1;
    for i in prod'range loop                    -- subtract n2*quot from n1
      diff := integer32(n1(i+inc)) - integer32(prod(i)) - integer32(carry);
      if diff >= 0
       then rest(i+inc) := natural32(diff);                       carry := 0;
       else rest(i+inc) := natural32(integer32(the_base) + diff); carry := 1;
      end if; 
    end loop;
   -- put("diff carry : "); put(carry,1);
    if carry > 0 then
     -- put_line("   OVERFLOW in Aligned_Div");
      quot := quot - 1;                       -- decrease quot
      carry := 0;                             -- compute n2*quot - n1
      for i in n2'range loop
        diff := integer32(prod(i)) - integer32(n1(i+inc)) - integer32(carry);
        if diff >= 0
         then rest(i+inc) := natural32(diff);                       carry := 0;
         else rest(i+inc) := natural32(integer32(the_base) + diff); carry := 1;
        end if; 
      end loop;
     -- put("diff carry : "); put(carry,1); new_line;
     -- put("n2*quot - n1 : "); put(rest); new_line;
      carry := 0;                             -- rest = n2 - (n2*quot-n1)
      for i in n2'range loop                  -- add n2 to the rest
        diff := integer32(n2(i)) - integer32(rest(i)) - integer32(carry);
        if diff >= 0
         then n1(i) := natural32(diff);                       carry := 0;
         else n1(i) := natural32(integer32(the_base) + diff); carry := 1;
        end if;
      end loop;
      for i in n2'last+1..n1'last loop
        n1(i) := rest(i);
      end loop;
    else -- new_line;
      n1 := rest;
    end if;
   -- put("quot : "); put(quot,1); new_line;
   -- put("n1 : "); put(n1); new_line;
    q := quot;
  end Aligned_Acc_Div;

  procedure Uneven_Acc_Div
               ( n1 : in out Array_of_Naturals; leadn1 : in natural32;
                 n2 : in Array_of_Naturals; q,shifts : out natural32 ) is

  -- DESCRIPTION :
  --   This accumulated division is auxiliary in the division routine of
  --   two long coefficient vectors.  This version assumes that n1 and n2
  --   are so that n1(leadn1) < n2(n2'last).

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

    head,tail,inc,shiftctr,quot,acc,carry,lostprod : natural32;
    prod : Array_of_Naturals(n2'range);
    rest : Array_of_Naturals(n1'range);
    diff : integer32;

  begin
   -- put_line("Entered Uneven Acc Div");
    head := n1(leadn1);                   -- head will gain digits from tail
    tail := n1(leadn1-1);                 --        to become >= n2(n2'last)
    shiftctr := 0;                        -- counts the number of shifts
    loop
      shiftctr := shiftctr + 1;
      inc := tail/the_base;               -- increment is first digit of tail
      head := head*fact + inc;            -- shift and add increment to head
      exit when (head >= n2(n2'last));
      tail := tail - inc*the_base;        -- remove increment from tail
      tail := tail*fact;                  -- shift tail up
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
   -- put("     n1 : "); put(n1); new_line;
   -- put("     n2 : "); put(n2); new_line;
   -- put("n2*quot : "); put(lostprod); put(prod);
    inc := leadn1-n2'last;                       -- subtract n2*quot from n1
   -- put("  inc : "); put(inc,1); new_line;
    rest := n1;
    diff := integer32(n1(inc-1)) - integer32(lostprod);
    if diff >= 0
     then rest(inc-1) := natural32(diff);                       carry := 0;
     else rest(inc-1) := natural32(diff + integer32(the_base)); carry := 1;
    end if;
    for i in n2'range loop
      diff := integer32(n1(i+inc)) - integer32(prod(i)) - integer32(carry);
      if diff >= 0
       then rest(i+inc) := natural32(diff);                       carry := 0;
       else rest(i+inc) := natural32(integer32(the_base) + diff); carry := 1;
      end if; 
    end loop;
   -- put("diff carry : "); put(carry,1);
    if carry > 0 then
     -- put_line("   OVERFLOW in Uneven_Acc_Div");
      quot := quot - 1;    -- decrease quot and compute n2*quot - n1
      diff := integer32(lostprod) - integer32(n1(inc-1)); 
      if diff >= 0
       then rest(inc-1) := natural32(diff);                       carry := 0;
       else rest(inc-1) := natural32(diff + integer32(the_base)); carry := 1;
      end if;
      for i in n2'range loop
        diff := integer32(prod(i)) - integer32(n1(i+inc)) - integer32(carry);
        if diff >= 0
         then rest(i+inc) := natural32(diff);                       carry := 0;
         else rest(i+inc) := natural32(integer32(the_base) + diff); carry := 1;
        end if; 
      end loop;
     -- put("diff carry : "); put(carry,1); new_line;
     -- put("n2*quot - n1 : "); put(rest); new_line;
      carry := 0;                             -- rest = n2 - (n2*quot-n1)
      for i in n2'range loop                  -- add n2 to the rest
        diff := integer32(n2(i)) - integer32(rest(i)) - integer32(carry);
        if diff >= 0
         then n1(i) := natural32(diff);                       carry := 0;
         else n1(i) := natural32(integer32(the_base) + diff); carry := 1;
         end if;
      end loop;
      for i in n2'last+1..n1'last loop
        n1(i) := rest(i);
      end loop;
    else -- new_line;
      n1 := rest;
    end if;
   -- put("quot : "); put(quot,1); new_line;
   -- put("n1      : ");  put(n1); new_line;
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
    dom := false;   -- correct default ?
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
    dn1size,multnorm,scn2 : natural32;
    quot : Array_of_Naturals(q'range) := (q'range => 0);

  begin
   -- put_line("Entering Div");
   -- put(" n1 : "); put(n1); new_line;
   -- put(" n2 : "); put(n2); new_line;
    dn1(n1'range) := n1;                             -- dn1 is shifted n1
    dn1(dn1'last) := 0;
    dn2 := n2;
   -- multnorm := 0;
    scn2 := Digits_to_Normal(dn2);
   -- while dn2(dn2'last) < sub_base loop              -- dn2 is normalized n2
   --   Mul_Fact(dn1,fact);
   --   Mul_Fact(dn2,fact);
   --   multnorm := multnorm+1;
   -- end loop;
    if scn2 /= 0 then
      Digits_to_Left(dn1,scn2);
      Digits_to_Left(dn2,scn2);
    end if;
   -- put("multnorm = "); put(multnorm,1);
   -- put("  digits to left = "); put(scn2,1); new_line;
   -- if multnorm /= scn2 then raise CONSTRAINT_ERROR; end if;
    multnorm := scn2;
    if dn1(dn1'last) = 0                             -- determine work size
     then dn1size := n1'last;                        -- of the rest of n1
     else dn1size := dn1'last;
    end if;
   -- put_line("In the Div on Array_of_Naturals : ");
   -- put("dn1 : "); put(dn1); new_line;
   -- put("dn2 : "); put(dn2); new_line;
   -- put("multnorm : "); put(multnorm,1); new_line;
    declare
      rest : Array_of_Naturals(0..dn1size);
      leadn1,quotfac,quotind,sum,carry,shifts : natural32;
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
          Dominating(rest,dn2,leadn1-1,dn2'last-1,dom,eq);
         -- put("Checking dominating : ");
         -- if dom
         --  then put_line("okay.");
         --  else put_line("not okay.");
         -- end if;
          if not dom then
           -- put("rest : "); put(rest); new_line;
           -- put("dn2 : "); put(dn2); new_line;
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
              quotind := quot'last+1;  -- stop the iteration, only if eq !!!
            end if;
           -- put("rest(leadn1)  : "); put(rest(leadn1),1); new_line;
           -- put("dn2(dn2'last) : "); put(dn2(dn2'last),1); new_line;
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
       -- put("accquot : "); put(quot); new_line;
       -- put("n2*quot : "); put(Mul(quot,n2)); new_line;
        while (leadn1 >= dn2'last) and (rest(leadn1) = 0) loop
          leadn1 := leadn1-1;
        end loop;
      end loop;
      if multnorm > 0 then                              -- shift the rest 
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
    n1 := q;
  end Div;

--  function Est_Quot ( n1 : Array_of_Naturals; lead : natural;
--                      n2f : double_float ) return natural is
--
--  -- DESCRIPTION :
--  --   Returns the estimate for the quotient as subroutine of Div.
--
--  -- REQUIRED : n1(0..lead) > n2.
--
--  -- ON ENTRY :
--  --   n1       numerator, number to be divided with n1(lead) /= 0;
--  --   lead     largest index in n1 so n1(lead) /= 0;
--  --   n2       denominator, only n2(1) and n2(0) play a role;
--  --   n2f      representation of n2 as a double float.
--
--    res : natural;
--    f1 : constant double_float := Double(n1,lead);
--    r : constant double_float := f1/n2f;
--
--  begin
--    res := natural(r);
--    if r < double_float(res) then
--      if res > 0
--       then res := res - 1;
--      end if;
--    end if;
--    return res;
--  end Est_Quot;

  procedure Make_Prod ( p : in out Array_of_Naturals;
                        q : in natural64; n2 : in Array_of_Naturals ) is

  -- DESCRIPTION :
  --   Computes in p the product of q with n2.

  -- REQUIRED : p'range = 0..n2'last+1, q < Base.

    prod,carry : natural64;

  begin
    carry := 0;
    for i in n2'range loop
      prod := natural64(n2(i))*q + carry;
      p(i) := natural32(prod mod base64);
      carry := prod/base64;
    end loop;
    p(p'last) := natural32(carry);
  end Make_Prod;

  function Alignment ( w : Array_of_Naturals; lead : natural32;
                       p : Array_of_Naturals ) return integer32 is

  -- DESCRIPTION :
  --   Determines the alignment of p with respect to w
  --   before subtracting p from w, auxiliary for Div2.

  -- REQUIRED : w(lead) /= 0 and lead > p'last.

  begin
   -- put("In alignment, w'last = "); put(w'last,1);
   -- put(", lead = "); put(lead,1);
   -- put(", p'last = "); put(p'last,1);
   -- put(", p(p'last) = "); put(p(p'last),1); new_line;
   -- put("w = "); put(w); new_line;
   -- put("p = "); put(p); new_line;
    if p(p'last) /= 0 then
      return integer32(lead) - integer32(p'last);
    elsif w(lead) < p(p'last-1) then
      return integer32(lead) - integer32(p'last);
    else
      return integer32(lead) - integer32(p'last) + 1;
    end if;
  end Alignment;

  procedure Sub_Div ( w : in out Array_of_Naturals;
                      p : in Array_of_Naturals; least : in natural32;
                      fail : out boolean ) is

  -- DESCRIPTION :
  --   Subtracts from w the product of quotient times divisor in p.
  --   This is a subroutine of Div2.
  --   If w - p turns out negative, then fail is true 
  --   and the estimate for the quotient was too high.

  -- REQUIRED :
  --   w(lead) is highest nonzero number in w.

  -- ON ENTRY :
  --   w         accumulated remainder;
  --   p         product of quotient with divisor;
  --   least     result of Alignment(w,lead,p).

  -- ON RETURN :
  --   w         result of w - p, provided fail is false;
  --   fail      true if w - p < 0, false otherwise.

    backup : constant Array_of_Naturals(w'range) := w;
    diff : integer32;
    carry : natural32 := 0;

  begin
   -- put("Sub w = "); put(w); new_line;
   -- put("  - p = "); put(p); new_line;
    fail := false;
    for i in 0..p'last-1 loop
      diff := integer32(w(least+i)) - integer32(p(i)) - integer32(carry);
      if diff >= 0 then
        w(least+i) := natural32(diff);
        carry := 0;
      else
        w(least+i) := natural32(integer32(the_base) + diff);
        carry := 1;
      end if;
    end loop;
    if least+p'last <= w'last then
      diff := integer32(w(least+p'last))
               - integer32(p(p'last)) - integer32(carry);
      if diff >= 0 then
        w(least+p'last) := natural32(diff);
      else
        fail := true;
      end if;
    elsif p(p'last) /= 0 or carry /= 0 then
      fail := true;
    end if;
    if fail
     then w := backup;
    end if;
   -- put("Aft w = "); put(w); new_line;
  end Sub_Div;

  procedure Div2 ( n1,n2 : in Array_of_Naturals;
                   q,r : out Array_of_Naturals ) is

    leadn1 : natural32 := n1'last;
    prevleast : natural32 := n1'last+1;
    qlast : natural32 := 1;
   -- n2f : constant double_float := Double(n2,n2'last);
    n2d : constant natural64 := To_Natural64(n2,n2'last);
    work : Array_of_Naturals(n1'range) := n1;
    quot : natural32;
    least,undershot : integer32;
    prod : Array_of_Naturals(0..n2'last+1);
    overshot : boolean;

  begin
   -- new_line;
   -- put_line("new call of Div2 ...");
   -- put("n1 = "); put(n1); new_line;
   -- put("n2 = "); put(n2); new_line;
    q := (q'range => 0);
    r := (r'range => 0);
    while n1(leadn1) = 0 and leadn1 > 0 loop
      leadn1 := leadn1 - 1;
    end loop;
    if n1(0..leadn1) < n2 then
      r(n1'range) := n1;
    else
      quot := Est_Quot2(work,leadn1,n2d); --Est_Quot2(work,leadn1,n2,n2d);
      undershot := 0;
      for i in 0..10*n1'last loop
       -- put("******** STEP "); put(i,1); put_line(" ***********");
       -- put("estimate for quotient : "); put(quot,1); new_line;
        for i in n1'range loop
          Make_Prod(prod,natural64(quot),n2);
          least := Alignment(work,leadn1,prod);
         -- put("least = "); put(least,1); new_line;
          if least <  0
           then overshot := true;
           else Sub_Div(work,prod,natural32(least),overshot);
          end if;
         -- put("work = "); put(work); new_line;
          exit when not overshot;
         -- put("overshot the quotient, quot = "); put(quot,1); new_line;
         -- put("overshot the quotient, work = "); put(work); new_line;
          if quot > 10 then
            quot := natural32(double_float(quot)*0.75);
         -- elsif quot > 0 then
          else
            quot := quot - 1;
          end if;
          exit when ((quot = 0) and (least = 0));
          if quot = 0 and least > 0 then
            if least <= integer32(prevleast)
             then prevleast := natural32(least);
            end if;
            least := least - 1;
            quot := the_base-1;
            Make_Prod(prod,natural64(quot),n2);
            Sub_Div(work,prod,natural32(least),overshot);
          end if;
          exit when not overshot;
        end loop;
        exit when (quot = 0);
       -- put("prevleast = "); put(prevleast,1);
       -- put(", least = "); put(least,1); new_line;
        if prevleast <= leadn1
         then undershot := integer32(prevleast) - least;
        end if;
        if undershot > 0
         then qlast := qlast + natural32(undershot);
        end if;
        Add_Quot2(q,qlast,quot,undershot);
        while work(leadn1) = 0 and leadn1 > 0 loop
          leadn1 := leadn1 - 1;
        end loop;
        exit when (work(leadn1) = 0);
        exit when (work(0..leadn1) < n2);
        if least <= integer32(prevleast)
         then prevleast := natural32(least);
        end if;
        quot := Est_Quot2(work,leadn1,n2d); --Est_Quot2(work,leadn1,n2,n2d);
      end loop;
      if least > 0 then
        for i in 0..(qlast-1) loop
          exit when (i+natural32(least) > q'last);
          q(i+natural32(least)) := q(i);
        end loop;
        for i in 0..least-1 loop
          q(natural32(i)) := 0;
        end loop;
      end if;
      r := work(r'range);
    end if;
 -- exception
 --   when others => put_line("exception raised in Div2");
 --                  put("n1 = "); put(n1); new_line;
 --                  put("n2 = "); put(n2); new_line; raise;
  end Div2;

  procedure Div2 ( n1 : in out Array_of_Naturals; n2 : in Array_of_Naturals;
                   r : out Array_of_Naturals ) is

    q : Array_of_Naturals(n1'range);

  begin
    Div2(n1,n2,q,r);
    n1 := q;
  end Div2;

end Multprec_Natural_Coefficients;
