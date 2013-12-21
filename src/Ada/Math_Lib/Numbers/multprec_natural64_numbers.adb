with unchecked_deallocation;
--with text_io,integer_io;  use text_io,integer_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;

package body Multprec_Natural64_Numbers is

-- NOTES ON THE CHOICE OF REPRESENTATION AND IMPLEMENTATION :
--   1) The base is decimal, for convenience of input/output.
--      Everything should remain correct if another base is chosen,
--      except of course the function "Decimal_Places".
--      The base is such that factor times the base is still a standard
--      natural number, so that standard arithmetic can be used.
--   2) A natural number is implemented by means of a coefficient array.
--      The advantage over a linked-list representation is that the
--      coefficients can be accessed faster.  The disadvantage of some
--      wasted storage is countered by the aim of implementing higher
--      precision floating-point numbers, where the mantissa is fixed.
--   3) The usage of a pointer to access the coefficient vector permits
--      to have uniformity in input/output, since on input it is not
--      known in advance how long a natural number may be.
--   4) A natural number n for which Empty(n) holds is considered as zero.
--   5) In producing the result of a binary operation, except in cases
--      where the result is zero, the result should be represented as
--      a number of the size equal to that of the operands.  This is to
--      enable accurate mixed precision floating-point arithmetic.

  fact : constant natural64 := Multprec_Natural64_Coefficients.Radix;
  expo : constant natural32 := Multprec_Natural64_Coefficients.Exponent;
  sqrt_base : constant natural64 := fact**(natural(expo/2));
  the_base : constant natural64 := sqrt_base*sqrt_base;

-- DATASTRUCTURES :

  type Natural_Number_Rep ( size : natural32 ) is record
    numb : Array_of_Naturals(0..size);
  end record;

-- AUXILIARY :

  function Size ( n : natural64 ) return natural32 is

  -- DESCRIPTION :
  --   Determines the size of the natural number to represent n.

    acc : natural64 := the_base;

  begin
    for i in 0..n loop
      if n < acc
       then return natural32(i);
       else acc := acc*the_base;
      end if;
    end loop;
    return natural32(n);
  end Size;

-- CREATORS :

  function Create ( n : natural64 ) return Array_of_Naturals is

    res : Array_of_Naturals(0..Size(n));
    remainder : natural64 := n;

  begin
    for i in res'range loop
      res(i) := remainder mod the_base;
      remainder := remainder/the_base;
    end loop;
    return res;
  end Create;

  function Create ( n : natural64 ) return Natural_Number is

    res : constant Natural_Number := Create(Create(n));

  begin
    return res;
  end Create;

  function Create ( n : Array_of_Naturals ) return Natural_Number is

    res : Natural_Number;
    res_rep : Natural_Number_Rep(n'last);

  begin
    res_rep.numb := n;
    res := new Natural_Number_Rep'(res_rep);
    return res;
  end Create;

  function Create ( n : Natural_Number ) return natural64 is

    res : natural64;

  begin
    if Empty(n) then
      res := 0;
    else
      res := Coefficient(n,Size(n));
      for i in reverse 0..Size(n)-1 loop
        res := res*the_base + Coefficient(n,i);
      end loop;
    end if;
    return res;
  end Create;

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

  function Empty ( n : Natural_Number ) return boolean is
  begin
    return ( n = null );
  end Empty;

  function Size ( n : Natural_Number ) return natural32 is
  begin
    if n = null
     then return 0;
     else return n.size;
    end if;
  end Size;

  function Decimal_Places ( n : natural64 ) return natural32 is

    acc : natural64 := 1;

  begin
    for i in 0..n loop
      if n < acc
       then return natural32(i);
       else acc := acc*10;
      end if;
    end loop;
    return natural32(n);
  end Decimal_Places;

  function Decimal_Places ( n : Natural_Number ) return natural32 is
  
    res : natural32 := 0;

  begin
    if not Empty(n) then
      for i in reverse 0..Size(n) loop
        if n.numb(i) /= 0
         then res := Decimal_Places(n.numb(i)) + i*expo; exit;
        end if;
      end loop;
    end if;
    return res;
  end Decimal_Places;

  function Coefficient ( n : Natural_Number; i : natural32 ) return natural64 is
  begin
    if (n = null or else i > n.size)
     then return 0;
     else return n.numb(i);
    end if;
  end Coefficient;

  function Coefficients ( n : Natural_Number ) return Array_of_Naturals is

    na : constant Array_of_Naturals(0..0) := (0..0 => 0);

  begin
    if Empty(n)
     then return na;
     else return n.numb;
    end if;
  end Coefficients;

-- COMPARISON AND COPYING :
--   Note that these implementations are all independent of the
--   actual representation of Natural_Number as they use the selectors.

  function Equal ( n1 : Natural_Number; n2 : natural64 ) return boolean is
  begin
    if Empty(n1) then
      if n2 = 0
       then return true;
       else return false;
      end if;
    else
      declare
        n2cff : constant Array_of_Naturals := Create(n2);
      begin
        if n2cff'last > Size(n1) then
          return false;
        else
          for i in n2cff'range loop
            if Coefficient(n1,i) /= n2cff(i)
             then return false;
            end if;
          end loop;
          for i in n2cff'last+1..Size(n1) loop
            if Coefficient(n1,i) /= 0
             then return false;
            end if;
          end loop;
       end if;
      end;
      return true;
    end if;
  end Equal;

  function Equal ( n1,n2 : Natural_Number ) return boolean is

    min_size : natural32;

  begin
    if Empty(n1) then
      return Equal(n2,0);
    else
      if Empty(n2) then
        return Equal(n1,0);
      else
        if Size(n1) < Size(n2) then
          for i in Size(n1)+1..Size(n2) loop
            if Coefficient(n2,i) /= 0
             then return false;
            end if;
          end loop;
          min_size := Size(n1);
        elsif Size(n1) > Size(n2) then
          for i in Size(n2)+1..Size(n1) loop
            if Coefficient(n1,i) /= 0
             then return false;
            end if;
          end loop;
          min_size := Size(n2);
        else
          min_size := Size(n1);
        end if;
        return Equal(Coefficients(n1)(0..min_size),
                     Coefficients(n2)(0..min_size));
      end if;
    end if;
  end Equal;

  function "<" ( n1 : Natural_Number; n2 : natural64 ) return boolean is
  begin
    if Empty(n1) then
      if n2 > 0
       then return true;
       else return false;
      end if;
    else
      declare
        n2cff : constant Array_of_Naturals := Create(n2);
      begin
        if n2cff'last > Size(n1) then
          return true;
        else
          if n2cff'last < Size(n1) then
            for i in reverse n2cff'last+1..Size(n1) loop
              if Coefficient(n1,i) /= 0
               then return false;
              end if;
            end loop;
          end if;
          for i in reverse n2cff'range loop
            if Coefficient(n1,i) >= n2cff(i)
             then return false;
            end if;
          end loop;
        end if;
      end;
      return true;
    end if;
  end "<";

  function "<" ( n1 : natural64; n2 : Natural_Number ) return boolean is
  begin
    if Empty(n2) then
      return false;
    else
      declare
        n1cff : constant Array_of_Naturals := Create(n1);
      begin
        if n1cff'last > Size(n2) then
          return false;
        else
          if n1cff'last < Size(n2) then
            for i in reverse n1cff'last+1..Size(n2) loop
              if Coefficient(n2,i) /= 0
               then return true;
              end if;
            end loop;
          end if;
          for i in reverse n1cff'range loop
            if n1cff(i) >= Coefficient(n2,i)
             then return false;
            end if;
          end loop;
        end if;
      end;
      return true;
    end if;
  end "<";

  function "<" ( n1,n2 : Natural_Number ) return boolean is

    min_size : natural32;

  begin
    if Empty(n1) then
      if Empty(n2)
       then return false;
       else return true;
      end if;
    elsif Empty(n2) then
      return false;
    else
      if Size(n1) < Size(n2) then
        for i in Size(n1)+1..Size(n2) loop
          if Coefficient(n2,i) /= 0
           then return true;
          end if;
        end loop;
        min_size := Size(n1);
      elsif Size(n1) > Size(n2) then
        for i in Size(n2)+1..Size(n1) loop
          if Coefficient(n1,i) /= 0
           then return false;
          end if;
        end loop;
        min_size := Size(n2);
      else
        min_size := Size(n1);
      end if;
      for i in reverse 0..min_size loop
        if Coefficient(n1,i) < Coefficient(n2,i) then
          return true;
        elsif Coefficient(n1,i) > Coefficient(n2,i) then
          return false;
        end if;
      end loop;
      return false;      -- n1 = n2
    end if;
  end "<";

  function ">" ( n1 : Natural_Number; n2 : natural64 ) return boolean is
  begin
    if Empty(n1) then
      return false;
    else
      declare
        n2cff : constant Array_of_Naturals := Create(n2);
      begin
        if n2cff'last > Size(n1) then
          return false;
        else
          if n2cff'last < Size(n1) then
            for i in n2cff'last+1..Size(n1) loop
              if Coefficient(n1,i) /= 0
               then return true;
              end if;
            end loop;
          end if;
          for i in reverse n2cff'range loop
            if Coefficient(n1,i) <= n2cff(i)
             then return false;
            end if;
          end loop;
        end if;
      end;
      return true;
    end if;
  end ">";
 
  function ">" ( n1 : natural64; n2 : Natural_Number ) return boolean is
  begin
    if Empty(n2) then
      if n1 > 0
       then return true;
       else return false;
      end if;
    else
      declare
        n1cff : constant Array_of_Naturals := Create(n1);
      begin
        if n1cff'last > Size(n2) then
          return true;
        else
          if n1cff'last < Size(n2) then
            for i in n1cff'last+1..Size(n2) loop
              if Coefficient(n2,i) /= 0
               then return false;
              end if;
            end loop;
          end if;
          for i in reverse n1cff'range loop
            if n1cff(i) <= Coefficient(n2,i)
             then return false;
            end if;
          end loop;
        end if;
      end;
      return true;
    end if;
  end ">";

  function ">" ( n1,n2 : Natural_Number ) return boolean is

    min_size : natural32;

  begin
    if Empty(n2) then
      if Empty(n1)
       then return false;
       else return true;
      end if;
    elsif Empty(n1) then
      return false;
    else
      if Size(n1) < Size(n2) then
        for i in Size(n1)+1..Size(n2) loop
          if Coefficient(n2,i) /= 0
           then return false;
          end if;
        end loop;
        min_size := Size(n1);
      elsif Size(n1) > Size(n2) then
        for i in Size(n2)+1..Size(n1) loop
          if Coefficient(n1,i) /= 0
           then return true;
          end if;
        end loop;
        min_size := Size(n2);
      else
        min_size := Size(n1);
      end if;
      for i in reverse 0..min_size loop
        if Coefficient(n1,i) > Coefficient(n2,i) then
          return true;
        elsif Coefficient(n1,i) < Coefficient(n2,i) then
          return false;
        end if;
      end loop;
      return false;  -- n1 = n2
    end if;
  end ">";

  procedure Copy ( n1 : in natural64; n2 : in out Natural_Number ) is
  begin
    Clear(n2);
    n2 := Create(n1);
  end Copy;

  procedure Copy ( n1 : in Natural_Number; n2 : in out Natural_Number ) is
  begin
    Clear(n2);
    if not Empty(n1) then
      declare
        n2rep : Natural_Number_Rep(n1.size);
      begin
        for i in n2rep.numb'range loop
          n2rep.numb(i) := n1.numb(i);
        end loop;
        n2 := new Natural_Number_Rep'(n2rep);
      end;
    end if;
  end Copy;

-- SHIFTS :

  procedure Shift_Left ( n : in out Natural_Number; ns : out natural32 ) is
  begin
    if n = null
     then ns := 0;
     else Shift_Left(n.numb,ns);
    end if;
  end Shift_Left;

  procedure Shift_Right ( n : in out Natural_Number; ns : out natural32 ) is
  begin
    if n = null
     then ns := 0;
     else Shift_Right(n.numb,ns);
    end if;
  end Shift_Right;

-- AUXILIARIES FOR ARITHMETIC OPERATIONS :

  function Extend ( n : Array_of_Naturals; carry : natural64 )
                  return Natural_Number_Rep is

  -- DESCRIPTION :
  --   Extends the array of naturals as much as need to store the carry
  --   obtained from adding a number to a natural number.

  begin
    if carry = 0 then
      declare
        nrep : Natural_Number_Rep(n'last);
      begin
        nrep.numb := n;
        return nrep;
      end;
    else
      declare
        np1 : Array_of_Naturals(n'first..n'last+1);
        carry1 : natural64;
      begin
        np1(n'range) := n;
        np1(np1'last) := carry mod the_base;
        carry1 := carry/the_base;
        return Extend(np1,carry1);
      end;
    end if;
  end Extend;

  procedure Extend ( n : in out Natural_Number; carry : in natural64 ) is

  -- DESCRIPTION :
  --   Extends the coefficient vector of the natural number to store
  --   the carry-over.
 
  begin
    if carry /= 0 then
      declare
        res : constant Natural_Number
            := new Natural_Number_Rep'(Extend(n.numb,carry));
      begin
        Clear(n);
        n := res;
      end;
    end if;
  end Extend;

 -- procedure Mul_Fact ( n : in out Natural_Number; f : in natural ) is

  -- DESCRIPTION :
  --   Multiplies the number n with 0 < f <= fact and stores the result in n.

 --   prod,carry : natural;

 -- begin
 --   if not Empty(n) then
 --     carry := 0;
 --     for i in n.numb'range loop
 --       prod := n.numb(i)*f + carry;
 --       n.numb(i) := prod mod the_base;
 --       carry := prod/the_base;
 --     end loop;
 --     Extend(n,carry);
 --   end if;
 -- end Mul_Fact;

  procedure Mul_Base ( n : in out Natural_Number; k : in natural32 ) is

  -- DESCRIPTION :
  --   Multiplies the number n with the_base**k.

  begin
    if k > 0 and not Empty(n) then
      declare
        npk : Natural_Number_Rep(n.size+k);
      begin
        npk.numb(0..k-1) := (0..k-1 => 0);
        npk.numb(k..npk.numb'last) := n.numb(n.numb'range);
        Clear(n);
        n := new Natural_Number_Rep'(npk);
      end;
    end if;
  end Mul_Base;

--  procedure Div_Base ( n : in out Natural_Number; k : in natural ) is

  -- DESCRIPTION :
  --   Divides the number n by the_base**k.

--  begin
--    if k > 0 and not Empty(n) then
--      if k > n.size then
--        Clear(n);
--      else
--        declare
--          nmk : Natural_Number_Rep(n.size-k);
--        begin
--          for i in nmk.numb'range loop
--            nmk.numb(i) := n.numb(k+i);
--          end loop;
--          Clear(n);
--          n := new Natural_Number_Rep'(nmk);
--        end;
--      end if;
--    end if;
--  end Div_Base;

  procedure Mul_Rdx ( n : in out Natural_Number; k : in natural32 ) is

  -- DESCRIPTION : 
  --   Multiplies n with the radix, assuming k < expo and n /= 0.

    multfac : constant natural64 := fact**natural(k);
    divisor : constant natural64 := fact**natural(expo-k);
    rest,quot : natural64;
    resnumb : Array_of_Naturals(0..Size(n)+1);

  begin
    quot := 0;
    for i in 0..Size(n) loop
      rest := n.numb(i) mod divisor;
      resnumb(i) := rest*multfac + quot;
      quot := n.numb(i)/divisor;
    end loop;
    if quot = 0 then
      n.numb := resnumb(0..Size(n));
    else
      Clear(n);
      resnumb(resnumb'last) := quot;
      n := Create(resnumb);
    end if;
  end Mul_Rdx;

  procedure Mul_Radix ( n : in out Natural_Number; k : in natural32 ) is

  -- DESCRIPTION :
  --   Multiplies n k times with the radix.

    q : natural32;

  begin
    if not Equal(n,0) then
      if k < expo then
        Mul_Rdx(n,k);
      else
        Mul_Rdx(n,(k mod expo));
        q := k/expo;
        if q > 0
         then Mul_Base(n,q);    
        end if;
      end if;
    end if;
  end Mul_Radix;

-- ARITHMETIC OPERATIONS :

  function "+" ( n1 : Natural_Number; n2 : natural64 ) return Natural_Number is

    res : Natural_Number;

  begin
    if Empty(n1) then
      res := Create(n2);
    else
      declare
        cff : Array_of_Naturals(n1.numb'range);
        sum,carry : natural64;
      begin
        carry := n2;
        for i in cff'range loop
          sum := n1.numb(i) + carry;
          cff(i) := sum mod the_base;
          carry := sum/the_base;
        end loop;
        res := new Natural_Number_Rep'(Extend(cff,carry));
      end;
    end if;
    return res;
  end "+";

  function "+" ( n1 : natural64; n2 : Natural_Number ) return Natural_Number is
  begin
    return n2+n1;
  end "+";

  function "+" ( n1,n2 : Natural_Number ) return Natural_Number is

    res : Natural_Number;
    sum,carry : natural64;

  begin
    if Empty(n1) then
      Copy(n2,res);
    elsif Empty(n2) then
      Copy(n1,res);
    else
      carry := 0;
      if n1.size > n2.size then
        declare
          res_rep : Natural_Number_Rep(n1.size) := n1.all;
        begin
          for i in n2.numb'range loop
            sum := res_rep.numb(i) + n2.numb(i) + carry;
            res_rep.numb(i) := sum mod the_base;
            carry := sum/the_base;
          end loop;
          if carry /= 0 then
            for i in n2.numb'last+1..n1.numb'last loop
              sum := res_rep.numb(i) + carry;
              res_rep.numb(i) := sum mod the_base;
              carry := sum/the_base;
              exit when (carry = 0);
            end loop;
          end if;
          res := new Natural_Number_Rep'(res_rep);
        end;
      else
        declare
          res_rep : Natural_Number_Rep(n2.size) := n2.all;
        begin
          for i in n1.numb'range loop
            sum := res_rep.numb(i) + n1.numb(i) + carry;
            res_rep.numb(i) := sum mod the_base;
            carry := sum/the_base;
          end loop;
          if carry /= 0 then
            for i in n1.numb'last+1..n2.numb'last loop
              sum := res_rep.numb(i) + carry;
              res_rep.numb(i) := sum mod the_base;
              carry := sum/the_base;
              exit when (carry = 0);
            end loop;
          end if;
          res := new Natural_Number_Rep'(res_rep);
        end;
      end if;
      Extend(res,carry);
    end if;
    return res;
  end "+";

  procedure Add ( n1 : in out Natural_Number; n2 : in natural64 ) is

    sum,carry : natural64;

  begin
    if Empty(n1) then
      n1 := Create(n2);
    else
      carry := n2;
      for i in n1.numb'range loop
        sum := n1.numb(i) + carry;
        n1.numb(i) := sum mod the_base;
        carry := sum/the_base;
      end loop;
      Extend(n1,carry);
    end if;
  end Add;

  procedure Add ( n1 : in out Natural_Number; n2 : in Natural_Number ) is

    sum,carry : natural64;

  begin
    if Empty(n1) then
      Copy(n1 => n2,n2 => n1);
    elsif not Empty(n2) then
      carry := 0;
      if n1.size >= n2.size then
        for i in n2.numb'range loop
          sum := n1.numb(i) + n2.numb(i) + carry;
          n1.numb(i) := sum mod the_base;
          carry := sum/the_base;
        end loop;
        if carry /= 0 then
          for i in n2.numb'last+1..n1.numb'last loop
            sum := n1.numb(i) + carry;
            n1.numb(i) := sum mod the_base;
            carry := sum/the_base;
            exit when (carry = 0);
          end loop;
          Extend(n1,carry);
        end if;
      else
        declare
          res_rep : Natural_Number_Rep(n2.size) := n2.all;
        begin
          for i in n1.numb'range loop
            sum := res_rep.numb(i) + n1.numb(i) + carry;
            res_rep.numb(i) := sum mod the_base;
            carry := sum/the_base;
          end loop;
          if carry /= 0 then
            for i in n1.numb'last+1..n2.numb'last loop
              sum := res_rep.numb(i) + carry;
              res_rep.numb(i) := sum mod the_base;
              carry := sum/the_base;
              exit when (carry = 0);
            end loop;
          end if;
          Clear(n1);
          n1 := new Natural_Number_Rep'(res_rep);
          if carry /= 0
           then Extend(n1,carry);
          end if;
        end;
      end if;
    end if;
  end Add;

  function "-" ( n1 : Natural_Number; n2 : natural64 ) return Natural_Number is

    res : Natural_Number;
    diff,carry : integer64;

  begin
    if not Empty(n1) then
      Copy(n1,res);
      declare
        n2cff : constant Array_of_Naturals := Create(n2);
        index : natural32 := n2cff'first;
      begin
        carry := integer64(n2cff(index));
        for i in n1.numb'range loop
          diff := integer64(n1.numb(i)) - carry;
          if diff >= 0 then
            res.numb(i) := natural64(diff mod integer64(the_base));
            carry := diff/integer64(the_base);
          else
            diff := integer64(the_base) + diff;
            res.numb(i) := natural64(diff mod integer64(the_base));
            carry := 1;
          end if;
          if index < n2cff'last then
            index := index+1;
            carry := carry + integer64(n2cff(index));
          end if;
          exit when (carry = 0);
        end loop;
      end;
    end if;
    return res;
  end "-";

  function "-" ( n1 : natural64; n2 : Natural_Number ) return Natural_Number is
  
    tmp : Natural_Number := Create(n1);
    res : constant Natural_Number := tmp - n2;

  begin
    Clear(tmp);
    return res;
  end "-";

  function "-" ( n1,n2 : Natural_Number ) return Natural_Number is

    res : Natural_Number;

  begin
    Copy(n1,res);
    if not Empty(n2)
     then Sub(res.numb,n2.numb);
    end if;
    return res;
  end "-";

  function "+" ( n : Natural_Number ) return Natural_Number is

    res : Natural_Number;

  begin
    Copy(n,res);
    return res;
  end "+";

  function "-" ( n : Natural_Number ) return Natural_Number is

    res : Natural_Number;

  begin
    Copy(n,res);
    return res;
  end "-";

  procedure Min ( n : in out Natural_Number ) is
  begin
    null;
  end Min;

  procedure Sub ( n1 : in out Natural_Number; n2 : in natural64 ) is

    diff,carry : integer64;

  begin
    if not Empty(n1) and n2 > 0 then
      declare
        n2cff : constant Array_of_Naturals := Create(n2);
        index : natural32 := n2cff'first;
      begin
        carry := integer64(n2cff(index));
        for i in n1.numb'range loop
          diff := integer64(n1.numb(i)) - carry;
          if diff >= 0 then
            n1.numb(i) := natural64(diff);
            carry := 0;
          else
            diff := integer64(the_base) + diff;
            n1.numb(i) := natural64(diff mod integer64(the_base));
            carry := 1;
          end if;
          if index < n2cff'last then
            index := index+1;
            carry := carry + integer64(n2cff(index));
          end if;
          exit when (carry = 0);
        end loop;
      end;
    end if;
  end Sub;

  procedure Sub ( n1 : in out Natural_Number; n2 : in Natural_Number ) is
  begin
    if not Empty(n1) and not Empty(n2)
     then Sub(n1.numb,n2.numb);
    end if;
  end Sub;

  function "*" ( n1 : Natural_Number; n2 : natural64 ) return Natural_Number is

    res : Natural_Number;
    acc0,acc1,n10,n11,n20,n21,prod,sum0,sum1,cr1 : natural64;
    sizeres : natural32;
    resnumb : Array_of_Naturals(0..Size(n1)+2);

  begin
    if n2 /= 0 and not Equal(n1,0) then
      n20 := n2 mod sqrt_base;
      n21 := n2/sqrt_base;
      resnumb := (resnumb'range => 0);
      for i in n1.numb'range loop
        n10 := n1.numb(i) mod sqrt_base;
        n11 := n1.numb(i)/sqrt_base;
        acc1 := n21*n11;
        acc0 := n20*n10;
        Acc_Add(resnumb,acc1,acc0,i,0);
        prod := n21*n10;
        if prod < sqrt_base
         then acc0 := prod*sqrt_base;
              acc1 := 0;
         else acc0 := (prod mod sqrt_base)*sqrt_base;
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
            if sum1 < the_base then
              acc1 := sum1;
            else
              acc1 := sum1 mod the_base;
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
        Acc_Add(resnumb,acc1,acc0,i,cr1);
      end loop;
      sizeres := 0;
      for i in reverse resnumb'range loop
        if resnumb(i) /= 0
         then sizeres := i;
        end if;
        exit when (sizeres /= 0);
      end loop;
      if sizeres < Size(n1)        -- do not truncate !!
       then sizeres := Size(n1);
      end if;
      res := Create(resnumb(0..sizeres));
    end if;
    return res;
  end "*";

  function "*" ( n1 : natural64; n2 : Natural_Number ) return Natural_Number is
  begin
    return n2*n1;
  end "*";

  function "*" ( n1,n2 : Natural_Number ) return Natural_Number is

    res : Natural_Number;

  begin
    if (Empty(n1) or else Empty(n2)) then
      return res;
    elsif n1.size >= n2.size then
      res := Create(Mul(n1.numb,n2.numb));
    else
      res := Create(Mul(n2.numb,n1.numb));
    end if;
    return res;
  end "*";

  procedure Mul ( n1 : in out Natural_Number; n2 : in natural64 ) is

    acc0,acc1,n10,n11,n20,n21,prod,sum0,sum1,cr1 : natural64;
    sizeres : natural32;
    resnumb : Array_of_Naturals(0..Size(n1)+2);

  begin
    if n2 = 0 then
      Clear(n1);
    elsif not Equal(n1,0) then
      n20 := n2 mod sqrt_base; n21 := n2/sqrt_base;
      resnumb := (resnumb'range => 0);
      for i in n1.numb'range loop
        n10 := n1.numb(i) mod sqrt_base;
        n11 := n1.numb(i)/sqrt_base;
        acc1 := n21*n11; acc0 := n20*n10;
        Acc_Add(resnumb,acc1,acc0,i,0);
        prod := n21*n10;
        if prod < sqrt_base
         then acc0 := prod*sqrt_base;                 acc1 := 0;
         else acc0 := (prod mod sqrt_base)*sqrt_base; acc1 := prod/sqrt_base;
        end if;
        cr1 := 0;
        prod := n20*n11;
        if prod < sqrt_base then
          sum0 := acc0 + prod*sqrt_base;
          if sum0 < the_base then
            acc0 := sum0;
          else
            acc0 := sum0 mod the_base; sum1 := acc1 + sum0/the_base;
            if sum1 < the_base
             then acc1 := sum1;
             else acc1 := sum1 mod the_base; cr1 := sum1/the_base;
            end if;
          end if;
        else
          sum0 := acc0 + (prod mod sqrt_base)*sqrt_base;
          if sum0 < the_base then
            acc0 := sum0;
          else
            acc0 := sum0 mod the_base; sum1 := acc1 + sum0/the_base;
            if sum1 < the_base
             then acc1 := sum1;
             else acc1 := sum1 mod the_base; cr1 := sum1/the_base;
            end if;
          end if;
          sum1 := acc1 + prod/sqrt_base;
          if sum1 < the_base
           then acc1 := sum1;
           else acc1 := sum1 mod the_base; cr1 := cr1 + sum1/the_base;
          end if;
        end if;
        Acc_Add(resnumb,acc1,acc0,i,cr1);
      end loop;
      sizeres := 0;
      for i in reverse resnumb'range loop
        if resnumb(i) /= 0
         then sizeres := i;
        end if;
        exit when (sizeres /= 0);
      end loop;
      if sizeres < Size(n1)        -- do not truncate !!
       then sizeres := Size(n1);
      end if;
      Clear(n1);
      n1 := Create(resnumb(0..sizeres));
    end if;
  end Mul;

  procedure Mul ( n1 : in out Natural_Number; n2 : in Natural_Number ) is
 
    res : Natural_Number;

  begin
    if not Empty(n1) then
      if Empty(n2) then
        Clear(n1);
      else
        if n1.size >= n2.size
         then res := Create(Mul(n1.numb,n2.numb));
         else res := Create(Mul(n2.numb,n1.numb));
        end if;
        Clear(n1);
        n1 := res;
      end if;
    end if;
  end Mul;

  function "**" ( n1 : Natural_Number;
                  n2 : natural64 ) return Natural_Number is

    res : Natural_Number;

  begin
    if n2 = 0 then
      res := Create(1);
    elsif not Empty(n1) then
      Copy(n1,res);
      for i in 1..n2-1 loop
        Mul(res,n1);
      end loop;
    end if;
    return res;
  end "**";

  function "**" ( n1 : natural64;
                  n2 : Natural_Number ) return Natural_Number is

    res,cnt : Natural_Number;

  begin
    if Equal(n2,0) then
      res := Create(1);
    else
      res := Create(n1);
      if n1 /= 0 then
        cnt := Create(1);
        while not Equal(cnt,n2) loop
          Mul(res,n1);
          Add(cnt,1);
        end loop;
        Clear(cnt);
      end if;
    end if;
    return res;
  end "**";

  function "**" ( n1,n2 : Natural_Number ) return Natural_Number is

    res,cnt : Natural_Number;

  begin
    if Equal(n2,0) then
      res := Create(1);
    elsif not Empty(n1) then
      Copy(n1,res);
      cnt := Create(1);
      while not Equal(cnt,n2) loop
        Mul(res,n1);
        Add(cnt,1);
      end loop;
      Clear(cnt);
    end if;
    return res;
  end "**";

  procedure Small_Div ( n1 : in Natural_Number; n2 : in natural64;
                        q : out Natural_Number; r : out natural64 ) is

  -- DESCRIPTION :
  --   n1 = n2*q + r, only applies when n2 <= fact.

    carry : natural64 := 0;

  begin
    if n2 /= 0 then
      if Equal(n1,0) then
        q := Create(0);
        r := 0;
      else
        q := new Natural_Number_Rep(n1.size);
        for i in reverse 1..n1.size loop
          q.numb(i) := (n1.numb(i)+carry)/n2;
          carry := (n1.numb(i)+carry) mod n2;
          carry := carry*the_base;
        end loop;
        q.numb(0) := (n1.numb(0)+carry)/n2;
        carry := (n1.numb(0)+carry) mod n2;
      end if;
    else
       raise NUMERIC_ERROR;
    end if;
    r := carry;
  end Small_Div;

--  procedure Small_Div ( n1 : in out Natural_Number; n2 : in natural;
--                        r : out natural ) is
--
--    q : Natural_Number;
--
--  begin
--    Small_Div(n1,n2,q,r);
--    Copy(q,n1); Clear(q);
--  end Small_Div;

  procedure Big_Div ( n1 : in Natural_Number; n2 : in natural64;
                      q : out Natural_Number; r : out natural64 ) is

  -- DESCRIPTION :
  --   This division has to be applied when n2 > fact.

  begin
    if n2 = 0 then
      raise Numeric_Error;
    elsif Equal(n1,0) then
      q := Create(0);
      r := 0;
    elsif n1 < n2 then
      q := Create(0);
      r := n1.numb(0);
    else
      declare
        qrep : Array_of_Naturals(n1.numb'range);
      begin
        Big_Div(n1.numb,n2,qrep,r);
        q := Create(qrep);
      end;
    end if;
  end Big_Div;

  function "/" ( n1 : Natural_Number; n2 : natural64 ) return Natural_Number is

    res : Natural_Number;
    r : natural64;

  begin
    if not Empty(n1) then
      if n2 <= fact
       then Small_Div(n1,n2,res,r);
       else Big_Div(n1,n2,res,r);
      end if;
    end if;
    return res;
  end "/";

  function "/" ( n1 : natural64; n2 : Natural_Number ) return natural64 is

    res : natural64;

  begin
    if n1 < n2 then
      res := 0;
    elsif not Empty(n2) then
      res := Create(n2);
      res := n1/res;
    else
      raise NUMERIC_ERROR;
    end if;
    return res;
  end "/";

  function "/" ( n1,n2 : Natural_Number ) return Natural_Number is

    res,rrr : Natural_Number;

  begin 
    Div(n1,n2,res,rrr);
    Clear(rrr);
    return res;
  end "/";

  function Rmd ( n1 : Natural_Number; n2 : natural64 ) return natural64 is

    res : natural64;
    q : Natural_Number;

  begin
    if n2 <= fact
     then Small_Div(n1,n2,q,res);
     else Big_Div(n1,n2,q,res);
    end if;
    Clear(q);
    return res;
  end Rmd;

  function Rmd ( n1 : natural64; n2 : Natural_Number ) return natural64 is

    res : natural64;

  begin
    if n1 < n2 then
      res := n1;
    elsif not Empty(n2) then
      res := Create(n2);
      res := n1 mod res;
    else
      raise NUMERIC_ERROR;
    end if;
    return res;
  end Rmd;

  function Rmd ( n1,n2 : Natural_Number ) return Natural_Number is

    res,q : Natural_Number;

  begin
    Div(n1,n2,q,res);
    Clear(q);
    return res;
  end Rmd;

  procedure Div ( n1 : in out Natural_Number; n2 : in natural64 ) is

    r : natural64;

  begin
    Div(n1,n2,r);
  end Div;

  procedure Div ( n1 : in out Natural_Number; n2 : in Natural_Number ) is

    r : Natural_Number;

  begin
    Div(n1,n2,r);
    Clear(r);
  end Div;

  procedure Div ( n1 : in Natural_Number; n2 : in natural64;
                  q : out Natural_Number; r : out natural64 ) is
  begin
    if n2 <= fact
     then Small_Div(n1,n2,q,r);
     else Big_Div(n1,n2,q,r);
    end if;
  end Div;

  procedure Div ( n1 : in out Natural_Number; n2 : in natural64;
                  r : out natural64 ) is

    q : Natural_Number;

  begin
    Div(n1,n2,q,r);
    Copy(q,n1); Clear(q);
  end Div;

  procedure Div ( n1,n2 : in Natural_Number;
                  q : out Natural_Number; r : out Natural_Number ) is

    leadn2 : natural32;
    rr : natural64;

  begin
    if (Empty(n2) or else Equal(n2,0)) then
      raise NUMERIC_ERROR;
    elsif Equal(n1,0) then
      q := null;
      r := null;
    elsif n1 < n2 then
      q := null;
      Copy(n1,r);
    elsif n2.size = 0 then
      Div(n1,n2.numb(0),q,rr);
      r := Create(rr);
    else 
      leadn2 := 0;             -- leading coefficient in n2
      for i in reverse n2.numb'range loop
        if n2.numb(i) /= 0
         then leadn2 := i; exit;
        end if;
      end loop;
      if leadn2 = 0 then
        Div(n1,n2.numb(0),q,rr);
        r := Create(rr);
      else
        declare
          qrep : Array_of_Naturals(n1.numb'range);
          rest : Array_of_Naturals(0..leadn2);
          sizequot : natural32 := 0;
        begin
         -- put_line("Calling Div on Array_of_Naturals");
         -- put("n1.numb : ");
         -- for i in n1.numb'range loop
         --   put(n1.numb(i));
         -- end loop;
         -- new_line;
         -- put("leadn2 : "); put(leadn2,1); new_line;
         -- put("n2.numb : ");
         -- for i in n2.numb'range loop
         --   put(n2.numb(i));
         -- end loop;
         -- new_line;
          Div(n1.numb,n2.numb(0..leadn2),qrep,rest);
          for i in reverse qrep'range loop
            if qrep(i) /= 0
             then sizequot := i; exit;
            end if;
          end loop;
          q := Create(qrep(0..sizequot));
          r := Create(rest);
        end;
      end if;
    end if;
  end Div;

  procedure Div ( n1 : in out Natural_Number; n2 : in Natural_Number;
                  r : out Natural_Number ) is
  begin
    if (Empty(n2) or else Equal(n2,0)) then
      raise NUMERIC_ERROR;
    elsif not Equal(n1,0) then
      declare
        leadn2 : natural32 := 0;
      begin
        for i in reverse n2.numb'range loop
          if n2.numb(i) /= 0
           then leadn2 := i; exit;
          end if;
        end loop;
        if leadn2 = 0 then
          declare
            rr : natural64;
            sizequot : natural32;
            quot : Array_of_Naturals(n1.numb'range);
          begin
            if n2.numb(0) <= fact
             then Small_Div(n1.numb,n2.numb(0),quot,rr);
             else Big_Div(n1.numb,n2.numb(0),quot,rr);
            end if;
            sizequot := 0;
            for i in reverse quot'range loop
              if quot(i) /= 0
               then sizequot := i; exit;
              end if;
            end loop;
            Clear(n1);
            n1 := Create(quot(0..sizequot));
            r := Create(rr);
          end;
        else
          declare
            rest : Array_of_Naturals(0..leadn2);
            quot : Array_of_Naturals(n1.numb'range);
            sizequot : natural32 := 0;
          begin
           -- put_line("Calling Div on Arrays_of_Naturals");
           -- put("n1.numb : ");
           -- for i in n1.numb'range loop
           --   put(n1.numb(i));
           -- end loop;
           -- new_line;
           -- put("leadn2 : "); put(leadn2,1); new_line;
           -- put("n2.numb : ");
           -- for i in n2.numb'range loop
           --   put(n2.numb(i));
           -- end loop;
           -- new_line;
            Div(n1.numb,n2.numb(0..leadn2),quot,rest);
            for i in reverse quot'range loop
              if quot(i) /= 0
               then sizequot := i; exit;
              end if;
            end loop;
            Clear(n1);
            n1 := Create(quot(0..sizequot));
            r := Create(rest);
          end;
        end if;
      end;
    end if;
  end Div;

-- DESTRUCTOR :

  procedure Clear ( n : in out Natural_Number ) is

    procedure free is
      new unchecked_deallocation(Natural_Number_Rep,Natural_Number);

  begin
    if not Empty(n)
     then free(n);
    end if;
  end Clear;

end Multprec_Natural64_Numbers;
