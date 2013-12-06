with unchecked_deallocation;

package body Multprec_Integer_Numbers is

-- NOTES ON THE CHOICE OF REPRESENTATION AND IMPLEMENTATION :
--   0) See also the notes in the body of Multprec_Natural_Numbers.
--      This package inherits the operations on natural numbers, with
--      additionally the tests on signs.
--      Integer numbers are in fact signed natural numbers.
--   1) The construction of tagged records was judged not appropriate to
--      extend the natural numbers, as this construction only applies to
--      records, it would have changed the privacy of the implementation.

-- DATA STRUCTURE :

  type Integer_Number_Rep is record
    plus : boolean;
    numb : Natural_Number;
  end record;

  procedure free is
      new unchecked_deallocation(Integer_Number_Rep,Integer_Number);

-- CREATORS :

  function Natural_Create ( n : natural32 ) return Integer_Number is

    res : Integer_Number;
    res_rep : Integer_Number_Rep;

  begin
    res_rep.plus := true;
    res_rep.numb := Create(n);
    res := new Integer_Number_Rep'(res_rep);
    return res;
  end Natural_Create;

  function Create ( n : Array_of_Naturals ) return Integer_Number is

    res : Integer_Number;
    res_rep : Integer_Number_Rep;

  begin
    res_rep.plus := true;
    res_rep.numb := Create(n);
    res := new Integer_Number_Rep'(res_rep);
    return res;
  end Create;

  function Shallow_Create ( n : Natural_Number ) return Integer_Number is

    res : Integer_Number;
    res_rep : Integer_Number_Rep;

  begin
    res_rep.plus := true;
    res_rep.numb := n;
    res := new Integer_Number_Rep'(res_rep);
    return res;
  end Shallow_Create;

  function Create ( n : Natural_Number ) return Integer_Number is

    res : Integer_Number;
    res_rep : Integer_Number_Rep;

  begin
    res_rep.plus := true;
    res_rep.numb := +n;       --Copy(n,res_rep.numb);
    res := new Integer_Number_Rep'(res_rep);
    return res;
  end Create;

  function Create ( i : integer ) return Integer_Number is
  begin
    return Create(integer32(i));
  end Create;

  function Create ( i : integer32 ) return Integer_Number is

    res : Integer_Number;
    n : natural32;

  begin
    if i >= 0 then
      n := natural32(i);
      res := Natural_Create(n);
      res.plus := true;
    else
      n := natural32(-i);
      res := Natural_Create(n);
      res.plus := false;
    end if;
    return res;
  end Create;

  function Convert ( n : Natural_Number ) return Integer_Number is

    res : Integer_Number;
    res_rep : Integer_Number_Rep;

  begin
    res_rep.numb := n;
    res_rep.plus := true;
    res := new Integer_Number_Rep'(res_rep);
    return res;
  end Convert;

  function Create ( i : Integer_Number ) return integer32 is

    res : integer32;
    nres : natural32;

  begin
    if (Empty(i) or else Empty(i.numb)) then
      res := 0;
    else
      nres := Create(i.numb);
      if i.plus
       then res := integer32(nres);
       else res := -integer32(nres);
      end if;
    end if;
    return res;
  end Create;

-- SELECTORS :

  function Empty ( i : Integer_Number ) return boolean is
  begin
    return (i=null);
  end Empty;

  function Size ( i : Integer_Number ) return natural32 is
  begin
    if Empty(i)
     then return 0;
     else return Size(i.numb);
    end if;
  end Size;

  function Coefficient ( i : Integer_Number; k : natural32 ) return natural32 is
  begin
    if (Empty(i) or else (k > Size(i)))
     then return 0;
     else return Coefficient(i.numb,k);
    end if;
  end Coefficient;

  function Coefficients ( i : Integer_Number ) return Array_of_Naturals is

    nullres : constant Array_of_Naturals(0..0) := (0..0 => 0);

  begin
    if not Empty(i)
     then return Coefficients(i.numb);
     else return nullres;
    end if;
  end Coefficients;

  function Decimal_Places ( i : Integer_Number ) return natural32 is
  begin
    if Empty(i)
     then return 0;
     else return Decimal_Places(i.numb);
    end if;
  end Decimal_Places;

  function Positive ( i : Integer_Number ) return boolean is
  begin
    if Empty(i) then
      return false;
    elsif Empty(i.numb) then
      return false;
    --  elsif Equal(i.numb,0)   -- whatever sign you wish to give to 0
    --      then return false;  -- convenient to work with for input
    else
      return i.plus;
    end if;
  end Positive;

  function Negative ( i : Integer_Number ) return boolean is
  begin
    if Empty(i) then
      return false;
    elsif Empty(i.numb) then
      return false;
    --  elsif Equal(i.numb,0)    -- for input of floating-point numbers
    --      then return false;   -- convenient for reading -0.01
    else
      return not i.plus;
    end if;
  end Negative;

  function Sign ( i : Integer_Number ) return integer32 is
  begin
    if Empty(i) or Equal(i,0) then
      return 0;
    elsif Positive(i) then
      return +1;
    else
      return -1;
    end if;
  end Sign;

  function Unsigned ( i : Integer_Number ) return Natural_Number is

    res : Natural_Number;

  begin
    if not Empty(i)
     then res := i.numb;
    end if;
    return res;
  end Unsigned;

-- COMPARISON AND COPYING :

  function Equal ( i1 : Integer_Number; i2 : integer32 ) return boolean is
  begin
    if Empty(i1) then
      return (i2 = 0);
    elsif ((i1.plus and i2 < 0) or else (not i1.plus and i2 > 0)) then
      return false;
    elsif Empty(i1.numb) then
      if i2 = 0
       then return true;
       else return false;
      end if; 
    elsif i2 >= 0 then
      return Equal(i1.numb,natural32(i2));
    else
      return Equal(i1.numb,natural32(-i2));
    end if;
  end Equal;

  function Equal ( i1,i2 : Integer_Number ) return boolean is
  begin
    if Empty(i1) then
      return Equal(i2,0);
    elsif Empty(i2) then
      return Equal(i1,0);
    else if (Positive(i1) and Negative(i2))
            or else (Negative(i1) and Positive(i2))
          then return false;
          else return Equal(i1.numb,i2.numb);
         end if;
    end if;
  end Equal;

  function "<" ( i1 : Integer_Number; i2 : integer32 ) return boolean is
  begin
    if Empty(i1) then
      return (i2 > 0);
    elsif Positive(i1) then
      if i2 <= 0
       then return false;
       else return (i1.numb < natural32(i2));
      end if;
    elsif Negative(i1) then
      if i2 >= 0
       then return true;
       else return (i1.numb > natural32(-i2));
      end if;
    else
      return (i2 > 0);
    end if;
  end "<";

  function "<" ( i1 : integer32; i2 : Integer_Number ) return boolean is
  begin
    if Empty(i2) then
      return (i1 < 0);
    else
      if Positive(i2) then
        if i1 <= 0
         then return true;
         else return (natural32(i1) < i2.numb);
        end if;
      elsif Negative(i2) then
        if i1 >= 0
         then return false;
         else return (natural32(-i1) > i2.numb);
        end if;
      else 
        return (i1 < 0);
      end if;
    end if;
  end "<";

  function "<" ( i1,i2 : Integer_Number ) return boolean is
  begin
    if Empty(i1) then
      return Positive(i2);
    elsif Empty(i2) then
      return Negative(i1);
    elsif Positive(i1) then
      if Negative(i2)
       then return false;
       else return (i1.numb < i2.numb);
      end if;
    elsif Negative(i1) then
      if Positive(i2)
       then return true;
       else return (i1.numb > i2.numb);
      end if;
    else
      return Positive(i2);
    end if;
  end "<";

  function ">" ( i1 : Integer_Number; i2 : integer32 ) return boolean is
  begin
    if Empty(i1) then
      return (i2 < 0);
    elsif Negative(i1) then
      if i2 >= 0
       then return false;
       else return (i1.numb < natural32(-i2));
      end if;
    elsif Positive(i1) then
      if i2 <= 0
       then return true;
       else return (i1.numb > natural32(i2));
      end if;
    else
      return (i2 < 0);
    end if;
  end ">";

  function ">" ( i1 : integer32; i2 : Integer_Number ) return boolean is
  begin
    if Empty(i2) then
      return (i1 > 0);
    elsif Positive(i2) then
      if i1 <= 0
       then return false;
       else return (natural32(i1) > i2.numb);
      end if;
    elsif Negative(i2) then
      if i1 >= 0
       then return true;
       else return (natural32(-i1) < i2.numb);
      end if;
    else
      return (i1 > 0);
    end if;
  end ">";

  function ">" ( i1,i2 : Integer_Number ) return boolean is
  begin
    if Empty(i1) then
      return Negative(i2);
    elsif Empty(i2) then
      return Positive(i1);
    elsif Positive(i1) then
      if Negative(i2)
       then return true;
       else return (i1.numb > i2.numb);
      end if;
    elsif Negative(i1) then
      if Positive(i2)
       then return false;
       else return (i1.numb < i2.numb);
      end if;
    else
      return Negative(i2);
    end if;
  end ">";

  procedure Copy ( i1 : in integer32; i2 : in out Integer_Number ) is
  begin
    Clear(i2);
    i2 := Create(i1);
  end Copy;

  procedure Copy ( i1 : in Integer_Number; i2 : in out Integer_Number ) is
  begin
    Clear(i2);
    if not Empty(i1) then
      declare
        i2rep : Integer_Number_Rep;
      begin
        i2rep.plus := i1.plus;
        i2rep.numb := +i1.numb;
        i2 := new Integer_Number_Rep'(i2rep);
      end;
    end if;
  end Copy;

-- SHIFTS :

  procedure Shift_Left ( i : in out Integer_Number; ns : out natural32 ) is
  begin
    if Empty(i)
     then ns := 0;
     else Shift_Left(i.numb,ns);
    end if;
  end Shift_Left;

  procedure Shift_Right ( i : in out Integer_Number; ns : out natural32 ) is
  begin
    if Empty(i)
     then ns := 0;
     else Shift_Right(i.numb,ns);
    end if;
  end Shift_Right;

-- ARITHMETIC OPERATIONS as functions :

  function "+" ( i1 : Integer_Number; i2 : integer32 ) return Integer_Number is

    res : Integer_Number;
    res_rep : Integer_Number_Rep;
    n : natural32;

  begin
    if (Empty(i1) or else Empty(i1.numb)) then
      res := Create(i2);
    else
      if i1.plus then
        if i2 >= 0 then
          n := natural32(i2);
          res_rep.plus := true;
          res_rep.numb := i1.numb + n;
          res := new Integer_Number_Rep'(res_rep);
        else
          n := natural32(-i2);
          if not Equal(i1.numb,n) then
            if i1.numb > n then
              res_rep.plus := true;
              res_rep.numb := i1.numb - n;
            else
              res_rep.plus := false;
              res_rep.numb := n - i1.numb;
            end if;
            res := new Integer_Number_Rep'(res_rep);
          end if;
        end if;
      else
        if i2 <= 0 then
          n := natural32(-i2);
          res_rep.plus := false;
          res_rep.numb := i1.numb + n;
          res := new Integer_Number_Rep'(res_rep);
        else
          n := natural32(i2);
          if not Equal(i1.numb,n) then
            if i1.numb < n then
              res_rep.plus := true;
              res_rep.numb := n - i1.numb;
            else
              res_rep.plus := false;
              res_rep.numb := i1.numb - n;
            end if;
            res := new Integer_Number_Rep'(res_rep);
          end if;
        end if;
      end if;
    end if;
    return res;
  end "+";

  function "+" ( i1 : integer32; i2 : Integer_Number ) return Integer_Number is
  begin
    return (i2+i1);
  end "+";

  function "+" ( i1,i2 : Integer_Number ) return Integer_Number is

    res : Integer_Number;
    res_rep : Integer_Number_Rep;

  begin
    if (Empty(i1) or else Empty(i1.numb)) then
      Copy(i2,res);
    else
      if (Empty(i2) or else Empty(i2.numb)) then
        Copy(i1,res);
      else
        if i1.plus then
          if i2.plus then
            res_rep.plus := true;
            res_rep.numb := i1.numb + i2.numb;
            res := new Integer_Number_Rep'(res_rep);
          else 
            if not Equal(i1.numb,i2.numb) then
              if i1.numb > i2.numb then
                res_rep.plus := true;
                res_rep.numb := i1.numb - i2.numb;
              else
                res_rep.plus := false;
                res_rep.numb := i2.numb - i1.numb;
              end if;
              res := new Integer_Number_Rep'(res_rep);
            end if;
          end if;
        else
          if not i2.plus then
            res_rep.plus := false;
            res_rep.numb := i1.numb + i2.numb;
            res := new Integer_Number_Rep'(res_rep);
          else
            if not Equal(i1.numb,i2.numb) then
              if i1.numb < i2.numb then
                res_rep.plus := true;
                res_rep.numb := i2.numb - i1.numb;
              else
                res_rep.plus := false;
                res_rep.numb := i1.numb - i2.numb;
              end if;
              res := new Integer_Number_Rep'(res_rep);
            end if;
          end if;
        end if;
      end if;
    end if;
    return res;
  end "+";

  function "+" ( i : Integer_Number ) return Integer_Number is

    res : Integer_Number;

  begin
    Copy(i,res);
    return res;
  end "+";

  function "-" ( i : Integer_Number ) return Integer_Number is

    res : Integer_Number;
    res_rep : Integer_Number_Rep;

  begin
   -- if not Empty(i) then
    if not Equal(i,0) then
      res_rep.plus := not i.plus;
      res_rep.numb := +i.numb;    -- Copy(i.numb,res_rep.numb);
      res := new Integer_Number_Rep'(res_rep);
    end if;
    return res;
  end "-";

  function "-" ( i1 : Integer_Number; i2 : integer32 ) return Integer_Number is

    mini2 : constant integer32 := -i2;

  begin
    return (i1+mini2);
  end "-";

  function "-" ( i1 : integer32; i2 : Integer_Number ) return Integer_Number is

    res : Integer_Number := i2 - i1;

  begin
    Min(res);
    return res;
  end "-";

  function "-" ( i1,i2 : Integer_Number ) return Integer_Number is

    res,mini2 : Integer_Number;
    mini2rep : Integer_Number_Rep;

  begin
    if (Empty(i2) or else Empty(i2.numb)) then
      Copy(i1,res);
    else
      mini2rep.numb := i2.numb;
      mini2rep.plus := not i2.plus;
      mini2 := new Integer_Number_Rep'(mini2rep);
      res := i1 + mini2;
      free(mini2);
    end if;
    return res;
  end "-";

  function "*" ( i1 : Integer_Number; i2 : integer32 ) return Integer_Number is

    res : Integer_Number;
    res_rep : Integer_Number_Rep;
    n : natural32;

  begin
    if not ((i2 = 0) or else Empty(i1) or else Empty(i1.numb)) then
      if i2 > 0 then
        n := natural32(i2);
        res_rep.plus := i1.plus;
      else
        n := natural32(-i2);
        res_rep.plus := not i1.plus;
      end if;
      res_rep.numb := i1.numb*n;
      res := new Integer_Number_Rep'(res_rep);
    end if;
    return res;
  end "*";

  function "*" ( i1 : integer32; i2 : Integer_Number ) return Integer_Number is
  begin
    return (i2*i1);
  end "*";

  function "*" ( i1,i2 : Integer_Number ) return Integer_Number is

    res : Integer_Number;
    res_rep : Integer_Number_Rep;

  begin
    if (not (Empty(i1) or else Empty(i1.numb)))
       and then (not (Empty(i2) or else Empty(i2.numb)))
     then res_rep.numb := i1.numb*i2.numb;
          res_rep.plus := i1.plus;
          if not i2.plus
           then res_rep.plus := not res_rep.plus;
          end if;
          res := new Integer_Number_Rep'(res_rep);
    end if;
    return res;
  end "*";

  function "**" ( i : Integer_Number; n : natural32 ) return Integer_Number is

    res : Integer_Number;
    res_rep : Integer_Number_Rep;

  begin
    if n = 0 then
      res := Create(integer32(1));
    else
      if not (Empty(i) or else Empty(i.numb)) then
        res_rep.numb := i.numb**n;
        res_rep.plus := i.plus;
        if ((not i.plus) and then (n mod 2 = 1))
         then res_rep.plus := not res_rep.plus;
        end if;
        res := new Integer_Number_Rep'(res_rep);
      end if;
    end if;
    return res;
  end "**";

  function "**" ( i : integer32; n : Natural_Number ) return Integer_Number is

    res : Integer_Number;
    res_rep : Integer_Number_Rep;
    ni : natural32;

  begin
    if (Empty(n) or else Equal(n,0)) then
      res := Create(integer32(1));
    else
      if i /= 0 then
        if i > 0 then
          ni := natural32(i);
          res_rep.plus := true;
        else
          ni := natural32(-i);
          res_rep.plus := false;
        end if;
        res_rep.numb := ni**n;
        if (i < 0 and then (Rmd(n,2) = 0))
         then res_rep.plus := not res_rep.plus;
        end if;
        res := new Integer_Number_Rep'(res_rep);
      end if;
    end if; 
    return res;
  end "**";

  function "**" ( i : Integer_Number; n : Natural_Number )
                return Integer_Number is

    res : Integer_Number;
    res_rep : Integer_Number_Rep;

  begin
    if (Empty(n) or else Equal(n,0)) then
      res := Create(integer32(1));
    else
      if not (Empty(i) or else Empty(i.numb)) then
        res_rep.numb := i.numb**n;
        res_rep.plus := i.plus;
        if ((not i.plus) and then (Rmd(n,2) = 0))
         then res_rep.plus := not res_rep.plus;
        end if;
        res := new Integer_Number_Rep'(res_rep);
      end if;
    end if;
    return res;
  end "**";

  function "/" ( i1 : Integer_Number; i2 : integer32 ) return Integer_Number is

    res : Integer_Number;
    i2n : natural32;
    res_rep : Integer_Number_Rep;

  begin
    if i2 /= 0 then
      if not (Empty(i1) or else Empty(i1.numb)) then
        if i2 > 0
         then i2n := natural32(i2);
         else i2n := natural32(-i2);
        end if;
        res_rep.numb := i1.numb/i2n;
        if (i1.plus and (i2 > 0)) or ((not i1.plus) and (i2 < 0))
         then res_rep.plus := true;
         else res_rep.plus := false;
        end if;
        res := new Integer_Number_Rep'(res_rep);
      end if;
    else
      raise NUMERIC_ERROR;
    end if;
    return res;
  end "/";

  function "/" ( i1 : integer32; i2 : Integer_Number ) return integer32 is

    res : integer32;
    i1n,nres : natural32;

  begin
    if (Empty(i2) or else Empty(i2.numb)) then
      raise NUMERIC_ERROR;
    else
      if i1 > 0
       then i1n := natural32(i1);
       else i1n := natural32(-i1);
      end if;
      nres := i1n/i2.numb;
      if ((i1 > 0) and i2.plus) or ((i1 < 0) and (not i2.plus))
       then res := integer32(nres);
       else res := -integer32(nres);
      end if;
    end if;
    return res;
  end "/";

  function "/" ( i1,i2 : Integer_Number ) return Integer_Number is

    res : Integer_Number;
    res_rep : Integer_Number_Rep;

  begin
    if not (Empty(i1) or else Empty(i1.numb)) then
      if (Empty(i2) or else Empty(i2.numb)) then
        raise NUMERIC_ERROR;
      else
        res_rep.numb := i1.numb/i2.numb;
        if (i1.plus and i2.plus) or ((not i1.plus) and (not i2.plus))
         then res_rep.plus := true;
         else res_rep.plus := false;
        end if;
        res := new Integer_Number_Rep'(res_rep);
      end if;
    end if;
    return res;
  end "/";

  function Rmd ( i1 : Integer_Number; i2 : integer32 ) return integer32 is

    res : integer32;
    i2n,nres : natural32;

  begin
    if i2 /= 0 then
      if (Empty(i1) or else Empty(i1.numb)) then
        res := 0;
      else 
        if i2 > 0
         then i2n := natural32(i2);
         else i2n := natural32(-i2);
        end if;
        nres := Rmd(i1.numb,i2n);
        if i1.plus
         then res := integer32(nres);
         else res := -integer32(nres);
        end if;
      end if;
    else
      raise NUMERIC_ERROR;
    end if;
    return res;
  end Rmd;

  function Rmd ( i1 : integer32; i2 : Integer_Number ) return integer32 is

    res : integer32;
    i1n,nres : natural32;

  begin
    if i1 = 0 then
      res := 0;
    else
      if (Empty(i2) or else Empty(i2.numb)) then
        raise NUMERIC_ERROR;
      else
        if i1 > 0
         then i1n := natural32(i1);
         else i1n := natural32(-i1);
        end if; 
        nres := Rmd(i1n,i2.numb);
        if i1 > 0
         then res := integer32(nres);
         else res := -integer32(nres);
        end if;
      end if;
    end if;
    return res;
  end Rmd;

  function Rmd ( i1,i2 : Integer_Number ) return Integer_Number is

    res : Integer_Number;
    res_rep : Integer_Number_Rep;

  begin
    if not (Empty(i1) or else Empty(i1.numb)) then
      if (Empty(i2) or else Empty(i2.numb)) then
        raise NUMERIC_ERROR;
      else
        res_rep.numb := Rmd(i1.numb,i2.numb);
        res_rep.plus := i1.plus;
        res := new Integer_Number_Rep'(res_rep);
      end if;
    end if;
    return res;
  end Rmd;

-- SPECIAL ARITHMETIC OPERATION : multiply with radix = shift left

  procedure Mul_Radix ( i : in out Integer_Number; k : in natural32 ) is
  begin
    Mul_Radix(i.numb,k);
  end Mul_Radix;

-- ARITHMETIC OPERATIONS as procedures for memory management :

  procedure Add ( i1 : in out Integer_Number; i2 : in integer32 ) is

    n : natural32;
    nn : Natural_Number;

  begin
    if (Empty(i1) or else Empty(i1.numb)) then
      i1 := Create(i2);
    else
      if i1.plus then
        if i2 >= 0 then
          n := natural32(i2);
          Add(i1.numb,n);
        else
          n := natural32(-i2);
          if not Equal(i1.numb,n) then
            if i1.numb > n then
              Sub(i1.numb,n);
            else
              i1.plus := false;
              nn := Create(n);
              Sub(nn,i1.numb);
              Clear(i1.numb); i1.numb := nn;
            end if;
          else
            Clear(i1);
          end if;
        end if;
      else
        if i2 <= 0 then
          n := natural32(-i2);
          Add(i1.numb,n);
        else
          n := natural32(i2);
          if not Equal(i1.numb,n) then
            if i1.numb < n then
              i1.plus := true;
              nn := Create(n);
              Sub(nn,i1.numb);
              Clear(i1.numb); i1.numb := nn;
            else
              Sub(i1.numb,n);
            end if;
          else
            Clear(i1);
          end if;
        end if;
      end if;
    end if;
  end Add;

  procedure Add ( i1 : in out Integer_Number; i2 : in Integer_Number ) is

    nn : Natural_Number;

  begin
    if (Empty(i1) or else Empty(i1.numb)) then
      Copy(i1 => i2,i2 => i1);
    else
      if not (Empty(i2) or else Empty(i2.numb)) then
        if i1.plus then
          if i2.plus then
            Add(i1.numb,i2.numb);
          else
            if not Equal(i1.numb,i2.numb) then
              if i1.numb > i2.numb then
                Sub(i1.numb,i2.numb);
              else
                Copy(i2.numb,nn);
                Sub(nn,i1.numb);
                Clear(i1.numb);
                i1.plus := false;
                i1.numb := nn;
              end if;
            else
              Clear(i1);
            end if;
          end if;
        else
          if not i2.plus then
            Add(i1.numb,i2.numb);
          else
            if not Equal(i1.numb,i2.numb) then
              if i1.numb < i2.numb then
                Copy(i2.numb,nn);
                Sub(nn,i1.numb);
                Clear(i1.numb);
                i1.plus := true;
                i1.numb := nn;
              else
                Sub(i1.numb,i2.numb);
              end if;
            else
              Clear(i1);
            end if;
          end if;
        end if;
      end if;
    end if;
  end Add;

  procedure Min ( i : in out Integer_Number ) is
  begin
   -- if not Empty(i)
    if not Equal(i,0)
     then i.plus := not i.plus;
    end if;
  end Min;

  procedure Set_Min ( i : in out Integer_Number ) is
  begin
    if not Empty(i)
     then i.plus := not i.plus;
    end if;
  end Set_Min;

  procedure Sub ( i1 : in out Integer_Number; i2 : in integer32 ) is
  begin
    Add(i1,-i2);
  end Sub;

  procedure Sub ( i1 : in out Integer_Number; i2 : in Integer_Number ) is

    mini2 : Integer_Number;
    mini2rep : Integer_Number_Rep;

  begin
    if not (Empty(i2) or else Empty(i2.numb)) then
      mini2rep.numb := i2.numb;
      mini2rep.plus := not i2.plus;
      mini2 := new Integer_Number_Rep'(mini2rep);
      Add(i1,mini2);
      free(mini2);
    end if;
  end Sub;

  procedure Mul ( i1 : in out Integer_Number; i2 : in integer32 ) is

    n : natural32;

  begin
    if not (Empty(i1) or else Empty(i1.numb)) then
      if i2 = 0 then
        Clear(i1);
      else
        if i2 > 0 then
          n := natural32(i2);
        else
          n := natural32(-i2);
          i1.plus := not i1.plus;
        end if;
        Mul(i1.numb,n);
      end if;
    end if;
  end Mul;

  procedure Mul ( i1 : in out Integer_Number; i2 : in Integer_Number ) is
  begin
    if (not (Empty(i1) or else Empty(i1.numb))) then
      if (Empty(i2) or else Empty(i2.numb)) then
        Clear(i1);
      else
        Mul(i1.numb,i2.numb);
        if not i2.plus
         then i1.plus := not i1.plus;
        end if;
      end if;
    end if;
  end Mul;

  procedure Rmd ( i1 : in out Integer_Number; i2 : in integer32 ) is

    res : constant Integer_Number := Create(Rmd(i1,i2));

  begin
    Clear(i1); i1 := res;
  end Rmd;

  procedure Rmd ( i1 : in out Integer_Number; i2 : in Integer_Number ) is

    res : constant Integer_Number := Rmd(i1,i2);

  begin
    Clear(i1); i1 := res;
  end Rmd;

  procedure Div ( i1 : in out Integer_Number; i2 : in integer32 ) is

    r : integer32;

  begin
    Div(i1,i2,r);
  end Div;

  procedure Div ( i1 : in out Integer_Number; i2 : in Integer_Number ) is

    r : Integer_Number;

  begin
    Div(i1,i2,r);
    Clear(r);
  end Div;

  procedure Div ( i1 : in Integer_Number; i2 : in integer32;
                  q : out Integer_Number; r : out integer32 ) is

    qrep : Integer_Number_Rep;
    i2n,rn : natural32;

  begin
    if i2 /= 0 then
      if not (Empty(i1) or else Empty(i1.numb)) then
        if i2 > 0
         then i2n := natural32(i2);
         else i2n := natural32(-i2);
        end if; 
        Div(i1.numb,i2n,qrep.numb,rn);
        if (i1.plus and (i2 > 0)) or ((not i1.plus) and (i2 < 0))
         then qrep.plus := true;
         else qrep.plus := false;
        end if;
        q := new Integer_Number_Rep'(qrep);
        if i1.plus
         then r := integer32(rn);
         else r := -integer32(rn);
        end if;
      end if;
    else
      raise NUMERIC_ERROR;
    end if;
  end Div;

  procedure Div ( i1 : in out Integer_Number; i2 : in integer32;
                  r : out integer32 ) is

    i2n,rn : natural32;
 
  begin
    if i2 /= 0 then
      if not (Empty(i1) or else Empty(i1.numb)) then
        if i2 > 0
         then i2n := natural32(i2);
         else i2n := natural32(-i2);
        end if;
        Div(i1.numb,i2n,rn);
        if i1.plus
         then r := integer32(rn);
         else r := -integer32(rn);
        end if;
        if (i1.plus and (i2 > 0)) or ((not i1.plus) and (i2 < 0))
         then i1.plus := true;
         else i1.plus := false;
        end if;
      end if;
    else
      raise NUMERIC_ERROR;
    end if;
  end Div;

  procedure Div ( i1,i2 : in Integer_Number; q,r : out Integer_Number ) is

    qrep,rrep : Integer_Number_Rep;

  begin
    if not (Empty(i2) or else Empty(i2.numb)) then
      if not (Empty(i1) or else Empty(i1.numb)) then
        Div(i1.numb,i2.numb,qrep.numb,rrep.numb);
        if (i1.plus and (i2 > 0)) or ((not i1.plus) and (i2 < 0))
         then qrep.plus := true;
         else qrep.plus := false;
        end if;
        q := new Integer_Number_Rep'(qrep);
        rrep.plus := i1.plus;
        r := new Integer_Number_Rep'(rrep);
      end if;
    else
      raise NUMERIC_ERROR;
    end if;
  end Div;

  procedure Div ( i1 : in out Integer_Number; i2 : in Integer_Number;
                  r : out Integer_Number ) is

    rrep : Integer_Number_Rep;

  begin
    if not (Empty(i2) or else Empty(i2.numb)) then
      if not (Empty(i1) or else Empty(i1.numb)) then
        Div(i1.numb,i2.numb,rrep.numb);
        rrep.plus := i1.plus;
        r := new Integer_Number_Rep'(rrep);
        if (i1.plus and (i2 > 0)) or ((not i1.plus) and (i2 < 0))
         then i1.plus := true;
         else i1.plus := false;
        end if;
      end if;
    else
      raise NUMERIC_ERROR;
    end if;
  end Div;

-- DESTRUCTOR :

  procedure Shallow_Clear ( i : in out Integer_Number ) is
  begin
    if not Empty(i)
     then free(i);
    end if;
  end Shallow_Clear;

  procedure Clear ( i : in out Integer_Number ) is
  begin
    if not Empty(i) then
      Clear(i.numb);
      free(i);
      i := null;
    end if;
  end Clear;

end Multprec_Integer_Numbers;
