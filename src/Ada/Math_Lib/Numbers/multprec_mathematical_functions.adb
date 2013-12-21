with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;

package body Multprec_Mathematical_Functions is

-- AUXILIARIES TO IMPLEMENT LN(x) :

  procedure ln_one_plus_x ( x,eps : in Floating_Number;
                            nit : out natural32; max : in natural32;
                            y : out Floating_Number; fail : out boolean ) is

  -- DESCRIPTION :
  --   Uses the standard series expansion to approximate ln(1+x) with
  --   the prescribed precision given by eps.

    sum,acc,quot : Floating_Number;   
    nk : natural32 := 0;
    fk : Floating_Number := Create(natural32(1));
    plus_sign : boolean := false;

  begin
    Copy(x,sum);
    Copy(x,acc);
    for i in 1..max loop
      nk := nk + 1;
      Mul(acc,x);
      Add(fk,1.0);
      quot := acc/fk;
      if plus_sign
       then Add(sum,quot);
       else Sub(sum,quot);
      end if;
      if quot > 0.0
       then fail := (quot > eps);
       else fail := (quot < eps);
      end if;
      Clear(quot);
      exit when not fail;
      plus_sign := not plus_sign;
    end loop;
    Clear(acc);
    nit := nk;
    y := sum;
  end ln_one_plus_x;

  procedure Subtract_Multiple_of_ln2
              ( x : in out Floating_Number; a : in natural32 ) is

  -- DESCRIPTION :
  --   Multiplies x with a*ln(2).

    size : constant natural32 := Size_Fraction(x);
    two : Floating_Number := Create(natural32(2));
    ln2 : Floating_Number;
    flt_a : Floating_Number := Create(a);

  begin
    Set_Size(two,size);
    ln2 := LN(two);
    if a > 1 then
      flt_a := Create(a);
      Mul(ln2,flt_a);
      Clear(flt_a);
    end if;
    Sub(x,ln2);
    Clear(two); Clear(ln2);
  end Subtract_Multiple_of_ln2;

  function Restricted_Natural_Logarithm
             ( x : Floating_Number ) return Floating_Number is

  -- DESCRIPTION :
  --   Returns the natural logarithm of a nonzero positive number &lt;= 1.5.

    res,wrk : Floating_Number;
    size : constant natural32 := Size_Fraction(x);
    nb_digits : constant natural32 := Size_to_Decimal(size);
    dbl_eps : constant double_float := 10.0**(-integer(nb_digits));
    eps : Floating_Number := Create(dbl_eps);
    steps,exp_of_two : natural32;
    max_steps : constant natural32 := 5*nb_digits;
    fail : boolean;

  begin
    exp_of_two := 0;
    if x < 0.5 then
      Copy(x,wrk);
      while wrk < 0.5 loop
        Mul(wrk,2.0);
        exp_of_two := exp_of_two + 1;
      end loop;
      Sub(wrk,1.0);
    else
      wrk := x - 1.0;
    end if;
    Set_Size(wrk,size);
    ln_one_plus_x(wrk,eps,steps,max_steps,res,fail);
    Clear(wrk);
   -- if fail
   --  then put("The desired accuracy of ");
   --       put(dbl_eps,3); put(" has NOT been reached in ");
   --       put(steps,1); put_line(" steps!!");
   --  else put("The desired accuracy of ");
   --       put(dbl_eps,3); put(" has been reached in ");
   --       put(steps,1); put_line(" steps.");
          if exp_of_two > 0
           then Subtract_Multiple_of_ln2(res,exp_of_two);
          end if;
   -- end if;
    Clear(eps);
    return res;
  end Restricted_Natural_Logarithm;

-- EXPONENTIAL AND LOGARITHMIC FUNCTIONS :

  function EXP ( x : Floating_Number ) return Floating_Number is

    res : Floating_Number;
    sx : constant double_float := Round(x);
    dbl_expx : constant double_float := exp(sx);
    delta_y,abs_delta_y : Floating_Number;
    size : constant natural32 := Size_Fraction(x);
    nb_digits : constant natural32 := Size_to_Decimal(size);
    dbl_eps : constant double_float := 10.0**(-integer(nb_digits));
    eps : Floating_Number := Create(dbl_eps);
    steps : natural32 := 0;
    max_steps : constant natural32 := 5*nb_digits;
    fail : boolean;

  begin
    res := Create(dbl_expx);
    Set_Size(res,size);
    for i in 1..max_steps loop
      steps := steps + 1;
      delta_y := ln(res);
      Sub(delta_y,x);
      Mul(delta_y,res);
      Sub(res,delta_y);
      abs_delta_y := AbsVal(delta_y);
      Clear(delta_y);
     -- put("|dy| : "); put(abs_delta_y,3); new_line;
      fail := (abs_delta_y > eps);
      Clear(abs_delta_y);
      exit when not fail;
    end loop;
    Clear(eps);
   -- if not fail
   --  then put("Reached accuracy ");
   --  else put("Failed to reach accuracy ");
   -- end if;
   -- put(eps,3); put(" in ");
   -- put(steps,1); put_line(" steps.");
    return res;
  end EXP;

  function LN ( x : Floating_Number ) return Floating_Number is

    res : Floating_Number;
    y : Floating_Number;

  begin
    if x < 0.0 then
      raise numeric_error;
    elsif Equal(x,1.0) then
      res := Create(natural32(0));
    elsif x > 1.5 then
      y := 1.0/x;
      res := Restricted_Natural_Logarithm(y);
      Min(res); Clear(y);
    else
      res := Restricted_Natural_Logarithm(x);
    end if;
    return res;
  end LN;

  function "**" ( x,y : Floating_Number ) return Floating_Number is

    res : Floating_Number;
    yln_x : Floating_Number := LN(x);

  begin
    Mul(yln_x,y);
    res := EXP(yln_x);
    Clear(yln_x);
    return res;
  end "**";

  function "**" ( x : double_float;
                  y : Floating_Number ) return Floating_Number is

    res : Floating_Number;
    size : constant natural32 := Size_Fraction(y);
    flt_x : Floating_Number := Create(x);

  begin
    Set_Size(flt_x,size);
    res := flt_x**y;
    Clear(flt_x);
    return res;
  end "**";

  function "**" ( x : Floating_Number;
                  y : double_float ) return Floating_Number is

    res : Floating_Number;
    size : constant natural32 := Size_Fraction(x);
    flt_y : Floating_Number := Create(y);

  begin
    Set_Size(flt_y,size);
    res := x**flt_y;
    Clear(flt_y);
    return res;
  end "**";

  function LOG2 ( x : Floating_Number ) return Floating_Number is

    res : Floating_Number := LN(x);
    size : constant natural32 := Size_Fraction(x);
    two : Floating_Number := Create(natural32(2));
    ln_two : Floating_Number;

  begin
    Set_Size(two,size);
    ln_two := LN(two);
    Div(res,ln_two);
    Clear(ln_two);
    return res;
  end LOG2;

  function LOG10 ( x : Floating_Number ) return Floating_Number is

    res : Floating_Number := LN(x);
    size : constant natural32 := Size_Fraction(x);
    ten : Floating_Number := Create(natural32(10));
    ln_ten : Floating_Number;

  begin
    Set_Size(ten,size);
    ln_ten := LN(ten);
    Div(res,ln_ten);
    Clear(ln_ten);
    return res;
  end LOG10;

  function SQRT ( x : Floating_Number ) return Floating_Number is

    res : Floating_Number;
    sx : double_float := Round(x);
    sizex : constant natural32 := Size_Fraction(x);
    eps : constant double_float
        := 10.0**(-integer(Decimal_Places_Fraction(x)));
    nstep : constant natural := integer(sqrt(double_float(sizex)));
    sizeres : natural32;

    procedure Iterate ( a : in Floating_Number ) is

    -- DESCRIPTION :
    --   Applies Newton's method to the equation x^2 - a = 0, which yields
    --   the iteration res := res - (res^2 - a)/(2*res).   The value of res
    --   on entry is the square root of x in double precision.
    --   The iteration stops when we can be sure of all decimal places,
    --   although we also put a bound on the number of iterations,
    --   because infinity cycling occasionally happened otherwise.

      oneres,twores,absinc : Floating_Number;
      stop : boolean := false;

    begin
      for i in 1..nstep loop
        Copy(res,oneres);             -- back up value of current iterate
        twores := 2.0*oneres;         -- twores = 2*res
        Mul(res,oneres);              -- res = res^2
        Sub(res,a);                   -- res = res^2 - a, which is the residual
        Div(res,twores);              -- res = (res^2 - a)/(2*res), increment
        Clear(twores);
        absinc := AbsVal(res);        -- save the increment for stop criterium
        stop := (absinc < eps);
        Clear(absinc);
        twores := oneres - res;       -- perform the update
        Copy(twores,res);             -- clean up for next iteration
        Clear(twores); Clear(oneres);
        exit when stop;
      end loop;
    end Iterate;

  begin
    sx := SQRT(sx);
    res := Create(sx); 
    sizeres := Size_Fraction(res);
    if (sx /= 0.0) and (sizeres < sizex) then
      Expand(res,sizex-sizeres);
      if sx >= 1.0 then
        Iterate(x);
      else
        declare
          invx : Floating_Number := 1.0/x;
          invres : Floating_Number := 1.0/res;
        begin
          Copy(invres,res);  Clear(invres);
          Iterate(invx);
          invres := 1.0/res;
          Copy(invres,res);  Clear(invres);
          Clear(invx);
        end;
      end if;
    end if;
    return res;
  end SQRT;

-- TRIGONOMETRIC FUNCTIONS :

  function SIN ( x : Floating_Number ) return Floating_Number is

    res : Floating_Number;
    sx : double_float := Round(x);

  begin
    sx := SIN(sx);
    res := Create(sx);
    return res;
  end SIN;

  function COS ( x : Floating_Number ) return Floating_Number is

    res : Floating_Number;
    sx : double_float := Round(x);

  begin
    sx := COS(sx);
    res := Create(sx);
    return res;
  end COS;

  function TAN ( x : Floating_Number ) return Floating_Number is

    res : Floating_Number;
    sx : double_float := Round(x);

  begin
    sx := TAN(sx);
    res := Create(sx);
    return res;
  end TAN;

  function ARCSIN ( x : Floating_Number ) return Floating_Number is

    res : Floating_Number;
    sx : double_float := Round(x);

  begin
    sx := ARCSIN(sx);
    res := Create(sx);
    return res;
  end ARCSIN;

  function ARCCOS ( x : Floating_Number ) return Floating_Number is

    res : Floating_Number;
    sx : double_float := Round(x);

  begin
    sx := ARCCOS(sx);
    res := Create(sx);
    return res;
  end ARCCOS;

  function ARCTAN ( x : Floating_Number ) return Floating_Number is

    res : Floating_Number;
    sx : double_float := Round(x);

  begin
    sx := ARCTAN(sx);
    res := Create(sx);
    return res;
  end ARCTAN;

  function Radius ( x,y : Floating_Number ) return Floating_Number is

    res : Floating_Number;
    sx : double_float := Round(x);
    sy : constant double_float := Round(y);

  begin
    sx := Radius(sx,sy);
    res := Create(sx);
    return res;
  end Radius;

  function Angle  ( x,y : Floating_Number ) return Floating_Number is

    res : Floating_Number;
    sx : double_float := Round(x);
    sy : constant double_float := Round(y);

  begin
    sx := Angle(sx,sy);
    res := Create(sx);
    return res;
  end Angle;

end Multprec_Mathematical_Functions;
