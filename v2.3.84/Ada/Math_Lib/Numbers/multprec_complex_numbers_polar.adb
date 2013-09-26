with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_Polar;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Multprec_Mathematical_Functions;    use Multprec_Mathematical_Functions;
with Multprec_Complex_Number_Tools;      use Multprec_Complex_Number_Tools;

package body Multprec_Complex_Numbers_Polar is

  function Radius ( c : Complex_Number ) return Floating_Number is

    res,acc1,acc2 : Floating_Number;

  begin
    acc1 := REAL_PART(c)*REAL_PART(c);
    acc2 := IMAG_PART(c)*IMAG_PART(c);
    Add(acc1,acc2); Clear(acc2);
    res := SQRT(acc1);
    Clear(acc1);
    return res;
  end Radius;

  function Angle ( c : Complex_Number ) return Floating_Number is
  begin
    return Multprec_Mathematical_Functions.Angle(IMAG_PART(c),REAL_PART(c));
  end Angle;

  function Root ( c : Complex_Number; n,i : natural32 ) return Complex_Number is

    res : Complex_Number;
    stc,stres : Standard_Complex_Numbers.Complex_Number;
    rpc : Floating_Number;
    sizeres : natural32;

    procedure Iterate ( x : in out Complex_Number ) is
  
    -- DESCRIPTION :
    --   Performs x := x - (x^n-1)/(n*x^(n-1)) until sufficiently accurate.

      deci : constant natural32 := Size_to_Decimal(sizeres);
      eps : constant double_float := 10.0**(-integer(deci));
      maxstep : constant natural := integer(sqrt(double_float(sizeres)));
      prod,quot : Complex_Number;
      absacc : Floating_Number;
      mpn : Complex_Number := Create(n);
      stop : boolean;

    begin
      for k in 1..maxstep loop
        Copy(x,prod);
        for l in 1..(n-2) loop
          Mul(prod,x);
        end loop;                        -- prod is x^(n-1)
        quot := mpn*prod;                -- quot is n*x^(n-1)
        Mul(prod,x);                     -- prod is x^n
        Sub(prod,c);                     -- prod is x^n - c
        Div(prod,quot);                  -- prod contains update
        Sub(x,prod);                     -- update x
        absacc := AbsVal(prod);          -- convergence test
        stop := (absacc < eps);
        Clear(prod); Clear(quot);        -- clean up for next iteration
        Clear(absacc); 
        exit when stop;
      end loop;
      Clear(mpn);
    end Iterate;

  begin
    if n <= 1 then
      Copy(c,res);
    else
      stc := Round(c);
      stres := Standard_Complex_Numbers_Polar.Root(stc,n,i);
      rpc := REAL_PART(c);
      sizeres := Size_Fraction(rpc);
      res := Create(stres);
      Set_Size(res,sizeres);
      Iterate(res);
    end if;
    return res;
  end Root;

end Multprec_Complex_Numbers_Polar;
