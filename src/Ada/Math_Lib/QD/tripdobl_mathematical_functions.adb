with text_io;                            use text_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Mathematical_Functions;
with Triple_Double_Constants;            use Triple_Double_Constants;

package body TripDobl_Mathematical_Functions is

  function SQRT ( x : triple_double ) return triple_double is

  --  Perform the following Newton iteration:
  --    y' = y + (1 - x * y^2) * y / 2;
  --  which converges to 1/sqrt(x), starting with the
  --  double precision approximation to 1/sqrt(x).
  --  Since Newton's iteration more or less doubles the
  --  number of correct digits, we only need to perform it twice.

    res,h : triple_double;
    f : double_float;
    one : constant triple_double := create(1.0);

  begin
    if is_zero(x) then
      res := Create(0.0);
    elsif is_negative(x) then  -- original code returns NaN
      res := Create(-1.0);
    else
      f := Standard_Mathematical_Functions.SQRT(hi_part(x));
      res := Create(f);
      res := one/res;
      h := mul_pwr2(x,0.5);
      res := res + ((0.5 - h*(res*res))*res);
      res := res + ((0.5 - h*(res*res))*res);
      res := res*x;
    end if;
    return res;
  end SQRT;

-- AUXILIARY ROUTINES for SIN and COS :

  procedure sincos_taylor
              ( x : in triple_double; sin_x,cos_x : out triple_double ) is

  -- DESCRIPTION :
  --   Assuming |x| <= pi/2048, returns in sin_x and cos_x
  --   a Taylor series approximation respectivley for sin(x) and cos(x).

    f : double_float := to_double(x);
    thresh : constant double_float := 0.5*td_eps*f;
    p,s,t,y : triple_double; 
    i : integer;

  begin
    if is_zero(x) then
      sin_x := create(0.0);
      cos_x := create(1.0);
    else
      y := -sqr(x); s := x; p := x;
      i := 0;
      loop
        p := p*y;
        t := p*i_fac(i);
        s := s + t;
        i := i + 2;
        exit when (i >= n_inv_fact);
        f := to_double(t);
        if f < 0.0
         then f := -f;
        end if;
        exit when (f <= thresh);
      end loop;
    end if;
    sin_x := s;
    cos_x := SQRT(1.0 - sqr(s));
  end sincos_taylor;

  function sin_taylor ( x : triple_double ) return triple_double is

  -- DESCRIPTION :
  --   Assuming |x| <= pi/2048,
  --   returns a Taylor series approximation for sin(x).

    res : triple_double;
    f : double_float := to_double(x);
    thresh : constant double_float := 0.5*td_eps*f;
    p,t,y : triple_double; 
    i : integer;

  begin
    if is_zero(x) then
      res := create(0.0);
    else
      y := -sqr(x); res := x; p := x;
      i := 0;
      loop
        p := p*y;
        t := p*i_fac(i);
        res := res + t;
        i := i + 2;
        exit when (i >= n_inv_fact);
        f := to_double(t);
        if f < 0.0
         then f := -f;
        end if;
        exit when (f <= thresh);
      end loop;
    end if;
    return res;
  end sin_taylor;

  function cos_taylor ( x : triple_double ) return triple_double is

  -- DESCRIPTION :
  --   Assuming |x| <= pi/2048,
  --   returns a Taylor series approximation for cos(x).

    res : triple_double;
    f : double_float := to_double(x);
    thresh : constant double_float := 0.5*td_eps*f;
    p,t,y : triple_double; 
    i : integer;

  begin
    if is_zero(x) then
      res := create(1.0);
    else
      y := -sqr(x);
      res := 1.0 + mul_pwr2(y,0.5);
      p := y;
      i := 1;
      loop
        p := p*y;
        t := p*i_fac(i);
        res := res + t;
        i := i + 2;
        exit when (i >= n_inv_fact);
        f := to_double(t);
        if f < 0.0
         then f := -f;
        end if;
        exit when (f <= thresh);
      end loop;
    end if;
    return res;
  end cos_taylor;

  procedure Reduce_Modulo_2pi
              ( x : in triple_double; t : out triple_double;
                j,k,abs_k : out integer; fail : out boolean ) is

  -- DESCRIPTION :
  --   Reduces x modulo 2*pi, modulo pi/2, and then modulo pi/1024.
  --   If this reduction does not work, then an error message is printed
  --   and fail is true on return.

  -- ON ENTRY :
  --   x        some triple double number.

  -- ON RETURN :
  --   t        what remains of x;
  --   j        result after reduction modulo pi/2;
  --   k        result after reduction modulo pi/16;
  --   abs_k    absolute value of k;
  --   fail     true if reduction fails.

    z,r : triple_double;
    q : double_float;

  begin
    z := nint(x/twopi); -- approximately reduce modulo 2*pi
    r := x - twopi*z;
   -- approximately reduce modulo pi/2 and then modulo pi/1024
    q := double_float'floor(hi_part(r)/hi_part(pi2) + 0.5);
    t := r - pi2*q;
    j := integer(q);
    if j < -2 or j > 2 then
      put_line("qd_sin: cannot reduce modulo pi/2");
      fail := true;
    else
      q := double_float'floor(hi_part(t)/hi_part(pi1024) + 0.5);
      t := t - pi1024*q;
      k := integer(q);
      if k < 0
       then abs_k := -k;
       else abs_k := k;
      end if;
      if abs_k > 256 then       
        put_line("qd_sin: cannot reduce modulo pi/1024");
        fail := true;
      else
        fail := false;
      end if;
    end if;
  end Reduce_Modulo_2pi;

-- TARGET TRIGONOMETRIC FUNCTIONS :

  function SIN ( x : triple_double ) return triple_double is

  -- DESCRIPTION :
  --   To compute sin(x), we choose integers a, b so that
  --     x = s + a * (pi/2) + b * (pi/1024)
  --   and |s| <= pi/2048.  Using a precomputed table of
  --   sin(k pi / 1024) and cos(k pi / 1024), we can compute
  --   sin(x) from sin(s) and cos(s).  This greatly increases the
  --   convergence of the sine Taylor series. 

    res : triple_double;

  begin
    if is_zero(x) then
      res := create(0.0);
    else
      declare
        t : triple_double;
        j,k,abs_k : integer;
        fail : boolean;
      begin
        Reduce_Modulo_2pi(x,t,j,k,abs_k,fail);
        if fail then       
          res := create(-2.0); -- instead of NaN
        else
          if k = 0 then
            case j is
              when 0 => res := sin_taylor(t);
              when 1 => res := cos_taylor(t);
              when -1 => res := -cos_taylor(t);
              when others => res := -sin_taylor(t);
            end case;
          else
            declare
              u,v,t_sin,t_cos : triple_double;
            begin
              u := cos_table(abs_k-1);
              v := sin_table(abs_k-1);
              sincos_taylor(t,t_sin,t_cos);
              if j = 0 then
                if k > 0
                 then res := u*t_sin + v*t_cos;
                 else res := u*t_sin - v*t_cos;
                end if;
              elsif j = 1 then
                if k > 0
                 then res := u*t_cos - v*t_sin;
                 else res := u*t_cos + v*t_sin;
                end if;
              elsif j = -1 then
                if k > 0
                 then res :=  v*t_sin - u*t_cos;
                 else res := -u*t_cos - v*t_sin;
                end if;
              else
                if k > 0
                 then res := -u*t_sin - v*t_cos;
                 else res :=  v*t_cos - u*t_sin;
                end if;
              end if;
            end;
          end if;
        end if;
      end;
    end if;
    return res;
  end SIN;

  function COS ( x : triple_double ) return triple_double is

    res : triple_double;

  begin
    if is_zero(x) then
      res := create(1.0);
    else
      declare
        t : triple_double;
        j,k,abs_k : integer;
        fail : boolean;
      begin
        Reduce_Modulo_2pi(x,t,j,k,abs_k,fail);
        if fail then       
          res := create(-2.0); -- instead of NaN
        else
          if k = 0 then
            case j is
              when 0 => res := cos_taylor(t);
              when 1 => res := -sin_taylor(t);
              when -1 => res := sin_taylor(t);
              when others => res := -cos_taylor(t);
            end case;
          else
            declare
              u,v,t_sin,t_cos : triple_double;
            begin
              u := cos_table(abs_k-1);
              v := sin_table(abs_k-1);
              sincos_taylor(t,t_sin,t_cos);
              if j = 0 then
                if k > 0
                 then res := u*t_cos - v*t_sin;
                 else res := u*t_cos + v*t_sin;
                end if;
              elsif j = 1 then
                if k > 0
                 then res := -u*t_sin - v*t_cos;
                 else res :=  v*t_cos - u*t_sin;
                end if;
              elsif j = -1 then
                if k > 0
                 then res := u*t_sin + v*t_cos;
                 else res := u*t_sin - v*t_cos;
                end if;
              else
                if k > 0
                 then res :=  v*t_sin - u*t_cos;
                 else res := -u*t_cos - v*t_sin;
                end if;
              end if;
            end;
          end if;
        end if;
      end;
    end if;
    return res;
  end COS;

  procedure SINCOS ( x : in triple_double; sin_x,cos_x : out triple_double ) is
  begin
    if is_zero(x) then
      sin_x := Create(0.0);
      cos_x := Create(1.0);
    else
      declare
        t : triple_double;
        j,k,abs_k : integer;
        fail : boolean;
      begin
        Reduce_Modulo_2pi(x,t,j,k,abs_k,fail);
        if fail then
          sin_x := Create(-2.0); -- instead of NaN
          cos_x := Create(-2.0);
        else
          declare
            sin_t,cos_t,u,v : triple_double;
          begin
            sincos_taylor(t,sin_t,cos_t);
            if k = 0 then
              case j is 
                when 0 => sin_x := sin_t; cos_x := cos_t;
                when 1 => sin_x := cos_t; cos_x := -sin_t;
                when -1 => sin_x := -cos_t; cos_x := sin_t;
                when others => sin_x := -sin_t; cos_x := -cos_t;
              end case;
            else
              u := cos_table(abs_k-1);
              v := sin_table(abs_k-1);
              if j = 0 then
                if k > 0 then
                  sin_x := u * sin_t + v * cos_t;
                  cos_x := u * cos_t - v * sin_t;
                else
                  sin_x := u * sin_t - v * cos_t;
                  cos_x := u * cos_t + v * sin_t;
                end if;
              elsif j = 1 then
                if k > 0 then
                  cos_x := - u * sin_t - v * cos_t;
                  sin_x := u * cos_t - v * sin_t;
                else
                  cos_x := v * cos_t - u * sin_t;
                  sin_x := u * cos_t + v * sin_t;
                end if;
              elsif j = -1 then
                if k > 0 then
                  cos_x := u * sin_t + v * cos_t;
                  sin_x :=  v * sin_t - u * cos_t;
                else
                  cos_x := u * sin_t - v * cos_t;
                  sin_x := - u * cos_t - v * sin_t;
                end if;
              else
                if k > 0 then
                  sin_x := - u * sin_t - v * cos_t;
                  cos_x := v * sin_t - u * cos_t;
                else
                  sin_x := v * cos_t - u * sin_t;
                  cos_x := - u * cos_t - v * sin_t;
                end if;
              end if;
            end if;
          end;
        end if;
      end;
    end if;
  end SINCOS;

  function TAN ( x : triple_double ) return triple_double is

    s,c : triple_double;

  begin
    sincos(x,s,c);
    return s/c;
  end TAN;

  function ARCTAN ( x : triple_double ) return triple_double is

    one : constant triple_double := create(1.0);

  begin
    return ARCTAN(x,one);
  end ARCTAN;

  function ARCTAN ( y,x : triple_double ) return triple_double is

  -- Instead of using Taylor series to compute arctan,
  -- we instead use Newton's iteration to solve the equation
  -- sin(z) = y/r, or cos(z) = x/r, where r = sqrt(x^2 + y^2).
  -- The iteration is given by
  --   z' = z + (y - sin(z)) / cos(z)          (for equation 1)
  --   z' = z - (x - cos(z)) / sin(z)          (for equation 2)
  -- Here, x and y are normalized so that x^2 + y^2 = 1.
  -- If |x| > |y|, then first iteration is used since the
  -- denominator is larger.  Otherwise, the second is used.

  begin
    if is_zero(x) then
      if is_zero(y) then
        put_line("qd_arctan: both arguments zero");
        return create(0.0); -- instead of NaN
      else
        if is_positive(y)
         then return pi2;
         else return -pi2;
        end if;
      end if;
    elsif is_zero(y) then
      if is_positive(x)
       then return create(0.0);
       else return pi;
      end if;
    elsif x = y then
      if is_positive(y)
       then return pi4;
       else return -threepi4;
      end if;
    elsif x = -y then
      if is_positive(y)
       then return threepi4;
       else return pi4;
      end if;
    else
      declare
        r : constant triple_double := SQRT(sqr(x) + sqr(y));
        xx : constant triple_double := x/r;
        yy : constant triple_double := y/r;
        f : constant double_float
          := Standard_Mathematical_Functions.Angle(to_double(y),to_double(x));
        z : triple_double := create(f);
        sin_z,cos_z : triple_double;
      begin
        sincos(z,sin_z,cos_z);
        if abs(hi_part(xx)) > abs(hi_part(yy)) then
        -- Use Newton iteration 1.  z' = z + (y - sin(z)) / cos(z) 
          z := z + (yy - sin_z) / cos_z;
          sincos(z, sin_z, cos_z);
          z := z + (yy - sin_z) / cos_z;
          sincos(z, sin_z, cos_z);
          z := z + (yy - sin_z) / cos_z;
        else 
        -- Use Newton iteration 2.  z' = z - (x - cos(z)) / sin(z) 
          z := z - (xx - cos_z) / sin_z;
          sincos(z, sin_z, cos_z);
          z := z - (xx - cos_z) / sin_z;
          sincos(z, sin_z, cos_z);
          z := z - (xx - cos_z) / sin_z;
        end if;
        return z;
      end;
    end if;
  end ARCTAN;

  function ARCSIN ( x : triple_double ) return triple_double is
  begin
    if x < -1.0 or x > 1.0 then
      put_line("qd_arcsin: argument out of domain");
      return x; -- instead of NaN
    elsif is_one(x) then
      if is_positive(x)
       then return pi2;
       else return -pi2;
      end if;
    else
      declare
        one : constant triple_double := create(1.0);
        s1x : constant triple_double := SQRT(one - sqr(x));
      begin
        return ARCTAN(x,s1x);
      end;
    end if;
  end ARCSIN;

  function ARCCOS ( x : triple_double ) return triple_double is
  begin
    if x < -1.0 or x > 1.0 then
      put_line("qd_arccos: argument out of domain");
      return x; -- instead of NaN
    elsif is_one(x) then
      if is_positive(x)
       then return create(0.0);
       else return pi;
      end if;
    else
      declare
        one : constant triple_double := create(1.0);
        s1x : constant triple_double := SQRT(one - sqr(x));
      begin
        return ARCTAN(s1x,x);
      end;
    end if;
  end ARCCOS;

  function Radius ( x,y : triple_double ) return triple_double is
  begin
    return SQRT(x*x + y*y);
  end Radius;

  function Angle ( x,y : triple_double ) return triple_double is

    first : constant triple_double := x;
    second : constant triple_double := y;

  begin
    return ARCTAN(first,second);
  end Angle;

end TripDobl_Mathematical_Functions;
