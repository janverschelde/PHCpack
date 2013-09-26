with text_io;                            use text_io;
with Standard_Mathematical_Functions;
with Double_Double_Constants;            use Double_Double_Constants;

package body DoblDobl_Mathematical_Functions is

  function SQRT ( x : double_float ) return double_double is

    dd_x : constant double_double := Create(x);

  begin
    return SQRT(dd_x);
  end SQRT;

  function SQRT ( x : double_double ) return double_double is

  -- Use Karp's trick: if x is an approximation to sqrt(a), then
  --   sqrt(a) = a*x + [a - (a*x)^2] * x / 2   (approx)
  -- The approximation is accurate to twice the accuracy of x.
  -- Also, the multiplication (a*x) and [-]*x can be done with
  -- only half the precision.

    res : double_double;
    wrk,tmp : double_float;

  begin
    if is_zero(x) then
      res := Create(0.0);
    elsif is_negative(x) then  -- original code returns NaN
      res := Create(-1.0);
    else
      wrk := Standard_Mathematical_Functions.SQRT(hi_part(x));
      wrk := 1.0/wrk;
      tmp := hi_part(x)*wrk;
      res := x - sqr(tmp);
      res := res*wrk;
      res := res*0.5;
      res := res + tmp;
    end if;
    return res;
  end SQRT;

-- AUXILIARIES FOR TRIG FUNCTIONS :

  function sin_taylor ( x : double_double ) return double_double is

  -- DESCRIPTION :
  --   Assuming |x| <= pi/32, computes sin(x) using Taylor series.

    res : double_double;
    thresh : constant double_float := 0.5*abs(to_double(x))*dd_eps;
    r,t,y : double_double;
    i : integer;

  begin
    if is_zero(x) then
      res := Create(0.0);
    else
      y := -sqr(x); res := x; r := x;
      i := 0;
      loop
        r := r*y;
        t := r*i_fac(i);
        res := res + t;
        i := i + 2;
        exit when (i >= n_inv_fact);
        exit when (abs(to_double(t)) <= thresh);
      end loop;
    end if;
    return res;
  end sin_taylor;

  function cos_taylor ( x : double_double ) return double_double is

  -- DESCRIPTION :
  --   Computes sin(x) using Taylor series.

    res : double_double;
    thresh : constant double_float := 0.5*dd_eps;
    r,t,y : double_double;
    i : integer;

  begin
    if is_zero(x) then
      res := Create(1.0);
    else
      y := -sqr(x); r := y;
      res := 1.0 + mul_pwr2(r,0.5);
      i := 1;
      loop
        r := r*y;
        t := r*i_fac(i);
        res := res + t;
        i := i + 2;
        exit when (i >= n_inv_fact);
        exit when (abs(to_double(t)) <= thresh);
      end loop;
    end if;
    return res;
  end cos_taylor;

  procedure sincos_taylor
              ( x : in double_double; sin_x,cos_x : out double_double ) is

  -- DESCRIPTION :
  --   Returns in sin_x and cos_x the Taylor approximations of sine
  --   and cosine of x.

  begin
    if is_zero(x) then
      sin_x := Create(0.0);
      cos_x := Create(1.0);
    else
      sin_x := sin_taylor(x);
      cos_x := SQRT(1.0 - sqr(sin_x));
    end if;
  end sincos_taylor;

  procedure Reduce_Modulo_2pi
              ( x : in double_double; t : out double_double;
                j,k,abs_k : out integer; fail : out boolean ) is

  -- DESCRIPTION :
  --   Reduces x modulo 2*pi, modulo pi/2, and then modulo pi/16.
  --   If this reduction does not work, then an error message is printed
  --   and fail is true on return.

  -- ON ENTRY :
  --   x        some double double number.

  -- ON RETURN :
  --   t        what remains of x;
  --   j        result after reduction modulo pi/2;
  --   k        result after reduction modulo pi/16;
  --   abs_k    absolute value of k;
  --   fail     true if reduction fails.

    z,r : double_double;
    q : double_float;

  begin
    z := nint(x/twopi);   -- aproximately reduce modulo 2*pi
    r := x - twopi*z;
   -- approximately reduce modulo pi/2 and then modulo pi/16
    q := double_float'floor(hi_part(r)/hi_part(pi2) + 0.5);
    t := r - pi2*q;
    j := integer(q);
    if j < -2 or j > 2 then
      put_line("dd_sin: cannot reduce modulo pi/2"); 
      fail := true;
    else
      q := double_float'floor(hi_part(t)/hi_part(pi16) + 0.5);
      t := t - pi16*q;
      k := integer(q);
      if k < 0
       then abs_k := -k;
       else abs_k := k;
      end if;
      if abs_k > 4 then
        put_line("dd_sin: cannot reduce modulo pi/16");
        fail := true;
      else
        fail := false;
      end if;
    end if;
  end Reduce_Modulo_2pi;

-- TARGET TRIG FUNCTIONS :

  function SIN ( x : double_double ) return double_double is

  -- To compute sin(x), we choose integers a, b so that
  -- x = s + a * (pi/2) + b * (pi/16) and |s| <= pi/32.
  -- Using sin(pi/16) = 0.5 * sqrt(2 - sqrt(2 + sqrt(2)))
  -- we can compute sin(x) from sin(s), cos(s).  This greatly
  -- increases the convergence of the sine Taylor series.

   res : double_double;

  begin
    if is_zero(x) then
      res := create(0.0);
    else
      declare
        t : double_double;
        j,k,abs_k : integer;
        fail : boolean;
      begin
        Reduce_Modulo_2pi(x,t,j,k,abs_k,fail);
        if fail then
          res := create(-2.0); -- instead of nan
        elsif k = 0 then
          case j is
            when 0 => res := sin_taylor(t);
            when 1 => res := cos_taylor(t);
            when -1 => res := -cos_taylor(t);
            when others => res := -sin_taylor(t);
          end case;
        else
          declare
            u,v,t_sin,t_cos : double_double;
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
      end;
    end if;
    return res;
  end SIN;

  function COS ( x : double_double ) return double_double is

   res : double_double;

  begin
    if is_zero(x) then
      res := create(1.0);
    else
      declare
        t : double_double;
        j,k,abs_k : integer;
        fail : boolean;
      begin
        Reduce_Modulo_2pi(x,t,j,k,abs_k,fail);
        if fail then
          res := create(-2.0); -- instead of nan
        elsif k = 0 then
          case j is
            when 0 => res := cos_taylor(t);
            when 1 => res := -sin_taylor(t);
            when -1 => res := sin_taylor(t);
            when others => res := -cos_taylor(t);
          end case;
        else
          declare
            u,v,t_sin,t_cos : double_double;
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
      end;
    end if;
    return res;
  end COS;

  procedure SINCOS ( x : in double_double; 
                     sin_x,cos_x : out double_double ) is

  begin
    if is_zero(x) then
      sin_x := create(0.0);
      cos_x := create(1.0);
    else
      declare
        t : double_double;
        j,k,abs_k : integer;
        fail : boolean;
      begin
        Reduce_Modulo_2pi(x,t,j,k,abs_k,fail);
        if fail then
          sin_x := create(-2.0);  -- instead of NaN
          cos_x := create(-2.0);
        else
          declare
            sin_t,cos_t,s,c,u,v : double_double;
          begin
            sincos_taylor(t,sin_t,cos_t);
            if abs_k = 0 then
              s := sin_t;
              c := cos_t;
            else
              u := cos_table(abs_k-1);
              v := sin_table(abs_k-1);
              if k > 0 then
                s := u*sin_t + v*cos_t;
                c := u*cos_t - v*sin_t;
              else
                s := u*sin_t - v*cos_t;
                c := u*cos_t + v*sin_t;
              end if;
            end if;
            case j is
              when 0 => sin_x := s; cos_x := c;
              when 1 => sin_x := c; cos_x := -s;
              when -1 => sin_x := -c; cos_x := s;
              when others => sin_x := -s; cos_x := -c;
            end case;
          end;
        end if;
      end;
    end if;
  end SINCOS;

  function TAN ( x : double_double ) return double_double is

    s,c : double_double;

  begin
    sincos(x,s,c);
    return s/c;
  end TAN;

  function ARCTAN ( x : double_double ) return double_double is

    one : constant double_double := create(1.0);

  begin
    return ARCTAN(x,one);
  end ARCTAN;

  function ARCTAN ( y,x : double_double ) return double_double is

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
        put_line("dd_arctan: both arguments are zero");
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
       else return -pi4;
      end if;
    else
      declare
        r : constant double_double := SQRT(sqr(x) + sqr(y));
        xx : constant double_double := x/r;
        yy : constant double_double := y/r;
        f : constant double_float 
          := Standard_Mathematical_Functions.Angle(to_double(y),to_double(x));
        z : double_double := create(f);
        sin_z,cos_z : double_double;
      begin
        sincos(z,sin_z,cos_z);
        if abs(hi_part(xx)) > abs(hi_part(yy)) then
         -- Use Newton iteration 1.  z' = z + (y - sin(z)) / cos(z)
          z := z + (yy - sin_z) / cos_z;
        else 
         -- Use Newton iteration 2.  z' = z - (x - cos(z)) / sin(z) 
          z := z - (xx - cos_z) / sin_z;
        end if;
        return z;
      end;
    end if;
  end ARCTAN;

  function ARCSIN ( x : double_double ) return double_double is
  begin
    if x < -1.0 or x > 1.0 then
      put_line("dd_arcsin: argument out of domain");
      return x; -- instead of NaN
    elsif is_one(x) then
      if is_positive(x)
       then return pi2;
       else return -pi2;
      end if;
    else
      declare
        one : constant double_double := create(1.0);
        s1x : constant double_double := SQRT(one - sqr(x));
      begin
        return ARCTAN(x,s1x);
      end;
    end if;
  end ARCSIN;

  function ARCCOS ( x : double_double ) return double_double is
  begin
    if x < -1.0 or x > 1.0 then
      put_line("dd_arccos: argument out of domain");
      return x; -- instead of NaN
    elsif is_one(x) then
      if is_positive(x)
       then return create(0.0);
       else return pi;
      end if;
    else
      declare
        one : constant double_double := create(1.0);
        s1x : constant double_double := SQRT(one - sqr(x));
      begin
        return ARCTAN(s1x,x);
      end;
    end if;
  end ARCCOS;

  function Radius ( x,y : double_double ) return double_double is
  begin
    return SQRT(sqr(x) + sqr(y));
  end Radius;

  function Angle ( x,y : double_double ) return double_double is
  begin
    return ARCTAN(x,y);
  end Angle;

end DoblDobl_Mathematical_Functions;
