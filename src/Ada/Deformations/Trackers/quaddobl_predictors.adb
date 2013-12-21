--with text_io; use text_io;
--with integer_io; use integer_io;
--with Standard_Floating_Numbers_io; use Standard_Floating_Numbers_io;
--with Standard_Complex_Numbers_io; use Standard_Complex_Numbers_io;
--with Quad_Double_Numbers_io; use Quad_Double_Numbers_io;
--with QuadDobl_Complex_Numbers_io;  use QuadDobl_Complex_Numbers_io;

with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Quad_Double_Constants;              use Quad_Double_Constants;
with QuadDobl_Mathematical_Functions;    use QuadDobl_Mathematical_Functions;
with Standard_Integer_Vectors;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with QUadDobl_Extrapolators;             use QuadDobl_Extrapolators;

package body QuadDobl_Predictors is

-- PREDICTORS for t :

  procedure Real_Predictor
              ( t : in out Complex_Number; target : in Complex_Number;
                h,tol : in double_float; pow : in natural32 := 1;
                hh : out double_float ) is

    nt : Complex_Number;
    dd_h,one,dpw,inv_dpw : quad_double;

  begin
    if pow = 1 then
      dd_h := Create(h);
      nt := t + dd_h;
    else 
      dpw := Create(pow);
      one := Create(1.0);
      inv_dpw := one/dpw;
      nt := Create((h+REAL_PART(t)**integer(pow))**inv_dpw);
    end if;
    if REAL_PART(nt) >= REAL_PART(target)
      or else abs( REAL_PART(nt) - REAL_PART(target) ) <= tol
     then dd_h := REAL_PART(target) - REAL_PART(t);
          t := target;
          hh := hihi_part(dd_h);
     else hh := h;
          t := nt;
    end if;
  --exception
  --  when others => put_line("Exception raised in Real_Predictor"); raise;
  end Real_Predictor;

  procedure Complex_Predictor
               ( t : in out Complex_Number; target : in Complex_Number;
                 h,tol : in double_float; hh : out double_float;
                 distance : in double_float; trial : in natural32 ) is

    nt,target_direction,step : Complex_Number;
    d,absnt,d_m_h : double_float;
    x,y,dre,dim,alfa : quad_double;
    zero : constant quad_double := Create(0.0);
    tr3 : natural32 range 0..2;

  begin
    target_direction := target - t;
    d := hihi_part(AbsVal(target_direction));
    if d < tol then
      nt := target - Create(distance);
      hh := hihi_part(AbsVal(nt - t));
    else
      tr3 := trial mod 3;
      case tr3 is
        when 0 =>
          target_direction := target_direction / Create(d);
          d_m_h := d-h; -- to resolve ambiguity in next > test
          if d_m_h > distance
           then step := Create(h) * target_direction;
           else step := (d - distance) * target_direction;
          end if;
          nt := t + step;
        when 1 =>
          dre := REAL_PART(target_direction);
          dim := IMAG_PART(target_direction);
          alfa := ARCTAN(dim/dre);
          x := h * COS(alfa + PI/4.0);
          y := h * SIN(alfa + PI/4.0);
          step := Create(x,y);
          nt := t + step;
          absnt := hihi_part(AbsVal(nt));
          if (absnt > AbsVal(target)) or else (AbsVal(nt - target) < distance)
           then step := Create(zero,y) * Create(h);
                nt := t + step;
          end if;
        when 2 =>
          dre := REAL_PART(target_direction);
          dim := IMAG_PART(target_direction);
          alfa := ARCTAN(dim/dre);
          x := h * COS(alfa - PI/4.0);
          y := h * SIN(alfa - PI/4.0);
          step := Create(x,y);
          nt := t + step;
          absnt := hihi_part(AbsVal(nt));
          if (absnt > AbsVal(target)) or else (AbsVal(nt - target) < distance)
           then step := Create(zero,-y) * Create(h);
                nt := t + step;
          end if;
      end case;
      hh := hihi_part(AbsVal(step));
    end if;
    t := nt;
  end Complex_Predictor;

  procedure Circular_Predictor
              ( t : in out Complex_Number; theta : in out double_float;
                t0_min_target,target : in Complex_Number; 
                h : in double_float ) is

    e_i_theta : Complex_Number;
    dd_theta : quad_double := create(theta);

  begin
    theta := theta + h;
    e_i_theta := Create(COS(dd_theta),SIN(dd_theta));
    t := target + t0_min_target * e_i_theta;
  end Circular_Predictor;

  procedure Geometric_Predictor
              ( t : in out Complex_Number; target : in Complex_Number;
                h,tol : in double_float ) is

    nt : Complex_Number;

  begin
    nt := target - Create(h)*(target - t);
    if abs( REAL_PART(nt) - REAL_PART(target) ) <= tol
     then t := target;
     else t := nt;
    end if;
  end Geometric_Predictor;

-- PREDICTORS FOR x :

  procedure Secant_Predictor
              ( x : in out Solution_Array; prev_x : in Solution_Array;
                fac : in Complex_Number; dist_x : in double_float ) is

    j : integer32;
    xx : Solution_Array(x'range);

  begin
    Copy(x,xx);
    for i in x'range loop
      x(i).v := x(i).v + fac * ( x(i).v - prev_x(i).v );
      j := xx'first;
      Equals(xx,x(i).v,i,dist_x,j);
      if j /= i
       then Copy(xx,x); exit;
      end if;
    end loop;
    Clear(xx);
  end Secant_Predictor;

-- PREDICTORS FOR t AND x :

  procedure Secant_Single_Real_Predictor
               ( x : in out Vector; prev_x : in Vector;
                 t : in out Complex_Number; prev_t,target : in Complex_Number;
                 h,tol : in double_float; pow : in natural32 := 1 ) is

    tmp : quad_double;
    hh : double_float;
    factor : Complex_Number;

  begin
    tmp := AbsVal(t - prev_t);
    Real_Predictor(t,target,h,tol,pow,hh);
    if tmp > tol
     then factor := Create(hh/tmp);
          x := x + factor * ( x - prev_x );
    end if;
  --exception 
  --  when others =>
  --    put_line("exception raised in secant_single_real_predictor");
  --    raise;
  end Secant_Single_Real_Predictor;

  procedure Quadratic_Single_Real_Predictor
                ( x : in out Vector; x1,x0 : in Vector;
                  t : in out Complex_Number;
                  t1,t0,target : in Complex_Number;
                  h,tol : in double_float; pow : in natural32 := 1 ) is

    tmp : quad_double;
    hh : double_float;
    t2 : constant Complex_Number := t;

  begin
    tmp := AbsVal(t - t1);
    Real_Predictor(t,target,h,tol,pow,hh);
    if tmp > tol then
      x := Extrapolate(t,t0,t1,t2,x0,x1,x);
    end if;
  end Quadratic_Single_Real_Predictor;

  procedure Cubic_Single_Real_Predictor
                ( x : in out Vector; x2,x1,x0 : in Vector;
                  t : in out Complex_Number;
                  t2,t1,t0,target : in Complex_Number;
                  h,tol : in double_float; pow : in natural32 := 1 ) is

    tmp : quad_double;
    hh : double_float;
    t3 : constant Complex_Number := t;

  begin
    tmp := AbsVal(t - t2);
    Real_Predictor(t,target,h,tol,pow,hh);
    if tmp > tol then
      x := Extrapolate(t,t0,t1,t2,t3,x0,x1,x2,x);
    end if;
  end Cubic_Single_Real_Predictor;

  procedure Quartic_Single_Real_Predictor
                ( x : in out Vector; x3,x2,x1,x0 : in Vector;
                  t : in out Complex_Number;
                  t3,t2,t1,t0,target : in Complex_Number;
                  h,tol : in double_float; pow : in natural32 := 1 ) is

    tmp : quad_double;
    hh : double_float;
    t4 : constant Complex_Number := t;

  begin
    tmp := AbsVal(t - t3);
    Real_Predictor(t,target,h,tol,pow,hh);
    if tmp > tol then
      x := Extrapolate(t,t0,t1,t2,t3,t4,x0,x1,x2,x3,x);
    end if;
  end Quartic_Single_Real_Predictor;

  procedure Quintic_Single_Real_Predictor
                ( x : in out Vector; x4,x3,x2,x1,x0 : in Vector;
                  t : in out Complex_Number;
                  t4,t3,t2,t1,t0,target : in Complex_Number;
                  h,tol : in double_float; pow : in natural32 := 1 ) is

    tmp : quad_double;
    hh : double_float;
    t5 : constant Complex_Number := t;

  begin
    tmp := AbsVal(t - t4);
    Real_Predictor(t,target,h,tol,pow,hh);
    if tmp > tol then
      x := Extrapolate(t,t0,t1,t2,t3,t4,t5,x0,x1,x2,x3,x4,x);
    end if;
  end Quintic_Single_Real_Predictor;

  procedure Secant_Multiple_Real_Predictor
               ( x : in out Solution_Array; prev_x : in Solution_Array;
                 t : in out Complex_Number; prev_t,target : in Complex_Number;
                 h,tol,dist_x : in double_float; pow : in natural32 := 1 ) is

    hh : double_float;
    dd_hh,tmp : quad_double;
    factor : Complex_Number;

  begin
    tmp := AbsVal(t - prev_t);
    Real_Predictor(t,target,h,tol,pow,hh);
    dd_hh := Create(hh);
    if tmp > tol
     then factor := Create(dd_hh/tmp);
          Secant_Predictor(x,prev_x,factor,dist_x);
    end if;
    for k in x'range loop
      x(k).t := t;
    end loop;
  end Secant_Multiple_Real_Predictor;

  procedure Secant_Single_Complex_Predictor
               ( x : in out Vector; prev_x : in Vector;
                 t : in out Complex_Number; prev_t,target : in Complex_Number;
                 h,tol,dist_t : in double_float; trial : in natural32 ) is

    hh : double_float;
    dd_hh,tmp : quad_double;
    factor : Complex_Number;

  begin
    tmp := AbsVal(t - prev_t);
    Complex_Predictor(t,target,h,tol,hh,dist_t,trial);
    dd_hh := Create(hh);
    if tmp > tol
     then factor := Create(dd_hh/tmp);
          x := x + factor * ( x - prev_x );
    end if;
  end Secant_Single_Complex_Predictor;

  procedure Secant_Multiple_Complex_Predictor
               ( x : in out Solution_Array; prev_x : in Solution_Array;
                 t : in out Complex_Number; prev_t,target : in Complex_Number;
                 h,tol,dist_x,dist_t : in double_float;
                 trial : in natural32 ) is

    hh : double_float;
    dd_hh,tmp : quad_double;
    factor : Complex_Number;

  begin
    tmp := AbsVal(t - prev_t);
    Complex_Predictor(t,target,h,tol,hh,dist_t,trial);
    dd_hh := Create(hh);
    if tmp > tol then
      factor := Create(dd_hh/tmp);
      Secant_Predictor(x,prev_x,factor,dist_x);
    end if;
    for k in x'range loop
      x(k).t := t;
    end loop;
  end Secant_Multiple_Complex_Predictor;

  procedure Secant_Circular_Predictor
               ( x : in out Vector; prev_x : in Vector;
                 t : in out Complex_Number; theta : in out double_float;
                 prev_t,t0_min_target,target : in Complex_Number;
                 h,tol : in double_float ) is

    factor : Complex_Number;
    tmp : constant quad_double := AbsVal(t-prev_t);

  begin
    if tmp <= tol then
      Circular_Predictor(t,theta,t0_min_target,target,h);
    else
      factor := Create(h/tmp);
      Circular_Predictor(t,theta,t0_min_target,target,h);
      x := x + factor * ( x - prev_x );
    end if;
  end Secant_Circular_Predictor;

  procedure Secant_Geometric_Predictor
               ( x : in out Vector; prev_x : in Vector;
                 t : in out Complex_Number; prev_t,target : in Complex_Number;
                 h,tol : in double_float ) is

    dist_prev,dist : quad_double;
    factor,tmp : Complex_Number;

  begin
    dist_prev := AbsVal(t - prev_t);     -- old stepsize
    tmp := t;
    Geometric_Predictor(t,target,h,tol);
    dist := AbsVal(t - tmp);             -- new stepsize
    if dist_prev > tol
     then factor := Create(dist/dist_prev);
          x := x + factor * ( x - prev_x );
    end if;
  end Secant_Geometric_Predictor;

  procedure Tangent_Single_Real_Predictor
               ( x : in out Vector; t : in out Complex_Number;
                 target : in Complex_Number; h,tol : in double_float;
                 pow : in natural32 := 1 ) is

    n : constant integer32 := x'last;
    info : integer32;
    hh : double_float;
    norm_tan,dd_hh : quad_double;
    rhs : Vector(x'range);
    j : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    prev_t : constant Complex_Number := t;

  begin
    Real_Predictor(t,target,h,tol,pow,hh);
    dd_hh := Create(hh);
    j := dH(x,prev_t);
    lufac(j,n,ipvt,info);
    if info = 0 then
      rhs := dH(x,prev_t);
      Min(rhs);
      lusolve(j,n,ipvt,rhs);
      norm_tan := Norm(rhs);
      if norm_tan > tol then
        dd_hh := dd_hh / norm_tan;
        Mul(rhs,Create(dd_hh));
        Add(x,rhs);
      end if;
    end if;
  end Tangent_Single_Real_Predictor;

  procedure Tangent_Multiple_Real_Predictor
               ( x : in out Solution_Array; t : in out Complex_Number;
                 target : in Complex_Number; h,tol,dist_x : in double_float;
                 nsys : in out natural32; pow : in natural32 := 1 ) is

    n : constant integer32 := x(x'first).n;
    hh : double_float;
    norm_tan,dd_hh : quad_double;
    rhs : Vector(1..n);
    info : integer32;
    j : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    xx : Solution_Array(x'range);
    jj : integer32;
    prev_t : constant Complex_Number := t;

  begin
    Real_Predictor(t,target,h,tol,pow,hh);
    dd_hh := Create(hh);
    Copy(x,xx);
    for i in x'range loop
      j := dH(x(i).v,prev_t);
      lufac(j,n,ipvt,info);
      if info = 0 then
        rhs := dH(x(i).v,prev_t);
        Min(rhs);
        lusolve(j,n,ipvt,rhs);
        nsys := nsys + 1;
        norm_tan := Norm(rhs);
        if norm_tan > tol then
          dd_hh := dd_hh / norm_tan;
          Mul(rhs,Create(dd_hh));
          Add(x(i).v,rhs);
        end if;
        jj := xx'first;
        Equals(xx,x(i).v,i,dist_x,jj);
        if jj /= i
         then Copy(xx,x); exit;
        end if;
      else
        Copy(xx,x); exit;
      end if;
    end loop;
    Clear(xx);
    for k in x'range loop
      x(k).t := t;
    end loop;
  end Tangent_Multiple_Real_Predictor;

  procedure Tangent_Single_Complex_Predictor 
               ( x : in out Vector; t : in out Complex_Number;
                 target : in Complex_Number;
                 h,tol,dist_t : in double_float; trial : in natural32 ) is

    n : constant integer32 := x'last;
    info : integer32;
    hh : double_float;
    norm_tan,dd_hh : quad_double;
    rhs : Vector(x'range);
    j : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    prev_t : constant Complex_Number := t;

  begin
    Complex_Predictor(t,target,h,tol,hh,dist_t,trial);
    dd_hh := Create(hh);
    j := dH(x,prev_t);
    lufac(j,n,ipvt,info);
    if info = 0 then
      rhs := dH(x,prev_t);
      Min(rhs);
      lusolve(j,n,ipvt,rhs);
      norm_tan := Norm(rhs);
      if norm_tan > tol then
        dd_hh := dd_hh / norm_tan;
        Mul(rhs,Create(dd_hh));
        Add(x,rhs);
      end if;
    end if;
  end Tangent_Single_Complex_Predictor;

  procedure Tangent_Multiple_Complex_Predictor
               ( x : in out Solution_Array; t : in out Complex_Number;
                 target : in Complex_Number;
                 h,tol,dist_x,dist_t : in double_float;
                 trial : in natural32; nsys : in out natural32 ) is

    n : constant integer32 := x(x'first).n;
    hh : double_float;
    norm_tan,dd_hh : quad_double;
    rhs : Vector(1..n);
    info : integer32;
    j : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    xx : Solution_Array(x'range);
    jj : integer32;
    prev_t : constant Complex_Number := t;

  begin
    Complex_Predictor(t,target,h,tol,hh,dist_t,trial);
    dd_hh := Create(hh);
    Copy(x,xx);
    for i in x'range loop
      j := dH(x(i).v,prev_t);
      lufac(j,n,ipvt,info);
      if info = 0 then
        rhs := dH(x(i).v,prev_t);
        Min(rhs);
        lusolve(j,n,ipvt,rhs);
        nsys := nsys + 1;
        norm_tan := Norm(rhs);
        if norm_tan > tol then
          dd_hh := dd_hh / norm_tan;
          Mul(rhs,Create(dd_hh));
          Add(x(i).v,rhs);
        end if;
        jj := xx'first;
        Equals(xx,x(i).v,i,dist_x,jj);
        if jj /= i
         then Copy(xx,x); exit;
        end if;
      else
        Copy(xx,x); exit;
      end if;
    end loop;
    Clear(xx);
    for k in x'range loop
      x(k).t := t;
    end loop;
  end Tangent_Multiple_Complex_Predictor;

  procedure Tangent_Circular_Predictor 
              ( x : in out Vector; t : in out Complex_Number;
                theta : in out double_float;
                t0_min_target,target : in Complex_Number;
                h,tol : in double_float ) is

    n : constant integer32 := x'last;
    info : integer32;
    norm_tan : quad_double;
    rhs : Vector(x'range);
    j : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);

  begin
    lufac(j,n,ipvt,info);
    if info = 0 then
      rhs := dH(x,t);
      Min(rhs);
      lusolve(j,n,ipvt,rhs);
    end if;
    Circular_Predictor(t,theta,t0_min_target,target,h);
    if info = 0 then
      norm_tan := Norm(rhs);
      if norm_tan > tol then
        Mul(rhs,Create(h/norm_tan));
        Add(x,rhs);
      end if;
    end if;
  end Tangent_Circular_Predictor;

  procedure Tangent_Geometric_Predictor
               ( x : in out Vector; t : in out Complex_Number;
                 target : in Complex_Number; h,tol : in double_float ) is

    n : constant integer32 := x'last;
    info : integer32;
    norm_tan,step : quad_double;
    rhs : Vector(x'range);
    j : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    prev_t : constant Complex_Number := t;

  begin
    Geometric_Predictor(t,target,h,tol);
    j := dH(x,prev_t);
    lufac(j,n,ipvt,info);
    if info = 0 then
      rhs := dH(x,prev_t);
      Min(rhs);
      lusolve(j,n,ipvt,rhs);
      norm_tan := Norm(rhs);
      if norm_tan > tol then
        step := AbsVal(t-prev_t) / norm_tan;
        Mul(rhs,Create(step));
        Add(x,rhs);
      end if;
    end if;
  end Tangent_Geometric_Predictor;

  function Hermite ( t0,t1,t,x0,x1,v0,v1 : Complex_Number )
                   return Complex_Number is

  -- DESCRIPTION :
  --   Returns the value of the third degree interpolating polynomial x(t),
  --   such that x(t0) = x0, x(t1) = x1, x'(t0) = v0 and x'(t1) = v1.

  -- REQUIRED : t0 /= t1.

  -- IMPLEMENTATION :
  --   x(t) = a3*t^3 + a2*t^2 + a1*t + a0,
  --   The four interpolation conditions lead to a linear system in
  --   the coefficients of x(t).  This system is first solved explicitly
  --   and then the polynomial x(t) is evaluated. 

    a0,a1,a2,a3,t10,v10 : Complex_Number;
    two : constant quad_double := Create(2.0);
    three : constant quad_double := Create(3.0);

  begin
    t10 := t1 - t0;
    v10 := (x1 - x0)/t10;
    a3 := (v1 + v0 - two*v10)/(t10**2);
    a2 := (v10 - v0 - (t1**2 + t1*t0 - two*t0**2)*a3)/t10;
    a1 := v0 - (three*a3*t0 + two*a2)*t0;
    a0 := x0 - ((a3*t0 + a2)*t0 + a1)*t0;
    return (((a3*t + a2)*t + a1)*t + a0);
 -- exception
 --   when others => put_line("Exception raised in Hermite 1"); raise;
  end Hermite;

  function Hermite ( t0,t1,t : Complex_Number; x0,x1,v0,v1 : Vector )
                   return Vector is

  -- DESCRIPTION :
  --   Returns the value of the third degree interpolating polynomial x(t),
  --   such that x(t0) = x0, x(t1) = x1, x'(t0) = v0 and x'(t1) = v1,
  --   for every component.

  -- REQUIRED : t0 /= t1.

    res : Vector(x0'range);

  begin
    for i in res'range loop
      res(i) := Hermite(t0,t1,t,x0(i),x1(i),v0(i),v1(i));
    end loop;
    return res;
 -- exception
 --   when others => put_line("Exception raised in Hermite 2"); raise;
  end Hermite;

  procedure Hermite_Single_Real_Predictor
                ( x : in out Vector; prev_x : in Vector;
                  t : in out Complex_Number; prev_t,target : in Complex_Number;
                  v : in out Vector; prev_v : in Vector;
                  h,tol : in double_float; pow : in natural32 := 1 ) is

    n : constant integer32 := x'last;
    info : integer32;
    hh : double_float;
    j : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    t1 : constant Complex_Number := t;

  begin
    Real_Predictor(t,target,h,tol,pow,hh);
    j := dH(x,t1);
    lufac(j,n,ipvt,info);
    if info = 0 then
      v := dH(x,t1);
      Min(v);
      lusolve(j,n,ipvt,v);
      if AbsVal(prev_t - t1) > 0.001 and AbsVal(t - t1) > 0.001
       then x := Hermite(prev_t,t1,t,prev_x,x,prev_v,v);
      end if;
    end if;
 -- exception
 --   when others => put_line("Exception raised in Hermite 3");
 --                  put("prev_t = "); put(prev_t); new_line;
 --                  put("t1 = "); put(t1); new_line;
 --                  put("t = "); put(t); new_line;
 --                  put("h = "); put(h); new_line; raise;
  end Hermite_Single_Real_Predictor;

end QuadDobl_Predictors;
