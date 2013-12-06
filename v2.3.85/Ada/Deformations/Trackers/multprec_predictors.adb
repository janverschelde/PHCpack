with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Mathematical_Functions;    use Multprec_Mathematical_Functions;
with Standard_Integer_Vectors;
with Multprec_Complex_Linear_Solvers;    use Multprec_Complex_Linear_Solvers;
with Multprec_Extrapolators;             use Multprec_Extrapolators;

package body Multprec_Predictors is

-- PREDICTORS for t :

  procedure Real_Predictor
              ( t : in out Complex_Number; target : in Complex_Number;
                h,tol : in Floating_Number; pow : in natural32 := 1;
                hh : out Floating_Number ) is

    nt : Complex_Number;

  begin
    --if pow = 1
    -- then
    nt := t + h;
    -- else nt := Create((h+REAL_PART(t)**pow)**(1.0/Floating_Number(pow)));
    --end if;
    if REAL_PART(nt) > REAL_PART(target)
      or else AbsVal( REAL_PART(nt) - REAL_PART(target) ) < tol
     then hh := REAL_PART(target) - REAL_PART(t);
          Copy(target,t);
     else Copy(h,hh);
          Clear(t); t := nt;
    end if;
  end Real_Predictor;

  procedure Complex_Predictor
               ( t : in out Complex_Number; target : in Complex_Number;
                 h,tol : in Floating_Number; hh : out Floating_Number;
                 distance : in Floating_Number; trial : in natural32 ) is

    nt,target_direction,step : Complex_Number;
    dre,dim,d,x,y,alfa,absnt,pi_d4,four,zero : Floating_Number;
    tr3 : natural32 range 0..2;

  begin
    target_direction := target - t;
    d := AbsVal(target_direction);
    if d < tol then
      nt := target - distance;
      step := nt - t;
      hh := AbsVal(step);
      Clear(step);
    else
      tr3 := trial mod 3;
      case tr3 is
        when 0 =>
          Div(target_direction,d);
          if (d - h) > distance
           then step := h * target_direction;
           else step := (d - distance) * target_direction;
          end if;
          nt := t + step;
        when 1 =>
          dre := REAL_PART(target_direction);
          dim := IMAG_PART(target_direction);
          alfa := ARCTAN(dim/dre);
          Clear(dre); Clear(dim);
          four := Create(integer(4));
          pi_d4 := PI/four;
          Add(alfa,pi_d4);
          x := h * COS(alfa);
          y := h * SIN(alfa);
          Clear(alfa); Clear(four); Clear(pi_d4);
          step := Create(x,y);
          nt := t + step;
          absnt := AbsVal(nt);
          if (absnt > AbsVal(target)) 
              or else (AbsVal(nt - target) < distance)
            then zero := Create(integer(0));
                 step := Create(zero,y) * Create(h);
                 Clear(zero);
                 Clear(nt); nt := t + step;
          end if;
          Clear(x); Clear(y);
        when 2 =>
          dre := REAL_PART(target_direction);
          dim := IMAG_PART(target_direction);
          alfa := ARCTAN(dim/dre);
          Clear(dre); Clear(dim);
          four := Create(integer(4));
          pi_d4 := PI/four;
          Sub(alfa,pi_d4);
          x := h * COS(alfa);
          y := h * SIN(alfa);
          Clear(alfa); Clear(four); Clear(pi_d4);
          step := Create(x,y);
          nt := t + step;
          absnt := AbsVal(nt);
          if (absnt > AbsVal(target)) 
              or else (AbsVal(nt - target) < distance)
           then zero := Create(integer(0));
                step := Create(zero,-y) * Create(h);
                Clear(zero);
                Clear(nt); nt := t + step;
          end if;
          Clear(x); Clear(y);
      end case;
      hh := AbsVal(step); Clear(step);
    end if;
    Clear(target_direction); Clear(d);
    Clear(t); t := nt;
  end Complex_Predictor;

  procedure Circular_Predictor
              ( t : in out Complex_Number; theta : in out Floating_Number;
                t0_min_target,target : in Complex_Number; 
                h : in Floating_Number ) is

    e_i_theta : Complex_Number;

  begin
    Add(theta,h);
    e_i_theta := Create(COS(theta),SIN(theta));
    t := target + t0_min_target * e_i_theta;
  end Circular_Predictor;

  procedure Geometric_Predictor
              ( t : in out Complex_Number; target : in Complex_Number;
                h,tol : in Floating_Number ) is

    nt,dif : Complex_Number;
    rep_nt,rep_tg : Floating_Number;

  begin
    dif := target - t;      -- nt := target - Create(h)*(target - t);
    Mul(dif,h);
    nt := target - dif;
    Clear(dif);
    rep_nt := REAL_PART(nt);
    rep_tg := REAL_PART(target);
    Sub(rep_nt,rep_tg);
    Clear(rep_tg);
    rep_tg := AbsVal(rep_nt);
    if rep_tg < tol    -- if AbsVal(REAL_PART(nt) - REAL_PART(target) ) < tol
     then Copy(target,t);
     else Clear(t); t := nt;
    end if;
    Clear(rep_tg); Clear(rep_nt);
  end Geometric_Predictor;

-- PREDICTORS FOR x :

  procedure Secant_Update ( x : in out Vector; prev_x : in Vector;
                            factor : in Complex_Number ) is

  -- DESCRIPTION :
  --   Updates the vector x as x := x + factor*(x-prev_x).

    inc : Vector(x'range) := x - prev_x;

  begin
    Mul(inc,factor);
    Add(x,inc);
    Clear(inc);
  end Secant_Update;

  procedure Secant_Predictor
              ( x : in out Solution_Array; prev_x : in Solution_Array;
                fac : in Complex_Number; dist_x : in Floating_Number ) is

    j : integer32;
    xx : Solution_Array(x'range);

  begin
    Copy(x,xx);
    for i in x'range loop
      Secant_Update(x(i).v,prev_x(i).v,fac);
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
                 h,tol : in Floating_Number; pow : in natural32 := 1 ) is

    hh,absdif : Floating_Number;
    factor,dif : Complex_Number;

  begin
    dif := t - prev_t;
    absdif := AbsVal(dif);
    Clear(dif);
    Real_Predictor(t,target,h,tol,pow,hh);
    if absdif > tol then
      Div(hh,absdif);
      factor := Create(hh);
      Secant_Update(x,prev_x,factor);
      Clear(factor);
    end if;
    Clear(hh); Clear(absdif);
  end Secant_Single_Real_Predictor;

  procedure Quadratic_Single_Real_Predictor
                ( x : in out Vector; x1,x0 : in Vector;
                  t : in out Complex_Number;
                  t1,t0,target : in Complex_Number;
                  h,tol : in Floating_Number; pow : in natural32 := 1 ) is

    hh,tmp : Floating_Number;
    t2 : Complex_Number;
    res : Vector(x'range);

  begin
    Copy(t,t2);
    tmp := AbsVal(t - t1);
    Real_Predictor(t,target,h,tol,pow,hh);
    if tmp > tol then
      res := Extrapolate(t,t0,t1,t2,x0,x1,x);
      Copy(res,x); Clear(res);
    end if;
    Clear(t2);
  end Quadratic_Single_Real_Predictor;

  procedure Cubic_Single_Real_Predictor
                ( x : in out Vector; x2,x1,x0 : in Vector;
                  t : in out Complex_Number;
                  t2,t1,t0,target : in Complex_Number;
                  h,tol : in Floating_Number; pow : in natural32 := 1 ) is

    hh,tmp : Floating_Number;
    t3 : Complex_Number;
    res : Vector(x'range);

  begin
    Copy(t,t3);
    tmp := AbsVal(t - t2);
    Real_Predictor(t,target,h,tol,pow,hh);
    if tmp > tol then
      res := Extrapolate(t,t0,t1,t2,t3,x0,x1,x2,x);
      Copy(res,x); Clear(res);
    end if;
    Clear(t3);
  end Cubic_Single_Real_Predictor;

  procedure Quartic_Single_Real_Predictor
                ( x : in out Vector; x3,x2,x1,x0 : in Vector;
                  t : in out Complex_Number;
                  t3,t2,t1,t0,target : in Complex_Number;
                  h,tol : in Floating_Number; pow : in natural32 := 1 ) is

    hh,tmp : Floating_Number;
    t4 : Complex_Number;
    res : Vector(x'range);

  begin
    Copy(t,t4);
    tmp := AbsVal(t - t3);
    Real_Predictor(t,target,h,tol,pow,hh);
    if tmp > tol then
      res := Extrapolate(t,t0,t1,t2,t3,t4,x0,x1,x2,x3,x);
      Copy(res,x); Clear(res);
    end if;
    Clear(t4);
  end Quartic_Single_Real_Predictor;

  procedure Quintic_Single_Real_Predictor
                ( x : in out Vector; x4,x3,x2,x1,x0 : in Vector;
                  t : in out Complex_Number;
                  t4,t3,t2,t1,t0,target : in Complex_Number;
                  h,tol : in Floating_Number; pow : in natural32 := 1 ) is

    hh,tmp : Floating_Number;
    t5 : Complex_Number;
    res : Vector(x'range);

  begin
    Copy(t,t5);
    tmp := AbsVal(t - t4);
    Real_Predictor(t,target,h,tol,pow,hh);
    if tmp > tol then
      res := Extrapolate(t,t0,t1,t2,t3,t4,t5,x0,x1,x2,x3,x4,x);
      Copy(res,x); Clear(res);
    end if;
    Clear(t5);
  end Quintic_Single_Real_Predictor;

  procedure Secant_Multiple_Real_Predictor
               ( x : in out Solution_Array; prev_x : in Solution_Array;
                 t : in out Complex_Number; prev_t,target : in Complex_Number;
                 h,tol,dist_x : in Floating_Number; pow : in natural32 := 1 ) is

    hh,absdif : Floating_Number;
    factor,dif : Complex_Number;

  begin
    dif := t - prev_t;
    absdif := AbsVal(dif);
    Clear(dif);
    Real_Predictor(t,target,h,tol,pow,hh);
    if absdif > tol then
      Div(hh,absdif);
      factor := Create(hh);
      Secant_Predictor(x,prev_x,factor,dist_x);
      Clear(factor);
    end if;
    Clear(hh);
    Clear(absdif);
    for k in x'range loop
      Copy(t,x(k).t);
    end loop;
  end Secant_Multiple_Real_Predictor;

  procedure Secant_Single_Complex_Predictor
               ( x : in out Vector; prev_x : in Vector;
                 t : in out Complex_Number; prev_t,target : in Complex_Number;
                 h,tol,dist_t : in Floating_Number; trial : in natural32 ) is

    hh,absdif : Floating_Number;
    factor,dif : Complex_Number;

  begin
    dif := t - prev_t;
    absdif := AbsVal(dif);
    Clear(dif);
    Complex_Predictor(t,target,h,tol,hh,dist_t,trial);
    if absdif > tol then
      Div(hh,absdif);
      factor := Create(hh);
      Secant_Update(x,prev_x,factor);
      Clear(factor);
    end if;
    Clear(absdif); Clear(hh);
  end Secant_Single_Complex_Predictor;

  procedure Secant_Multiple_Complex_Predictor
               ( x : in out Solution_Array; prev_x : in Solution_Array;
                 t : in out Complex_Number; prev_t,target : in Complex_Number;
                 h,tol,dist_x,dist_t : in Floating_Number;
                 trial : in natural32 ) is

    hh,absdif : Floating_Number;
    factor,dif : Complex_Number;

  begin
    dif := t - prev_t;
    absdif := AbsVal(dif);
    Clear(dif);
    Complex_Predictor(t,target,h,tol,hh,dist_t,trial);
    if absdif > tol then
      Div(hh,absdif);
      factor := Create(hh);
      Secant_Predictor(x,prev_x,factor,dist_x);
      Clear(factor);
    end if;
    Clear(absdif); Clear(hh);
    for k in x'range loop
      Copy(t,x(k).t);
    end loop;
  end Secant_Multiple_Complex_Predictor;

  procedure Secant_Circular_Predictor
               ( x : in out Vector; prev_x : in Vector;
                 t : in out Complex_Number; theta : in out Floating_Number;
                 prev_t,t0_min_target,target : in Complex_Number;
                 h,tol : in Floating_Number ) is

    factor,dif : Complex_Number;
    absdif : Floating_Number;

  begin
    dif := t - prev_t;
    absdif := AbsVal(dif);
    Clear(dif);
    if absdif < tol
     then Circular_Predictor(t,theta,t0_min_target,target,h);
     else factor := Create(h);
          Div(factor,absdif);
          Circular_Predictor(t,theta,t0_min_target,target,h);
          Secant_Update(x,prev_x,factor);
          Clear(factor);
    end if;
    Clear(absdif);
  end Secant_Circular_Predictor;

  procedure Secant_Geometric_Predictor
               ( x : in out Vector; prev_x : in Vector;
                 t : in out Complex_Number; prev_t,target : in Complex_Number;
                 h,tol : in Floating_Number ) is

    dist_prev,dist : Floating_Number;
    factor,cur_t,dif : Complex_Number;

  begin
    dif := t - prev_t;
    dist_prev := AbsVal(dif);                   -- old stepsize
    Clear(dif);
    Copy(t,cur_t);
    Geometric_Predictor(t,target,h,tol);
    dif := t - cur_t;
    Clear(cur_t);
    dist := AbsVal(dif);                        -- new stepsize
    Clear(dif);
    if dist_prev > tol then
      Div(dist,dist_prev);
      factor := Create(dist);
      Secant_Update(x,prev_x,factor);
      Clear(factor);
    end if;
    Clear(dist_prev);
    Clear(dist);
  end Secant_Geometric_Predictor;

  procedure Tangent_Single_Real_Predictor
               ( x : in out Vector; t : in out Complex_Number;
                 target : in Complex_Number; h,tol : in Floating_Number;
                 pow : in natural32 := 1 ) is

    n : constant integer32 := x'last;
    info : integer32;
    norm_tan,hh : Floating_Number;
    rhs : Vector(x'range);
    j : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    prev_t,fac : Complex_Number;

  begin
    Copy(t,prev_t);
    Real_Predictor(t,target,h,tol,pow,hh);
    j := dH(x,prev_t);
    lufac(j,n,ipvt,info);
    if info = 0 then
      rhs := dH(x,prev_t);
      Min(rhs);
      lusolve(j,n,ipvt,rhs);
      norm_tan := Norm(rhs);
      if norm_tan > tol then
        Div(hh,norm_tan);
        fac := Create(hh);
        Mul(rhs,fac);
        Clear(fac);
        Add(x,rhs);
      end if;
      Clear(norm_tan);
      Clear(rhs);
    end if;
    Clear(hh);
    Clear(j);
    Clear(prev_t);
  end Tangent_Single_Real_Predictor;

  procedure Tangent_Multiple_Real_Predictor
               ( x : in out Solution_Array; t : in out Complex_Number;
                 target : in Complex_Number; h,tol,dist_x : in Floating_Number;
                 nsys : in out natural32; pow : in natural32 := 1 ) is

    n : constant integer32 := x(x'first).n;
    norm_tan,hh : Floating_Number;
    rhs : Vector(1..n);
    info : integer32;
    j : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    xx : Solution_Array(x'range);
    jj : integer32;
    prev_t,fac : Complex_Number;

  begin
    Copy(t,prev_t);
    Real_Predictor(t,target,h,tol,pow,hh);
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
          Div(hh,norm_tan);
          fac := Create(hh);
          Mul(rhs,fac);
          Clear(fac);
          Add(x(i).v,rhs);
        end if;
        Clear(norm_tan);
        Clear(rhs);
        jj := xx'first;
        Equals(xx,x(i).v,i,dist_x,jj);
        if jj /= i
         then Copy(xx,x); exit;
        end if;
      else
        Clear(j); Copy(xx,x); exit;
      end if;
    end loop;
    Clear(xx);
    for k in x'range loop
      Copy(t,x(k).t);
    end loop;
    Clear(hh);
    Clear(prev_t);
  end Tangent_Multiple_Real_Predictor;

  procedure Tangent_Single_Complex_Predictor 
               ( x : in out Vector; t : in out Complex_Number;
                 target : in Complex_Number;
                 h,tol,dist_t : in Floating_Number; trial : in natural32 ) is

    n : constant integer32 := x'last;
    info : integer32;
    norm_tan,hh : Floating_Number;
    rhs : Vector(x'range);
    j : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    prev_t,fac : Complex_Number;

  begin
    Copy(t,prev_t);
    Complex_Predictor(t,target,h,tol,hh,dist_t,trial);
    j := dH(x,prev_t);
    lufac(j,n,ipvt,info);
    if info = 0 then
      rhs := dH(x,prev_t);
      Min(rhs);
      lusolve(j,n,ipvt,rhs);
      norm_tan := Norm(rhs);
      if norm_tan > tol then
        Div(hh,norm_tan);
        fac := Create(hh);
        Mul(rhs,fac);
        Add(x,rhs);
      end if;
      Clear(norm_tan);
      Clear(rhs);
    end if;
    Clear(j);
    Clear(prev_t);
  end Tangent_Single_Complex_Predictor;

  procedure Tangent_Multiple_Complex_Predictor
               ( x : in out Solution_Array; t : in out Complex_Number;
                 target : in Complex_Number;
                 h,tol,dist_x,dist_t : in Floating_Number;
                 trial : in natural32; nsys : in out natural32 ) is

    n : integer32 := x(x'first).n;
    norm_tan,hh : Floating_Number;
    rhs : Vector(1..n);
    info : integer32;
    j : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    xx : Solution_Array(x'range);
    jj : integer32;
    prev_t,fac : Complex_Number;

  begin
    Copy(t,prev_t);
    Complex_Predictor(t,target,h,tol,hh,dist_t,trial);
    Copy(x,xx);
    for i in x'range loop
      j := dH(x(i).v,prev_t);
      lufac(j,n,ipvt,info);
      if info = 0 then
        rhs := dH(x(i).v,prev_t);
        Min(rhs);
        lusolve(j,n,ipvt,rhs); Clear(j);
        nsys := nsys + 1;
        norm_tan := Norm(rhs);
        if norm_tan > tol then
          Div(hh,norm_tan);
          fac := Create(hh);
          Mul(rhs,fac);
          Clear(fac);
          Add(x(i).v,rhs);
        end if;
        jj := xx'first;
        Equals(xx,x(i).v,i,dist_x,jj);
        if jj /= i
         then Copy(xx,x); exit;
        end if;
      else
        Clear(j); Copy(xx,x); exit;
      end if;
    end loop;
    Clear(xx);
    for k in x'range loop
      Copy(t,x(k).t);
    end loop;
    Clear(prev_t);
  end Tangent_Multiple_Complex_Predictor;

  procedure Tangent_Circular_Predictor 
              ( x : in out Vector; t : in out Complex_Number;
                theta : in out Floating_Number;
                t0_min_target,target : in Complex_Number;
                h,tol : in Floating_Number ) is

    n : constant integer32 := x'last;
    info : integer32;
    norm_tan,refac : Floating_Number;
    fac : Complex_Number;
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
    Clear(j);
    Circular_Predictor(t,theta,t0_min_target,target,h);
    if info = 0 then
      norm_tan := Norm(rhs);
      if norm_tan > tol then
        refac := h/norm_tan;
        fac := Create(refac);
        Mul(rhs,fac);
        Clear(fac);
        Clear(refac);
        Add(x,rhs);
      end if;
      Clear(norm_tan);
      Clear(rhs);
    end if;
  end Tangent_Circular_Predictor;

  procedure Tangent_Geometric_Predictor
               ( x : in out Vector; t : in out Complex_Number;
                 target : in Complex_Number; h,tol : in Floating_Number ) is

    n : constant integer32 := x'last;
    info : integer32;
    rhs : Vector(x'range);
    j : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    prev_t,dif : Complex_Number;
    norm_tan,step : Floating_Number;

  begin
    Copy(t,prev_t);
    Geometric_Predictor(t,target,h,tol);
    j := dH(x,prev_t);
    lufac(j,n,ipvt,info);
    if info = 0 then
      rhs := dH(x,prev_t);
      Min(rhs);
      lusolve(j,n,ipvt,rhs);
      norm_tan := Norm(rhs);
      if norm_tan > tol then
        dif := t - prev_t;
        step := AbsVal(dif);
        Clear(dif);
        Div(step,norm_tan);
        dif := Create(step);
        Mul(rhs,dif);
        Clear(dif);
        Add(x,rhs);
      end if;
      Clear(norm_tan);
      Clear(rhs);
    end if;
    Clear(prev_t);
    Clear(j);
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

  begin
    t10 := t1 - t0;
    v10 := (x1 - x0)/t10;
    a3 := (v1 + v0 - Create(2.0)*v10)/(t10**2);
    a2 := (v10 - v0 - (t1**2 + t1*t0 - Create(2.0)*t0**2)*a3)/t10;
    a1 := v0 - (Create(3.0)*a3*t0 + Create(2.0)*a2)*t0;
    a0 := x0 - ((a3*t0 + a2)*t0 + a1)*t0;
    return (((a3*t + a2)*t + a1)*t + a0);
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
  end Hermite;

  procedure Hermite_Single_Real_Predictor
                ( x : in out Vector; prev_x : in Vector;
                  t : in out Complex_Number; prev_t,target : in Complex_Number;
                  v : in out Vector; prev_v : in Vector;
                  h,tol : in Floating_Number; pow : in natural32 := 1 ) is

    n : constant integer32 := x'last;
    info : integer32;
    j : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    t1,dif : Complex_Number;
    hh,absdif : Floating_Number;

  begin
    Copy(t,t1);
    Real_Predictor(t,target,h,tol,pow,hh);
    Clear(hh);
    j := dH(x,t1);
    lufac(j,n,ipvt,info);
    if info = 0 then
      v := dH(x,t1);
      Min(v);
      lusolve(j,n,ipvt,v);
      dif := prev_t - t1;
      absdif := AbsVal(dif);
      Clear(dif);
      if absdif > tol
       then x := Hermite(prev_t,t1,t,prev_x,x,prev_v,v);
      end if;
      Clear(absdif);
      Clear(v);
    end if;
    Clear(j); Clear(t1);
  end Hermite_Single_Real_Predictor;

end Multprec_Predictors;
