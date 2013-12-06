with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Numbers_Polar;     use Standard_Complex_Numbers_Polar;
with Standard_Random_Numbers;            use Standard_Random_Numbers;

package body Polynomial_Roots is

  function Eval ( p : Vector; x : Complex_Number ) return Complex_Number is

    y : Complex_Number := p(p'last);

  begin
    for i in reverse p'first+1..p'last loop
      y := y*x + p(i-1);
    end loop;
    return y;
  end Eval;

  function Diff ( p : Vector ) return Vector is

    res : Vector(0..p'last-1);
    ind : double_float;

  begin
    for i in 1..p'last loop
      ind := double_float(i);
      res(i-1) := Create(ind)*p(i);
    end loop;
    return res;
  end Diff;

  function Roots_of_Unity ( n : natural32 ) return Vector is

    res : Vector(1..integer32(n));
    one : constant Complex_Number := Create(1.0);

  begin
    for i in 1..n loop
      res(integer32(i)) := Root(one,n,i);
    end loop;
    return res;
  end Roots_of_Unity;

  procedure Scale ( s0,s1 : in out Complex_Number ) is

    abss1 : constant double_float := Radius(s1);

  begin
    if abss1 > 1.0 then
      s0 := s0/abss1;
      s1 := s1/abss1;
    end if;
  end Scale;

  procedure Scale ( p : in out Vector; s0 : in Complex_Number ) is

    d : constant integer32 := p'last;
    scaler : Complex_Number := s0;

  begin
    for i in reverse 0..d-1 loop
      p(i) := p(i)*scaler;
      scaler := scaler*s0;
    end loop;
  end Scale;

  function Scale ( p : Vector; s0 : Complex_Number ) return Vector is

    res : Vector(p'range) := p;
    d : constant integer32 := p'last;
    scaler : Complex_Number := s0;

  begin
    for i in reverse 0..d-1 loop
      res(i) := res(i)*scaler;
      scaler := scaler*s0;
    end loop;
    return res;
  end Scale;

  procedure Newton ( p,dp : in Vector; x : in out Complex_Number;
                     eps : in double_float; maxit : in natural32;
                     nit : out natural32; fail : out boolean ) is

    cnt : natural32 := 0;
    dx,y : Complex_Number;
    absdx,absy : double_float;

  begin
    fail := true;
    while cnt < maxit loop
      cnt := cnt+1;
      y := Eval(p,x);
      absy := AbsVal(y);
      dx := y/Eval(dp,x);
      x := x - dx;
      absdx := AbsVal(dx);
     -- put("  |dx| = "); put(absdx,3);
     -- put("  |p(x)| = "); put(absy,3); new_line;
      if ((absdx < eps) or (absy < eps))
       then fail := false; exit;
      end if;
    end loop;
    nit := cnt;
  end Newton;

  function Homotopy ( p,q : Vector; a,t : Complex_Number ) return Vector is

    res : Vector(p'range);
    one_min_t : constant Complex_Number := Create(1.0) - t;

  begin
    for i in p'range loop
      res(i) := p(i)*one_min_t + a*q(i)*t; 
    end loop;
    return res;
  end Homotopy;

  procedure Affine_Path_Tracker
               ( p,q : in Vector; a : in Complex_Number;
                 s : in out Complex_Number; nbsteps : out natural32 ) is

    t : Complex_Number := Create(1.0);
    dt : double_float := 0.1;
    prev_s,prev_t : Complex_Number;
    h : Vector(p'range);
    dh : Vector(0..p'last-1);
    nit,stepcnt : natural32 := 0;
    tol : constant double_float := 1.0E-8;
    fail : boolean;
    
  begin
    prev_s := s;
    prev_t := t;
    t := t - dt;
    loop
      h := Homotopy(p,q,a,t);
      dh := Diff(h);
     -- put("Newton at t = "); put(t); new_line;
      Newton(h,dh,s,tol,4,nit,fail);
      stepcnt := stepcnt + 1;
      if fail then --put_line("Newton failed to converge.");
        s := prev_s;
        t := prev_t;
        dt := dt/2.0;
        t := t - dt;
      else --put_line("Newton succeeded.");
        if REAL_PART(t) = 0.0
         then exit;
        end if;
        prev_s := s;
        prev_t := t;
        t := t - dt;
        dt := dt*1.5;
        if dt > 0.1
         then dt := 0.1;
        end if;
        if REAL_PART(t) < 0.0
         then t := Create(0.0);
        end if;
      end if;
    end loop;
    dh := Diff(p);
    Newton(p,dh,s,1.0E-13,4,nit,fail);
    nbsteps := stepcnt;
  end Affine_Path_Tracker;

  procedure Projective_Path_Tracker
               ( p,q : in Vector; a : in Complex_Number;
                 s0,s1 : in out Complex_Number; nbsteps : out natural32 ) is

    t : Complex_Number := Create(1.0);
    dt : double_float := 0.1;
    prev_s0,prev_s1,prev_t : Complex_Number;
    h : Vector(p'range);
    dh : Vector(0..p'last-1);
    nit,stepcnt : natural32 := 0;
    tol : constant double_float := 1.0E-8;
    fail : boolean;
    
  begin
    prev_s0 := s0;
    prev_s1 := s1;
    prev_t := t;
    t := t - dt;
    loop
      Scale(s0,s1);
      h := Homotopy(p,q,a,t);
      Scale(h,s0);
      dh := Diff(h);
     -- put("Newton at t = "); put(t); new_line;
      Newton(h,dh,s1,tol,4,nit,fail);
      stepcnt := stepcnt + 1;
      if fail then --put_line("Newton failed to converge.");
        s0 := prev_s0;
        s1 := prev_s1;
        t := prev_t;
        dt := dt/2.0;
        t := t - dt;
      else --put_line("Newton succeeded.");
        if REAL_PART(t) = 0.0
         then exit;
        end if;
        prev_s0 := s0;
        prev_s1 := s1;
        prev_t := t;
        t := t - dt;
        dt := dt*1.5;
        if dt > 0.1
         then dt := 0.1;
        end if;
        if REAL_PART(t) < 0.0
         then t := Create(0.0);
        end if;
      end if;
    end loop;
    Scale(s0,s1);
    h := Scale(p,s0);
    dh := Diff(h);
    Newton(h,dh,s1,1.0E-12,4,nit,fail);
    nbsteps := stepcnt;
  end Projective_Path_Tracker;

  procedure Affine_Continuation
              ( p,q : in Vector; a : in Complex_Number; s : in out Vector ) is

    nbsteps : natural32;

  begin
    for i in s'range loop
      Affine_Path_Tracker(p,q,a,s(i),nbsteps);
    end loop;
  end Affine_Continuation;

  procedure Affine_Continuation
              ( file : in file_type;
                p,q : in Vector; a : in Complex_Number; s : in out Vector ) is

    nbsteps : natural32;

  begin
    for i in s'range loop
      put(file,"Tracking path "); put(file,i,1); put(file," ... ");
      Affine_Path_Tracker(p,q,a,s(i),nbsteps);
      put(file," done "); put(file,nbsteps,1);
      put_line(file," predictor-corrector steps.");
    end loop;
  end Affine_Continuation;

  procedure Projective_Continuation
              ( p,q : in Vector; a : in Complex_Number;
                s0,s1 : in out Vector ) is

    nbsteps : natural32;

  begin
    for i in s1'range loop
      Projective_Path_Tracker(p,q,a,s0(i),s1(i),nbsteps);
    end loop;
  end Projective_Continuation;

  procedure Projective_Continuation
              ( file : in file_type;
                p,q : in Vector; a : in Complex_Number;
                s0,s1 : in out Vector ) is

    nbsteps : natural32;

  begin
    for i in s1'range loop
      put(file,"Tracking path "); put(file,i,1); put(file," ... ");
      Projective_Path_Tracker(p,q,a,s0(i),s1(i),nbsteps);
      put(file," done "); put(file,nbsteps,1);
      put_line(file," predictor-corrector steps.");
    end loop;
  end Projective_Continuation;

  procedure Affine_Solve ( p : in Vector; s : out Vector; 
                           maxres : out double_float ) is

    n : constant integer32 := p'last;
    q : Vector(p'range) := (0..n => Create(0.0));

  begin
    q(n) := Create(1.0);
    q(0) := -q(n);
    s := Roots_of_Unity(natural32(n));
    Affine_Continuation(p,q,Random1,s);
    maxres := Maximal_Residual(p,s);
  end Affine_Solve;

  procedure Affine_Solve ( file : in file_type; p : in Vector;
                           s : out Vector; maxres : out double_float ) is

    n : constant integer32 := p'last;
    q : Vector(p'range) := (0..n => Create(0.0));

  begin
    q(n) := Create(1.0);
    q(0) := -q(n);
    s := Roots_of_Unity(natural32(n));
    Affine_Continuation(p,q,Random1,s);
    Test_Affine_Roots(file,p,s,maxres);
  end Affine_Solve;

  procedure Projective_Solve 
               ( p : in Vector; s0,s1 : out Vector;
                 maxres : out double_float ) is

    n : constant integer32 := p'last;
    q : Vector(p'range) := (0..n => Create(0.0));

  begin
    q(n) := Create(1.0);
    q(0) := -q(n);
    s0 := (1..n => Create(1.0));
    s1 := Roots_of_Unity(natural32(n));
    Projective_Continuation(p,q,Random1,s0,s1);
    maxres := Maximal_Residual(p,s0,s1);
  end Projective_Solve;

  procedure Projective_Solve 
               ( file : in file_type; p : in Vector;
                 s0,s1 : out Vector; maxres : out double_float ) is

    n : constant integer32 := p'last;
    q : Vector(p'range) := (0..n => Create(0.0));

  begin
    q(n) := Create(1.0);
    q(0) := -q(n);
    s0 := (1..n => Create(1.0));
    s1 := Roots_of_Unity(natural32(n));
    Projective_Continuation(p,q,Random1,s0,s1);
    Test_Projective_Roots(file,p,s0,s1,maxres);
  end Projective_Solve;

  function Maximal_Residual ( p,s : Vector ) return double_float is

    eva : Complex_Number;
    abseva : double_float;
    maxres : double_float := 0.0;

  begin
    for i in s'range loop
      eva := Eval(p,s(i));
      abseva := AbsVal(eva);
      if abseva > maxres
       then maxres := abseva;
      end if;
    end loop;
    return maxres;
  end Maximal_Residual;

  function Maximal_Residual ( p,s0,s1 : Vector ) return double_float is

    sp : Vector(p'range);
    eva : Complex_Number;
    abseva : double_float;
    maxres : double_float := 0.0;

  begin
    for i in s1'range loop
      sp := Scale(p,s0(i));
      eva := Eval(sp,s1(i));
      abseva := AbsVal(eva);
      if abseva > maxres
       then maxres := abseva;
      end if;
    end loop;
    return maxres;
  end Maximal_Residual;

  procedure Test_Affine_Roots ( file : in file_type; p,s : in Vector;
                                maxres : out double_float ) is

    eva : Complex_Number;
    abseva : double_float;

  begin
    maxres := 0.0;
    for i in s'range loop
      put(file,s(i)); put(file," : ");
      eva := Eval(p,s(i));
      abseva := AbsVal(eva);
      put(file,abseva); new_line(file);
      if abseva > maxres
       then maxres := abseva;
      end if;
    end loop;
    put(file,"The maximal residual is ");
    put(file,maxres,3); new_line(file);
  end Test_Affine_Roots;

  procedure Test_Projective_Roots ( file : in file_type; p,s0,s1 : in Vector;
                                    maxres : out double_float ) is

    sp : Vector(p'range);
    eva : Complex_Number;
    abseva : double_float;

  begin
    maxres := 0.0;
    for i in s1'range loop
      put(file," [ "); put(file,s0(i),3); put(file," : ");
                       put(file,s1(i),3); put(file," ] :  ");
      sp := Scale(p,s0(i));
      eva := Eval(sp,s1(i));
      abseva := AbsVal(eva);
      put(file,abseva,3); new_line(file);
      if abseva > maxres
       then maxres := abseva;
      end if;
    end loop;
    put(file,"The maximal residual is ");
    put(file,maxres,3); new_line(file);
  end Test_Projective_Roots;

end Polynomial_Roots;
