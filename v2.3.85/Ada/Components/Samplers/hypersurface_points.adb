with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Numbers_Polar;    use Standard_Complex_Numbers_Polar;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Characters_and_Numbers;            use Characters_and_Numbers;
with Symbol_Table;                      use Symbol_Table;

package body Hypersurface_Points is

  function Affine_Eval ( p : Eval_Poly; b,v : Vector; t : Complex_Number )
                       return Complex_Number is

    point : constant Vector(b'range) := b + t*v;

  begin
    return Eval(p,point);
  end Affine_Eval;

  function Projective_Eval ( p : Eval_Poly; b,v : Vector; t : Complex_Number;
                             z : Vector ) return Complex_Number is

    point : Vector(1..2*z'last);

  begin
    point(b'range) := b + t*v;
    for i in z'range loop
      point(b'last+i) := z(i);
    end loop;
    return Eval(p,point); 
  end Projective_Eval;

  procedure Scale ( b,v : in Vector; t : in Complex_Number;
                    x : out Vector; z : in out Vector ) is

    absxi : double_float;

  begin
    x := b + t*v;
    for i in z'range loop
      absxi := Radius(x(i));
      if absxi > 1.0
       then x(i) := x(i)/absxi;
            z(i) := z(i)/absxi;
      end if; 
    end loop;
  end Scale;

  procedure Multihomogenization_Symbols ( n : in natural32 ) is
  begin
    Symbol_Table.Init(2*n);
    for i in 1..n loop
      declare
        xsb : Symbol;
        n1,n2 : character;
      begin
        xsb := (xsb'range => ' ');
        xsb(1) := 'x';
        if i < 10 
         then n1 := Convert_Decimal(i);
              xsb(2) := n1;
         else n1 := Convert_Decimal(i/10);
              xsb(2) := n1;
              n2 := Convert_Decimal(i mod 10);
              xsb(3) := n2;
        end if;
        Symbol_Table.Add(xsb);
      end;
    end loop;
    for i in 1..n loop
      declare
        zsb : Symbol;
        n1,n2 : character;
      begin
        zsb := (zsb'range => ' ');
        zsb(1) := 'z';
        if i < 10 
         then n1 := Convert_Decimal(i);
              zsb(2) := n1;
         else n1 := Convert_Decimal(i/10);
              zsb(2) := n1;
              n2 := Convert_Decimal(i mod 10);
              zsb(3) := n2;
        end if;
        Symbol_Table.Add(zsb);
      end;
    end loop;
  end Multihomogenization_Symbols;

  function Multihomogenize ( n : natural32; t : Term;
                             deg : Standard_Natural_Vectors.Vector )
                           return Term is

  -- DESCRIPTION :
  --   Returns the n-homogeneous term, where n is the number of variables
  --   and deg contains the degree of the polynomial in each variable.

    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..2*integer32(n));
    for i in 1..integer32(n) loop
      res.dg(i) := t.dg(i);
      res.dg(i+integer32(n)) := deg(i) - t.dg(i);
    end loop;
    return res;
  end Multihomogenize;

  function Multihomogenize ( n : natural32; p : Poly ) return Poly is

  -- DESCRIPTION :
  --   Returns the n-homogeneous polynomial.

    res : Poly := Null_Poly;
    deg : Standard_Natural_Vectors.Vector(1..integer32(n));

    procedure Multihomogenize_Term ( t : in Term; continue : out boolean ) is

      mht : Term := Multihomogenize(n,t,deg);

    begin
      Add(res,mht);
      continue := true;
    end Multihomogenize_Term;
    procedure Multihomogenize_Terms is
      new Visiting_Iterator(Multihomogenize_Term);

  begin
    for i in 1..integer32(n) loop
      deg(i) := natural32(Degree(p,i));
    end loop;
    Multihomogenize_Terms(p);
    return res;
  end Multihomogenize;

  procedure Affine_Newton_for_Line
                ( p : in Eval_Poly_Sys;
                  b,v : in Vector; t : in out Complex_Number;
                  ft,dt : out Complex_Number ) is

    x : constant Vector(b'range) := b + t*v;
    dp : Complex_Number := Create(0.0);

  begin
    ft := Eval(p(0),x);
    for i in 1..p'last loop                -- apply the chain rule
      dp := dp + v(i)*Eval(p(i),x);
    end loop;
    dt := -ft/dp;
    t := t + dt;
  end Affine_Newton_for_Line;

  procedure Projective_Newton_for_Line
                ( p : in Eval_Poly_Sys;
                  b,v : in Vector; t : in out Complex_Number;
                  z : in Vector; ft,dt : out Complex_Number ) is

    x : Vector(1..2*b'last);
    dp : Complex_Number := Create(0.0);

  begin
    x(b'range) := b + t*v;
    for i in z'range loop
      x(b'last+i) := z(i);
    end loop;
    ft := Eval(p(0),x);
    for i in 1..p'last loop                -- apply the chain rule
      dp := dp + v(i)*Eval(p(i),x);
    end loop;
    dt := -ft/dp;
    t := t + dt;
  end Projective_Newton_for_Line;

  procedure Affine_Newton_for_Surface
                ( p,q : in Eval_Poly_Sys; s : in Complex_Number;
                  b,v : in Vector; t : in out Complex_Number;
                  ft,dt : out Complex_Number ) is

    x : Vector(b'range) := b + t*v;
    dp : Complex_Number := Create(0.0);
    one_min_s : Complex_Number := Create(1.0) - s;
    dpq : Complex_Number;

  begin
    ft := Eval(p(0),x)*one_min_s + Eval(q(0),x)*s;
    for i in 1..p'last loop
      dpq := Eval(p(i),x)*one_min_s + Eval(q(i),x)*s;
      dp := dp + v(i)*dpq;
    end loop;
    dt := -ft/dp;
    t := t + dt;
  end Affine_Newton_for_Surface;

  procedure Projective_Newton_for_Surface
                ( p,q : in Eval_Poly_Sys; s : in Complex_Number;
                  b,v : in Vector; t : in out Complex_Number;
                  z : in Vector; ft,dt : out Complex_Number ) is

    x : Vector(1..2*b'last);
    dp : Complex_Number := Create(0.0);
    one_min_s : Complex_Number := Create(1.0) - s;
    dpq : Complex_Number;

  begin
    x(b'range) := b + t*v;
    for i in z'range loop
      x(b'last+i) := z(i);
    end loop;
    ft := Eval(p(0),x)*one_min_s + Eval(q(0),x)*s;
    for i in 1..p'last loop
      dpq := Eval(p(i),x)*one_min_s + Eval(q(i),x)*s;
      dp := dp + v(i)*dpq;
    end loop;
    dt := -ft/dp;
    t := t + dt;
  end Projective_Newton_for_Surface;

  procedure Affine_Correct_for_Line
                ( p : in Eval_Poly_Sys;
                  b,v : in Vector; t : in out Complex_Number;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean ) is

    ft,dt : Complex_Number;

  begin
    for k in 1..maxit loop
      Affine_Newton_for_Line(p,b,v,t,ft,dt);
      if (AbsVal(ft) < eps) or (AbsVal(dt) < eps)
       then numit := k; fail := false;
            return;
      end if;
    end loop;
    numit := maxit; fail := true;
  end Affine_Correct_for_Line;

  procedure Projective_Correct_for_Line
                ( p : in Eval_Poly_Sys;
                  b,v : in Vector; t : in out Complex_Number; z : in Vector;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean ) is

    ft,dt  : Complex_Number;

  begin
    for k in 1..maxit loop
      Projective_Newton_for_Line(p,b,v,t,z,ft,dt);
      if (AbsVal(ft) < eps) or (AbsVal(dt) < eps)
       then numit := k; fail := false;
            return;
      end if;
    end loop;
    numit := maxit; fail := true;
  end Projective_Correct_for_Line;

  procedure Affine_Correct_for_Line
                ( file : in file_type; p : in Eval_Poly_Sys;
                  b,v : in Vector; t : in out Complex_Number;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean ) is

    ft,dt : Complex_Number;
    aft,adt : double_float;

  begin
    for k in 1..maxit loop
      Affine_Newton_for_Line(p,b,v,t,ft,dt);
      aft := AbsVal(ft);
      adt := AbsVal(dt);
      put(file,"  Step "); put(file,k,1); put(file," :");
      put(file,"  |ft| : "); put(file,aft,3);
      put(file,"  |dt| : "); put(file,aft,3); new_line(file);
      if (aft < eps) or (adt < eps)
       then numit := k; fail := false;
            return;
      end if;
    end loop;
    numit := maxit; fail := true;
  end Affine_Correct_for_Line;

  procedure Projective_Correct_for_Line
                ( file : in file_type; p : in Eval_Poly_Sys;
                  b,v : in Vector; t : in out Complex_Number; z : in Vector;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean ) is

    ft,dt : Complex_Number;
    aft,adt : double_float;

  begin
    for k in 1..maxit loop
      Projective_Newton_for_Line(p,b,v,t,z,ft,dt);
      aft := AbsVal(ft);
      adt := AbsVal(dt);
      put(file,"  Step "); put(file,k,1); put(file," :");
      put(file,"  |ft| : "); put(file,aft,3);
      put(file,"  |dt| : "); put(file,aft,3); new_line(file);
      if (aft < eps) or (adt < eps)
       then numit := k; fail := false;
            return;
      end if;
    end loop;
    numit := maxit; fail := true;
  end Projective_Correct_for_Line;

  procedure Affine_Correct_for_Surface
                ( p,q : in Eval_Poly_Sys; s : in Complex_Number;
                  b,v : in Vector; t : in out Complex_Number;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean ) is

    ft,dt : Complex_Number;

  begin
    for k in 1..maxit loop
      Affine_Newton_for_Surface(p,q,s,b,v,t,ft,dt);
      if (AbsVal(ft) < eps) or (AbsVal(dt) < eps)
       then numit := k; fail := false;
            return;
      end if;
    end loop;
    numit := maxit; fail := true;
  end Affine_Correct_for_Surface;

  procedure Projective_Correct_for_Surface
                ( p,q : in Eval_Poly_Sys; s : in Complex_Number;
                  b,v : in Vector; t : in out Complex_Number; z : in Vector;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean ) is

    ft,dt : Complex_Number;

  begin
    for k in 1..maxit loop
      Projective_Newton_for_Surface(p,q,s,b,v,t,z,ft,dt);
      if (AbsVal(ft) < eps) or (AbsVal(dt) < eps)
       then numit := k; fail := false;
            return;
      end if;
    end loop;
    numit := maxit; fail := true;
  end Projective_Correct_for_Surface;

  procedure Affine_Correct_for_Surface
                ( file : in file_type;
                  p,q : in Eval_Poly_Sys; s : in Complex_Number;
                  b,v : in Vector; t : in out Complex_Number;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean ) is

    ft,dt : Complex_Number;
    aft,adt : double_float;

  begin
    for k in 1..maxit loop
      Affine_Newton_for_Surface(p,q,s,b,v,t,ft,dt);
      aft := AbsVal(ft);
      adt := AbsVal(dt);
      put(file,"  Step "); put(file,k,1); put(file," :");
      put(file,"  |ft| : "); put(file,aft,3);
      put(file,"  |dt| : "); put(file,aft,3); new_line(file);
      if (aft < eps) or (adt < eps)
       then numit := k; fail := false;
            return;
      end if;
    end loop;
    numit := maxit; fail := true;
  end Affine_Correct_for_Surface;

  procedure Projective_Correct_for_Surface
                ( file : in file_type;
                  p,q : in Eval_Poly_Sys; s : in Complex_Number;
                  b,v : in Vector; t : in out Complex_Number; z : in Vector;
                  maxit : in natural32; numit : out natural32;
                  eps : in double_float; fail : out boolean ) is

    ft,dt : Complex_Number;
    aft,adt : double_float;

  begin
    for k in 1..maxit loop
      Projective_Newton_for_Surface(p,q,s,b,v,t,z,ft,dt);
      aft := AbsVal(ft);
      adt := AbsVal(dt);
      put(file,"  Step "); put(file,k,1); put(file," :");
      put(file,"  |ft| : "); put(file,aft,3);
      put(file,"  |dt| : "); put(file,aft,3); new_line(file);
      if (aft < eps) or (adt < eps)
       then numit := k; fail := false;
            return;
      end if;
    end loop;
    numit := maxit; fail := true;
  end Projective_Correct_for_Surface;

  procedure Affine_Start_Hypersurface
                ( n,d,k : in natural32; q : out Poly; b,v,t : out Vector ) is

    tt : Term;
    one : Complex_Number := Create(1.0);

  begin
    tt.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    tt.cf := Create(-1.0);
    q := Create(tt);
    tt.dg(integer32(k)) := d;
    tt.cf := one;
    Add(q,tt);
    b := Random_Vector(1,integer32(n));
    b(integer32(k)) := Create(0.0);
    v := (1..integer32(n) => Create(0.0));
    v(integer32(k)) := Create(1.0);
    for i in 1..d loop
      t(integer32(i)) := Root(one,d,i);
    end loop;
  end Affine_Start_Hypersurface;

  procedure Degree_Start_Hypersurface
                ( deg : in Standard_Natural_Vectors.Vector; d : in natural32;
                  v : in Vector; q : out Poly; t : out Vector ) is

    t1,t2 : Term;
    lcf : Complex_Number := Create(1.0);
    arg : double_float;
    fd : constant double_float := double_float(d);

  begin
    t1.dg := new Standard_Natural_Vectors.Vector'(deg);
    t1.cf := Create(1.0);
    q := Create(t1);
    t2.dg := new Standard_Natural_Vectors.Vector'(deg'range => 0);
    t2.cf := Create(1.0);
    Sub(q,t2);
    Clear(t1); Clear(t2);
    for i in deg'range loop
      for j in 1..deg(i) loop
        Mul(lcf,v(i));
      end loop;
    end loop;
    arg := Angle(lcf);
   -- put("lead cff  = "); put(lcf); new_line;
   -- put("e^(i*arg) = "); put(Create(COS(arg),SIN(arg))); new_line;
    arg := -arg;
    for i in 1..d loop
      t(integer32(i)) := Create(COS(arg/fd),SIN(arg/fd));
      arg := arg + 2.0*PI;
    end loop;
  end Degree_Start_Hypersurface;

  procedure Step_Control ( fail : in boolean; maxstep : in double_float;
                           thres : in natural32; succnt : in out natural32;
                           step : in out double_float ) is

  -- DESCRIPTION :
  --   Implements a simple adaptive step control strategy.

  -- ON ENTRY :
  --   fail        if fail, then the step size will be multiplied with 0.75,
  --               otherwise, the step size will be multiplied with 1.5,
  --               provided the number of consecutive successes exceeds
  --               the thres value, also step <= maxstep;
  --   maxstep     maximal step size;
  --   thres       threshold to augment the step size;
  --   succnt      counts the number of successful corrections;
  --   step        current value of the step size.

  -- ON RETURN :
  --   succnt      updated counter;
  --   step        new step size. 

  begin
    if fail
     then succnt := 0;
          step := step*0.75;
     else succnt := succnt + 1;
          if succnt > thres
           then step := step*1.25;
                if step > maxstep
                 then step := maxstep;
                end if;
          end if;
    end if;
  end Step_Control;

  procedure Predict ( b0,v0,b1,v1 : in Vector; lambda : in Complex_Number;
                      b,v : out Vector ) is

  -- DESCRIPTION :
  --   Computes the line (b,v) := (b0,v0)*lambda + (b1,v1)*(1-lambda).

  -- ON ENTRY :
  --   b0             base point of line at the start;
  --   v0             direction of line at the start;
  --   b1             base point of line at the end;
  --   v1             direction of line at the end;
  --   lambda         continuation parameter moves from zero to one.

  -- ON RETURN :
  --   b              base point for new general line;
  --   v              direction for new general line.

    one_min_lambda : constant Complex_Number := Create(1.0) - lambda;

  begin
    b := lambda*b0 + one_min_lambda*b1;
    v := lambda*v0 + one_min_lambda*v1;
  end Predict;

  procedure Affine_Track_Moving_Line
                ( p : in Eval_Poly_Sys; b0,v0,b1,v1 : in Vector;
                  t : in out Complex_Number;
                  numbsteps : out natural32; fail : out boolean ) is

  -- NOTE :
  --   We use lambda as continuation parameter, going from 1 to 0.

    maxsteps : constant natural32 := 5000;
    maxstep : constant double_float := 0.05;
    maxit : constant natural32 := 4;
    eps : constant double_float := 1.0E-8;
    thres : constant natural32 := 3;
    succnt : natural32 := 0;
    step : double_float := maxstep;
    b,v : Vector(b0'range);
    numit : natural32;
    fail_newt : boolean := false;
    backup_t,t0 : Complex_Number := t;
    backup_lambda : Complex_Number := Create(1.0);
    lambda : Complex_Number := backup_lambda - step;

  begin
    backup_t := t;
    t0 := t;
    for i in 1..maxsteps loop
      Predict(b0,v0,b1,v1,lambda,b,v);
      if not fail_newt
       then t := t + Create(step)*(backup_t - t0);
      end if;
      Affine_Correct_for_Line(p,b,v,t,maxit,numit,eps,fail_newt);
      Step_Control(fail_newt,maxstep,thres,succnt,step);
      if fail_newt
       then t := backup_t;
            lambda := backup_lambda;
       else if (REAL_PART(lambda) = 0.0)
             then numbsteps := i; fail := false;
                  return;
             else t0 := backup_t;
                  backup_t := t;
                  backup_lambda := lambda;
                  lambda := lambda - step;
                  if REAL_PART(lambda) < 0.0
                   then lambda := Create(0.0);
                  end if;
            end if;
      end if;
    end loop;
    numbsteps := maxsteps; fail := true;
  end Affine_Track_Moving_Line;

  procedure Projective_Track_Moving_Line
                ( p : in Eval_Poly_Sys; b0,v0,b1,v1 : in Vector;
                  t : in out Complex_Number; z : in out Vector;
                  numbsteps : out natural32; fail : out boolean ) is

  -- NOTE :
  --   We use lambda as continuation parameter, going from 1 to 0.

    maxsteps : constant natural32 := 5000;
    maxstep : constant double_float := 0.05;
    maxit : constant natural32 := 4;
    eps : constant double_float := 1.0E-8;
    thres : constant natural32 := 3;
    succnt : natural32 := 0;
    step : double_float := maxstep;
    b,v,x,backup_z : Vector(b0'range);
    numit : natural32;
    fail_newt : boolean := false;
    backup_t,t0 : Complex_Number := t;
    backup_lambda : Complex_Number := Create(1.0);
    lambda : Complex_Number := backup_lambda - step;

  begin
    backup_t := t;
    backup_z := z;
    t0 := t;
    for i in 1..maxsteps loop
      Predict(b0,v0,b1,v1,lambda,b,v);
      if not fail_newt
       then t := t + Create(step)*(backup_t - t0);
      end if;
      Scale(b,v,t,x,z);
      Projective_Correct_for_Line(p,b,v,t,z,maxit,numit,eps,fail_newt);
      Step_Control(fail_newt,maxstep,thres,succnt,step);
      if fail_newt
       then t := backup_t;
            z := backup_z;
            lambda := backup_lambda;
       else if (REAL_PART(lambda) = 0.0)
             then numbsteps := i; fail := false;
                  return;
             else t0 := backup_t;
                  backup_t := t;
                  backup_z := z;
                  backup_lambda := lambda;
                  lambda := lambda - step;
                  if REAL_PART(lambda) < 0.0
                   then lambda := Create(0.0);
                  end if;
            end if;
      end if;
    end loop;
    numbsteps := maxsteps; fail := true;
  end Projective_Track_Moving_Line;

  procedure Affine_Track_Moving_Line
                ( file : in file_type;
                  p : in Eval_Poly_Sys; b0,v0,b1,v1 : in Vector;
                  t : in out Complex_Number;
                  numbsteps : out natural32; fail : out boolean ) is

  -- NOTE :
  --   We use lambda as continuation parameter, going from 1 to 0.

    maxsteps : constant natural32 := 5000;
    maxstep : constant double_float := 0.05;
    maxit : constant natural32 := 4;
    eps : constant double_float := 1.0E-8;
    thres : constant natural32 := 3;
    succnt : natural32 := 0;
    step : double_float := maxstep;
    b,v : Vector(b0'range);
    numit : natural32;
    fail_newt : boolean := false;
    backup_t,t0 : Complex_Number := t;
    backup_lambda : Complex_Number := Create(1.0);
    lambda : Complex_Number := backup_lambda - step;

  begin
    backup_t := t;
    t0 := t;
    for i in 1..maxsteps loop
      Predict(b0,v0,b1,v1,lambda,b,v);
      if not fail_newt
       then t := t + Create(step)*(backup_t - t0);
      end if;
      put(file,"Step "); put(file,i,1);
      put(file,"  with lambda = "); put(file,lambda); new_line(file);
      Affine_Correct_for_Line(file,p,b,v,t,maxit,numit,eps,fail_newt);
      Step_Control(fail_newt,maxstep,thres,succnt,step);
      if fail_newt
       then t := backup_t;
            lambda := backup_lambda;
       else if (REAL_PART(lambda) = 0.0)
             then numbsteps := i; fail := false;
                  return;
             else t0 := backup_t;
                  backup_t := t;
                  backup_lambda := lambda;
                  lambda := lambda - step;
                  if REAL_PART(lambda) < 0.0
                   then lambda := Create(0.0);
                  end if;
            end if;
      end if;
    end loop;
    numbsteps := maxsteps; fail := true;
  end Affine_Track_Moving_Line;

  procedure Projective_Track_Moving_Line
                ( file : in file_type;
                  p : in Eval_Poly_Sys; b0,v0,b1,v1 : in Vector;
                  t : in out Complex_Number; z : in out Vector;
                  numbsteps : out natural32; fail : out boolean ) is

  -- NOTE :
  --   We use lambda as continuation parameter, going from 1 to 0.

    maxsteps : constant natural32 := 5000;
    maxstep : constant double_float := 0.05;
    maxit : constant natural32 := 4;
    eps : constant double_float := 1.0E-8;
    thres : constant natural32 := 3;
    succnt : natural32 := 0;
    step : double_float := maxstep;
    b,v,x,backup_z : Vector(b0'range);
    numit : natural32;
    fail_newt : boolean := false;
    backup_t,t0 : Complex_Number := t;
    backup_lambda : Complex_Number := Create(1.0);
    lambda : Complex_Number := backup_lambda - step;

  begin
    backup_t := t;
    backup_z := z;
    t0 := t;
    for i in 1..maxsteps loop
      Predict(b0,v0,b1,v1,lambda,b,v);
      if not fail_newt
       then t := t + Create(step)*(backup_t - t0);
      end if;
      Scale(b,v,t,x,z);
      put(file,"Step "); put(file,i,1);
      put(file,"  with lambda = "); put(file,lambda); new_line(file);
      Projective_Correct_for_Line(file,p,b,v,t,z,maxit,numit,eps,fail_newt);
      Step_Control(fail_newt,maxstep,thres,succnt,step);
      if fail_newt
       then t := backup_t;
            z := backup_z;
            lambda := backup_lambda;
       else if (REAL_PART(lambda) = 0.0)
             then numbsteps := i; fail := false;
                  return;
             else t0 := backup_t;
                  backup_t := t;
                  backup_z := z;
                  backup_lambda := lambda;
                  lambda := lambda - step;
                  if REAL_PART(lambda) < 0.0
                   then lambda := Create(0.0);
                  end if;
            end if;
      end if;
    end loop;
    numbsteps := maxsteps; fail := true;
  end Projective_Track_Moving_Line;

  procedure Affine_Track_Moving_Surface
                ( p,q : in Eval_Poly_Sys; b,v : in Vector;
                  t : in out Complex_Number;
                  numbsteps : out natural32; fail : out boolean ) is

    maxsteps : constant natural32 := 5000;
    maxstep : constant double_float := 0.05;
    maxit : constant natural32 := 4;
    eps : constant double_float := 1.0E-8;
    thres : constant natural32 := 3;
    succnt : natural32 := 0;
    step : double_float := maxstep;
    numit : natural32;
    fail_newt : boolean := false;
    backup_t,t0 : Complex_Number := t;
    backup_lambda : Complex_Number := Create(1.0);
    lambda : Complex_Number := backup_lambda - step;

  begin
    backup_t := t;
    t0 := t;
    for i in 1..maxsteps loop
      if not fail_newt
       then t := t + Create(step)*(backup_t - t0);
      end if;
      Affine_Correct_for_Surface(p,q,lambda,b,v,t,maxit,numit,eps,fail_newt);
      Step_Control(fail_newt,maxstep,thres,succnt,step);
      if fail_newt
       then t := backup_t;
            lambda := backup_lambda;
       else if (REAL_PART(lambda) = 0.0)
             then numbsteps := i; fail := false;
                  return;
             else t0 := backup_t;
                  backup_t := t;
                  backup_lambda := lambda;
                  lambda := lambda - step;
                  if REAL_PART(lambda) < 0.0
                   then lambda := Create(0.0);
                  end if;
            end if;
      end if;
    end loop;
    numbsteps := maxsteps; fail := true;
  end Affine_Track_Moving_Surface;

  procedure Projective_Track_Moving_Surface
                ( p,q : in Eval_Poly_Sys; b,v : in Vector;
                  t : in out Complex_Number; z : in out Vector;
                  numbsteps : out natural32; fail : out boolean ) is

    maxsteps : constant natural32 := 5000;
    maxstep : constant double_float := 0.05;
    maxit : constant natural32 := 4;
    eps : constant double_float := 1.0E-8;
    thres : constant natural32 := 3;
    succnt : natural32 := 0;
    step : double_float := maxstep;
    numit : natural32;
    fail_newt : boolean := false;
    backup_t,t0 : Complex_Number := t;
    backup_lambda : Complex_Number := Create(1.0);
    lambda : Complex_Number := backup_lambda - step;
    x,backup_z : Vector(b'range);

  begin
    backup_t := t;
    backup_z := z;
    t0 := t;
    for i in 1..maxsteps loop
      if not fail_newt
       then t := t + Create(step)*(backup_t - t0);
      end if;
      Scale(b,v,t,x,z);
      Projective_Correct_for_Surface
        (p,q,lambda,b,v,t,z,maxit,numit,eps,fail_newt);
      Step_Control(fail_newt,maxstep,thres,succnt,step);
      if fail_newt
       then t := backup_t;
            z := backup_z;
            lambda := backup_lambda;
       else if (REAL_PART(lambda) = 0.0)
             then numbsteps := i; fail := false;
                  return;
             else t0 := backup_t;
                  backup_t := t;
                  backup_z := z;
                  backup_lambda := lambda;
                  lambda := lambda - step;
                  if REAL_PART(lambda) < 0.0
                   then lambda := Create(0.0);
                  end if;
            end if;
      end if;
    end loop;
    numbsteps := maxsteps; fail := true;
  end Projective_Track_Moving_Surface;

  procedure Affine_Track_Moving_Surface
                ( file : in file_type;
                  p,q : in Eval_Poly_Sys; b,v : in Vector;
                  t : in out Complex_Number;
                  numbsteps : out natural32; fail : out boolean ) is

    maxsteps : constant natural32 := 5000;
    maxstep : constant double_float := 0.05;
    maxit : constant natural32 := 4;
    eps : constant double_float := 1.0E-8;
    thres : constant natural32 := 3;
    succnt : natural32 := 0;
    step : double_float := maxstep;
    numit : natural32;
    fail_newt : boolean := false;
    backup_t,t0 : Complex_Number := t;
    backup_lambda : Complex_Number := Create(1.0);
    lambda : Complex_Number := backup_lambda - step;

  begin
    backup_t := t;
    t0 := t;
    for i in 1..maxsteps loop
      if not fail_newt
       then t := t + Create(step)*(backup_t - t0);
      end if;
      put(file,"Step "); put(file,i,1);
      put(file,"  with lambda = "); put(file,lambda); new_line(file);
      Affine_Correct_for_Surface
        (file,p,q,lambda,b,v,t,maxit,numit,eps,fail_newt);
      Step_Control(fail_newt,maxstep,thres,succnt,step);
      if fail_newt
       then t := backup_t;
            lambda := backup_lambda;
       else if (REAL_PART(lambda) = 0.0)
             then numbsteps := i; fail := false;
                  return;
             else t0 := backup_t;
                  backup_t := t;
                  backup_lambda := lambda;
                  lambda := lambda - step;
                  if REAL_PART(lambda) < 0.0
                   then lambda := Create(0.0);
                  end if;
            end if;
      end if;
    end loop;
    numbsteps := maxsteps; fail := true;
  end Affine_Track_Moving_Surface;

  procedure Projective_Track_Moving_Surface
                ( file : in file_type;
                  p,q : in Eval_Poly_Sys; b,v : in Vector;
                  t : in out Complex_Number; z : in out Vector;
                  numbsteps : out natural32; fail : out boolean ) is

    maxsteps : constant natural32 := 5000;
    maxstep : constant double_float := 0.05;
    maxit : constant natural32 := 4;
    eps : constant double_float := 1.0E-8;
    thres : constant natural32 := 3;
    succnt : natural32 := 0;
    step : double_float := maxstep;
    numit : natural32;
    fail_newt : boolean := false;
    backup_t,t0 : Complex_Number := t;
    backup_lambda : Complex_Number := Create(1.0);
    lambda : Complex_Number := backup_lambda - step;
    x,backup_z : Vector(b'range);

  begin
    backup_t := t;
    backup_z := z;
    t0 := t;
    for i in 1..maxsteps loop
      if not fail_newt
       then t := t + Create(step)*(backup_t - t0);
      end if;
      Scale(b,v,t,x,z);
      put(file,"Step "); put(file,i,1);
      put(file,"  with lambda = "); put(file,lambda); new_line(file);
      Projective_Correct_for_Surface
        (file,p,q,lambda,b,v,t,z,maxit,numit,eps,fail_newt);
      Step_Control(fail_newt,maxstep,thres,succnt,step);
      if fail_newt
       then t := backup_t;
            z := backup_z;
            lambda := backup_lambda;
       else if (REAL_PART(lambda) = 0.0)
             then numbsteps := i; fail := false;
                  return;
             else t0 := backup_t;
                  backup_t := t;
                  backup_z := z;
                  backup_lambda := lambda;
                  lambda := lambda - step;
                  if REAL_PART(lambda) < 0.0
                   then lambda := Create(0.0);
                  end if;
            end if;
      end if;
    end loop;
    numbsteps := maxsteps; fail := true;
  end Projective_Track_Moving_Surface;

  function Maximal_Affine_Residual
              ( p : Eval_Poly; b,v,roots : Vector ) return double_float is

    res : double_float := 0.0;
    eva : Complex_Number;
    abseva : double_float;

  begin
    for i in roots'range loop
      put(roots(i)); put(" : ");
      eva := Affine_Eval(p,b,v,roots(i));
      abseva := AbsVal(eva);
      put(abseva,3);
      new_line;
      if abseva > res
       then res := abseva;
      end if;
    end loop;
    return res;
  end Maximal_Affine_Residual;

  function Maximal_Projective_Residual
                ( p : Eval_Poly; b,v,roots : Vector; z : VecVec )
                return double_float is

    res : double_float := 0.0;
    eva : Complex_Number;
    abseva : double_float;

  begin
    for i in roots'range loop
      put(roots(i)); put(" : ");
      eva := Projective_Eval(p,b,v,roots(i),z(i).all);
      abseva := AbsVal(eva);
      put(abseva,3);
      new_line;
      put_line(z(i).all);
      if abseva > res
       then res := abseva;
      end if;
    end loop;
    return res;
  end Maximal_Projective_Residual;

end Hypersurface_Points;
