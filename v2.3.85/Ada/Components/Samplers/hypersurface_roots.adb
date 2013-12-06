with text_io;                           use text_io;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Polynomial_Roots;                  use Polynomial_Roots;

package body Hypersurface_Roots is

  procedure Substitute ( t : in Term; v : in Vector;
                         c : out Complex_Number; d : out natural32 ) is

  -- DESCRIPTION :
  --   Returns the term c*x^d after substituting x(i) by v(i)*x in t.

  begin
    c := t.cf;
    d := 0;
    for i in t.dg'range loop
      for j in 1..t.dg(i) loop
        Mul(c,v(i));
      end loop;
      d := d + t.dg(i);
    end loop;
  end Substitute;

  function Substitute ( p : Poly; v : Vector ) return Vector is

    d : constant integer32 := Degree(p);
    res : Vector(0..d) := (0..d => Create(0.0));

    procedure Substitute_Term ( t : in Term; cont : out boolean ) is

      ind : integer32;
      cff : Complex_Number;

    begin
      Substitute(t,v,cff,natural32(ind));
      res(ind) := res(ind) + cff;
      cont := true;
    end Substitute_Term;
    procedure Substitute_Terms is new Visiting_Iterator(Substitute_Term);

  begin
    Substitute_Terms(p);
    return res;
  end Substitute;

  procedure Scale ( v : in out Vector ) is

    nrm : double_float := Max_Norm(v);

  begin
    for i in v'range loop
      v(i) := v(i)/nrm;
    end loop;
  end Scale;

  function Homotopy ( p : Poly; v0,v1 : Vector; t : Complex_Number )
                    return Vector is

  -- DESCRIPTION :
  --   Returns the coefficient vector of the polynomial p(t*v0 + (1-t)*v1).

    v : Vector(v0'range) := t*v0;
    one_min_t : constant Complex_Number := Create(1.0)-t;
 
  begin
    v := v + one_min_t*v1;
   -- for i in v'range loop
   --   v(i) := v(i)/Radius(v(i));
   -- end loop;
    return Substitute(p,v);
  end Homotopy;

  procedure Affine_Path_Tracker
              ( d : in natural32; p : in Poly; v0,v1 : in Vector;
                s : in out Complex_Number; nbsteps : out natural32 ) is

    t : Complex_Number := Create(1.0);
    dt : double_float := 0.1;
    prev_s,prev_t : Complex_Number;
    h : Vector(0..integer32(d));
    dh : Vector(0..integer32(d)-1);
    nit,stepcnt : natural32 := 0;
    tol : constant double_float := 1.0E-8;
    fail : boolean;
    
  begin
    prev_s := s;
    prev_t := t;
    t := t - dt;
    for i in 1..10000 loop
      h := Homotopy(p,v0,v1,t);
     -- Scale(h);
      dh := Diff(h);
     -- put("Newton at t = "); put(t); new_line;
      Newton(h,dh,s,tol,4,nit,fail);
      stepcnt := stepcnt + 1;
      if fail then
        -- put_line("Newton failed to converge.");
        s := prev_s;
        t := prev_t;
        dt := dt/2.0;
        t := t - dt;
      else
        -- put_line("Newton succeeded.");
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
    h := Substitute(p,v1);
   -- Scale(h);
    dh := Diff(h);
    Newton(h,dh,s,1.0E-13,4,nit,fail);
    nbsteps := stepcnt;
  end Affine_Path_Tracker;

  procedure Projective_Path_Tracker
              ( d : in natural32; p : in Poly; v0,v1 : in Vector;
                s0,s1 : in out Complex_Number; nbsteps : out natural32 ) is

    t : Complex_Number := Create(1.0);
    dt : double_float := 0.1;
    prev_s0,prev_s1,prev_t : Complex_Number;
    h : Vector(0..integer32(d));
    dh : Vector(0..integer32(d)-1);
    nit,stepcnt : natural32 := 0;
    tol : constant double_float := 1.0E-8;
    fail : boolean;

  begin
    prev_s0 := s0;
    prev_s1 := s1;
    prev_t := t;
    t := t - dt;
    for i in 1..10000 loop
      Scale(s0,s1);
      h := Homotopy(p,v0,v1,t);
     -- Scale(h);
      Scale(h,s0);
      dh := Diff(h);
     -- put("Newton at t = "); put(t); new_line;
      Newton(h,dh,s1,tol,4,nit,fail);
      stepcnt := stepcnt + 1;
      if fail then
        --put_line("Newton failed to converge.");
        s0 := prev_s0;
        s1 := prev_s1;
        t := prev_t;
        dt := dt/2.0;
        t := t - dt;
      else
        --put_line("Newton succeeded.");
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
    h := Substitute(p,v1);
   -- Scale(h);
    Scale(h,s0);
    dh := Diff(h);
    Newton(h,dh,s1,1.0E-13,4,nit,fail);
    nbsteps := stepcnt;
  end Projective_Path_Tracker;

  procedure Affine_Track_Moving_Line
              ( p : in Poly; v0,v1 : in Vector; s : in out Vector ) is

    nbsteps : natural32 := 0;

  begin
    for i in s'range loop
      put("Tracking path "); put(i,1); put(" ... ");
      Affine_Path_Tracker(natural32(Degree(p)),p,v0,v1,s(i),nbsteps);
      put("done in "); put(nbsteps,1); put_line(" steps.");
    end loop;
  end Affine_Track_Moving_Line;

  procedure Projective_Track_Moving_Line
              ( p : in Poly; v0,v1 : in Vector; s0,s1 : in out Vector ) is

    nbsteps : natural32 := 0;

  begin
    for i in s0'range loop
      put("Tracking path "); put(i,1); put(" ... ");
      Projective_Path_Tracker
        (natural32(Degree(p)),p,v0,v1,s0(i),s1(i),nbsteps);
      put("done in "); put(nbsteps,1); put_line(" steps.");
    end loop;
  end Projective_Track_Moving_Line;

end Hypersurface_Roots;
