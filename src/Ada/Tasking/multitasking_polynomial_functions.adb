with text_io;                        use text_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multitasking;                   use Multitasking;
with Standard_Polynomial_Flatteners;
with DoblDobl_Polynomial_Flatteners;
with QuadDobl_Polynomial_Flatteners;
with Multitasking_Matrix_x_Vector;   use Multitasking_Matrix_x_Vector;

package body Multitasking_Polynomial_Functions is

-- JOBS to EVALUATE THE MONOMIALS :

  procedure Silent_Eval_Job
             ( i,n : in integer32;
               v : in Standard_Integer_VecVecs.VecVec;
               x : in Standard_Complex_Vectors.Vector;
               y : out Standard_Complex_Vectors.Vector ) is
  begin
    for r in v'range loop
      if r mod n = i-1 then
        y(r) := Standard_Polynomial_Flatteners.Eval(v(r).all,x);
      end if;
    end loop;
  end Silent_Eval_Job;

  procedure Silent_Eval_Job
             ( i,n : in integer32;
               v : in Standard_Integer_VecVecs.VecVec;
               x : in DoblDobl_Complex_Vectors.Vector;
               y : out DoblDobl_Complex_Vectors.Vector ) is
  begin
    for r in v'range loop
      if r mod n = i-1 then
        y(r) := DoblDobl_Polynomial_Flatteners.Eval(v(r).all,x);
      end if;
    end loop;
  end Silent_Eval_Job;

  procedure Silent_Eval_Job
             ( i,n : in integer32;
               v : in Standard_Integer_VecVecs.VecVec;
               x : in QuadDobl_Complex_Vectors.Vector;
               y : out QuadDobl_Complex_Vectors.Vector ) is
  begin
    for r in v'range loop
      if r mod n = i-1 then
        y(r) := QuadDobl_Polynomial_Flatteners.Eval(v(r).all,x);
      end if;
    end loop;
  end Silent_Eval_Job;

  procedure Reporting_Eval_Job
             ( i,n : in integer32;
               v : in Standard_Integer_VecVecs.VecVec;
               x : in Standard_Complex_Vectors.Vector;
               y : out Standard_Complex_Vectors.Vector ) is
  begin
    put_line("hello from task " & to_string(i) & " evaluating monomials");
    for r in v'range loop
      if r mod n = i-1 then
        put_line("task " & to_string(i) & " computes " & to_string(r));
        y(r) := Standard_Polynomial_Flatteners.Eval(v(r).all,x);
      end if;
    end loop;
  end Reporting_Eval_Job;

  procedure Reporting_Eval_Job
             ( i,n : in integer32;
               v : in Standard_Integer_VecVecs.VecVec;
               x : in DoblDobl_Complex_Vectors.Vector;
               y : out DoblDobl_Complex_Vectors.Vector ) is
  begin
    put_line("hello from task " & to_string(i) & " evaluating monomials");
    for r in v'range loop
      if r mod n = i-1 then
        put_line("task " & to_string(i) & " computes " & to_string(r));
        y(r) := DoblDobl_Polynomial_Flatteners.Eval(v(r).all,x);
      end if;
    end loop;
  end Reporting_Eval_Job;

  procedure Reporting_Eval_Job
             ( i,n : in integer32;
               v : in Standard_Integer_VecVecs.VecVec;
               x : in QuadDobl_Complex_Vectors.Vector;
               y : out QuadDobl_Complex_Vectors.Vector ) is
  begin
    put_line("hello from task " & to_string(i) & " evaluating monomials");
    for r in v'range loop
      if r mod n = i-1 then
        put_line("task " & to_string(i) & " computes " & to_string(r));
        y(r) := QuadDobl_Polynomial_Flatteners.Eval(v(r).all,x);
      end if;
    end loop;
  end Reporting_Eval_Job;

-- EVALUATION OF THE MONOMIALS :

  function Silent_Eval
             ( n : integer32;
               v : Standard_Integer_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(v'range);

    procedure Job ( i,n : in integer32 ) is
    begin
      for r in v'range loop
        if r mod n = i-1 then
          res(r) := Standard_Polynomial_Flatteners.Eval(v(r).all,x);
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Silent_Workers(Job);

  begin
    do_jobs(n);
    return res;
  end Silent_Eval;

  function Silent_Eval
             ( n : integer32;
               v : Standard_Integer_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(v'range);

    procedure Job ( i,n : in integer32 ) is
    begin
      for r in v'range loop
        if r mod n = i-1 then
          res(r) := DoblDobl_Polynomial_Flatteners.Eval(v(r).all,x);
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Silent_Workers(Job);

  begin
    do_jobs(n);
    return res;
  end Silent_Eval;

  function Silent_Eval
             ( n : integer32;
               v : Standard_Integer_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(v'range);

    procedure Job ( i,n : in integer32 ) is
    begin
      for r in v'range loop
        if r mod n = i-1 then
          res(r) := QuadDobl_Polynomial_Flatteners.Eval(v(r).all,x);
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Silent_Workers(Job);

  begin
    do_jobs(n);
    return res;
  end Silent_Eval;

  function Reporting_Eval
             ( n : integer32;
               v : Standard_Integer_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(v'range);

    procedure Job ( i,n : in integer32 ) is
    begin
      put_line("hello from task " & to_string(i));
      for r in v'range loop
        if r mod n = i-1 then
          put_line("task " & to_string(i) & " computes " & to_string(r));
          res(r) := Standard_Polynomial_Flatteners.Eval(v(r).all,x);
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Silent_Workers(Job);

  begin
    do_jobs(n);
    return res;
  end Reporting_Eval;

  function Reporting_Eval
             ( n : integer32;
               v : Standard_Integer_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(v'range);

    procedure Job ( i,n : in integer32 ) is
    begin
      put_line("hello from task " & to_string(i));
      for r in v'range loop
        if r mod n = i-1 then
          put_line("task " & to_string(i) & " computes " & to_string(r));
          res(r) := DoblDobl_Polynomial_Flatteners.Eval(v(r).all,x);
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Silent_Workers(Job);

  begin
    do_jobs(n);
    return res;
  end Reporting_Eval;

  function Reporting_Eval
             ( n : integer32;
               v : Standard_Integer_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(v'range);

    procedure Job ( i,n : in integer32 ) is
    begin
      put_line("hello from task " & to_string(i));
      for r in v'range loop
        if r mod n = i-1 then
          put_line("task " & to_string(i) & " computes " & to_string(r));
          res(r) := QuadDobl_Polynomial_Flatteners.Eval(v(r).all,x);
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Silent_Workers(Job);

  begin
    do_jobs(n);
    return res;
  end Reporting_Eval;

-- EVALUATION OF THE DENSE FLATTENED SYSTEM :

  function Silent_Eval
             ( n : integer32;
               A : Standard_Complex_Matrices.Matrix;
               v : Standard_Integer_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    vx : constant Standard_Complex_Vectors.Vector := Silent_Eval(n,v,x);

  begin
    return Silent_Multiply(n,A,vx);
  end Silent_Eval;

  function Silent_Eval
             ( n : integer32;
               A : DoblDobl_Complex_Matrices.Matrix;
               v : Standard_Integer_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector is

    vx : constant DoblDobl_Complex_Vectors.Vector := Silent_Eval(n,v,x);

  begin
    return Silent_Multiply(n,A,vx);
  end Silent_Eval;

  function Silent_Eval
             ( n : integer32;
               A : QuadDobl_Complex_Matrices.Matrix;
               v : Standard_Integer_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    vx : constant QuadDobl_Complex_Vectors.Vector := Silent_Eval(n,v,x);

  begin
    return Silent_Multiply(n,A,vx);
  end Silent_Eval;

  function Reporting_Eval
             ( n : integer32;
               A : Standard_Complex_Matrices.Matrix;
               v : Standard_Integer_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    vx : constant Standard_Complex_Vectors.Vector := Reporting_Eval(n,v,x);

  begin
    return Reporting_Multiply(n,A,vx);
  end Reporting_Eval;

  function Reporting_Eval
             ( n : integer32;
               A : DoblDobl_Complex_Matrices.Matrix;
               v : Standard_Integer_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector is

    vx : constant DoblDobl_Complex_Vectors.Vector := Reporting_Eval(n,v,x);

  begin
    return Reporting_Multiply(n,A,vx);
  end Reporting_Eval;

  function Reporting_Eval
             ( n : integer32;
               A : QuadDobl_Complex_Matrices.Matrix;
               v : Standard_Integer_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    vx : constant QuadDobl_Complex_Vectors.Vector := Reporting_Eval(n,v,x);

  begin
    return Reporting_Multiply(n,A,vx);
  end Reporting_Eval;

-- EVALUATION OF THE SPARSE FLATTENED SYSTEM :

  function Silent_Eval 
             ( n : integer32;
               c : Standard_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(c'range);
    vx : constant Standard_Complex_Vectors.Vector := Silent_Eval(n,v,x);

    procedure Job ( i,n : in integer32 ) is
    begin
      for r in c'range loop
        if r mod n = i-1 then
          res(r) := Standard_Polynomial_Flatteners.Eval(c(r).all,vx,k(r).all);
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Silent_Workers(Job);

  begin
    for i in res'range loop
      res(i) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    do_jobs(n);
    return res;
  end Silent_Eval;

  function Silent_Eval 
             ( n : integer32;
               c : DoblDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(c'range);
    vx : constant DoblDobl_Complex_Vectors.Vector := Silent_Eval(n,v,x);

    procedure Job ( i,n : in integer32 ) is
    begin
      for r in c'range loop
        if r mod n = i-1 then
          res(r) := DoblDobl_Polynomial_Flatteners.Eval(c(r).all,vx,k(r).all);
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Silent_Workers(Job);

  begin
    for i in res'range loop
      res(i) := DoblDobl_Complex_Numbers.Create(integer32(0));
    end loop;
    do_jobs(n);
    return res;
  end Silent_Eval;

  function Silent_Eval 
             ( n : integer32;
               c : QuadDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(c'range);
    vx : constant QuadDobl_Complex_Vectors.Vector := Silent_Eval(n,v,x);

    procedure Job ( i,n : in integer32 ) is
    begin
      for r in c'range loop
        if r mod n = i-1 then
          res(r) := QuadDobl_Polynomial_Flatteners.Eval(c(r).all,vx,k(r).all);
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Silent_Workers(Job);

  begin
    for i in res'range loop
      res(i) := QuadDobl_Complex_Numbers.Create(integer32(0));
    end loop;
    do_jobs(n);
    return res;
  end Silent_Eval;

  function Reporting_Eval 
             ( n : integer32;
               c : Standard_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(c'range);
    vx : constant Standard_Complex_Vectors.Vector := Reporting_Eval(n,v,x);

    procedure Job ( i,n : in integer32 ) is
    begin
      put_line("hello from task " & to_string(i));
      for r in res'range loop
        if r mod n = i-1 then
          put_line("task " & to_string(i) & " computes " & to_string(r));
          res(r) := Standard_Polynomial_Flatteners.Eval(c(r).all,vx,k(r).all);
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Silent_Workers(Job);

  begin
    for i in res'range loop
      res(i) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    do_jobs(n);
    return res;
  end Reporting_Eval;

  function Reporting_Eval 
             ( n : integer32;
               c : DoblDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(c'range);
    vx : constant DoblDobl_Complex_Vectors.Vector := Reporting_Eval(n,v,x);

    procedure Job ( i,n : in integer32 ) is
    begin
      put_line("hello from task " & to_string(i));
      for r in res'range loop
        if r mod n = i-1 then
          put_line("task " & to_string(i) & " computes " & to_string(r));
          res(r) := DoblDobl_Polynomial_Flatteners.Eval(c(r).all,vx,k(r).all);
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Silent_Workers(Job);

  begin
    for i in res'range loop
      res(i) := DoblDobl_Complex_Numbers.Create(integer32(0));
    end loop;
    do_jobs(n);
    return res;
  end Reporting_Eval;

  function Reporting_Eval 
             ( n : integer32;
               c : QuadDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(c'range);
    vx : constant QuadDobl_Complex_Vectors.Vector := Reporting_Eval(n,v,x);

    procedure Job ( i,n : in integer32 ) is
    begin
      put_line("hello from task " & to_string(i));
      for r in res'range loop
        if r mod n = i-1 then
          put_line("task " & to_string(i) & " computes " & to_string(r));
          res(r) := QuadDobl_Polynomial_Flatteners.Eval(c(r).all,vx,k(r).all);
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Silent_Workers(Job);

  begin
    for i in res'range loop
      res(i) := QuadDobl_Complex_Numbers.Create(integer32(0));
    end loop;
    do_jobs(n);
    return res;
  end Reporting_Eval;

-- EVALUATION OF THE SPARSE FLATTENED SYSTEM with LOOPING Workers :

  function Silent_Looping_Eval 
             ( n : integer32;
               c : Standard_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(c'range);
    first : boolean_array(1..n) := (1..n => true);
    vx : Standard_Complex_Vectors.Vector(v'range);

    use Standard_Polynomial_Flatteners;

    procedure Job ( i,n : in integer32; continue : out boolean ) is
    begin
      if first(i) then
        Silent_Eval_Job(i,n,v,x,vx);
        first(i) := false;
        continue := true;
      else
        for r in c'range loop
          if r mod n = i-1 then
            res(r) := Eval(c(r).all,vx,k(r).all);
          end if;
        end loop;
        continue := false;
      end if;
    end Job;
    procedure do_jobs is new Silent_Looping_Workers(Job);

  begin
    for i in res'range loop
      res(i) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    do_jobs(n);
    return res;
  end Silent_Looping_Eval;

  function Silent_Looping_Eval 
             ( n : integer32;
               c : DoblDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(c'range);
    first : boolean_array(1..n) := (1..n => true);
    vx : DoblDobl_Complex_Vectors.Vector(v'range);

    use DoblDobl_Polynomial_Flatteners;

    procedure Job ( i,n : in integer32; continue : out boolean ) is
    begin
      if first(i) then
        Silent_Eval_Job(i,n,v,x,vx);
        first(i) := false;
        continue := true;
      else
        for r in c'range loop
          if r mod n = i-1 then
            res(r) := Eval(c(r).all,vx,k(r).all);
          end if;
        end loop;
        continue := false;
      end if;
    end Job;
    procedure do_jobs is new Silent_Looping_Workers(Job);

  begin
    for i in res'range loop
      res(i) := DoblDobl_Complex_Numbers.Create(integer32(0));
    end loop;
    do_jobs(n);
    return res;
  end Silent_Looping_Eval;

  function Silent_Looping_Eval 
             ( n : integer32;
               c : QuadDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(c'range);
    first : boolean_array(1..n) := (1..n => true);
    vx : QuadDobl_Complex_Vectors.Vector(v'range);

    use QuadDobl_Polynomial_Flatteners;

    procedure Job ( i,n : in integer32; continue : out boolean ) is
    begin
      if first(i) then
        Silent_Eval_Job(i,n,v,x,vx);
        first(i) := false;
        continue := true;
      else
        for r in c'range loop
          if r mod n = i-1 then
            res(r) := Eval(c(r).all,vx,k(r).all);
          end if;
        end loop;
        continue := false;
      end if;
    end Job;
    procedure do_jobs is new Silent_Looping_Workers(Job);

  begin
    for i in res'range loop
      res(i) := QuadDobl_Complex_Numbers.Create(integer32(0));
    end loop;
    do_jobs(n);
    return res;
  end Silent_Looping_Eval;

  function Reporting_Looping_Eval 
             ( n : integer32;
               c : Standard_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(c'range);
    first : boolean_array(1..n) := (1..n => true);
    vx : Standard_Complex_Vectors.Vector(v'range);

    use Standard_Polynomial_Flatteners;

    procedure Job ( i,n : in integer32; continue : out boolean ) is
    begin
      if first(i) then
        put_line("hello from task " & to_string(i) & " at first stage");
        Reporting_Eval_Job(i,n,v,x,vx);
        first(i) := false;
        continue := true;
      else
        put_line("hello from task " & to_string(i) & " at second stage");
        for r in res'range loop
          if r mod n = i-1 then
            put_line("task " & to_string(i) & " computes " & to_string(r));
            res(r) := Eval(c(r).all,vx,k(r).all);
          end if;
        end loop;
        continue := false;
      end if;
    end Job;
    procedure do_jobs is new Reporting_Looping_Workers(Job);

  begin
    for i in res'range loop
      res(i) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    do_jobs(n);
    return res;
  end Reporting_Looping_Eval;

  function Reporting_Looping_Eval 
             ( n : integer32;
               c : DoblDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(c'range);
    first : boolean_array(1..n) := (1..n => true);
    vx : DoblDobl_Complex_Vectors.Vector(v'range);

    use DoblDobl_Polynomial_Flatteners;

    procedure Job ( i,n : in integer32; continue : out boolean ) is
    begin
      if first(i) then
        put_line("hello from task " & to_string(i) & " at first stage");
        Reporting_Eval_Job(i,n,v,x,vx);
        first(i) := false;
        continue := true;
      else
        put_line("hello from task " & to_string(i) & " at second stage");
        for r in res'range loop
          if r mod n = i-1 then
            put_line("task " & to_string(i) & " computes " & to_string(r));
            res(r) := Eval(c(r).all,vx,k(r).all);
          end if;
        end loop;
        continue := false;
      end if;
    end Job;
    procedure do_jobs is new Reporting_Looping_Workers(Job);

  begin
    for i in res'range loop
      res(i) := DoblDobl_Complex_Numbers.Create(integer32(0));
    end loop;
    do_jobs(n);
    return res;
  end Reporting_Looping_Eval;

  function Reporting_Looping_Eval 
             ( n : integer32;
               c : QuadDobl_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(c'range);
    first : boolean_array(1..n) := (1..n => true);
    vx : QuadDobl_Complex_Vectors.Vector(v'range);

    use QuadDobl_Polynomial_Flatteners;

    procedure Job ( i,n : in integer32; continue : out boolean ) is
    begin
      if first(i) then
        put_line("hello from task " & to_string(i) & " at first stage");
        Reporting_Eval_Job(i,n,v,x,vx);
        first(i) := false;
        continue := true;
      else
        put_line("hello from task " & to_string(i) & " at second stage");
        for r in res'range loop
          if r mod n = i-1 then
            put_line("task " & to_string(i) & " computes " & to_string(r));
            res(r) := Eval(c(r).all,vx,k(r).all);
          end if;
        end loop;
        continue := false;
      end if;
    end Job;
    procedure do_jobs is new Reporting_Looping_Workers(Job);

  begin
    for i in res'range loop
      res(i) := QuadDobl_Complex_Numbers.Create(integer32(0));
    end loop;
    do_jobs(n);
    return res;
  end Reporting_Looping_Eval;

end Multitasking_Polynomial_Functions;
