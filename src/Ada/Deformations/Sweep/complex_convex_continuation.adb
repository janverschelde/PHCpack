with Standard_Complex_Matrices;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Vector_Norms;     use DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Vector_Norms;     use QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Poly_Functions;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with Standard_IncFix_Continuation;
with DoblDobl_IncFix_Continuation;
with QuadDobl_IncFix_Continuation;

package body Complex_Convex_Continuation is

  function Interpolate ( a,b : Standard_Complex_Vectors.Vector;
                         t : Standard_Complex_Numbers.Complex_Number )
                       return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Numbers;

    res : Standard_Complex_Vectors.Vector(a'range);
    zero : constant Complex_Number := Create(0.0);
    one : Complex_Number := Create(1.0);

  begin
    if t = zero then
       res := a;
    elsif t = one then
       res := b;
    else
       one := one - t;
       for i in a'range loop
         res(i) := one*a(i) + t*b(i);  
       end loop;
    end if;
    return res;
  end Interpolate;

  function Interpolate ( a,b : DoblDobl_Complex_Vectors.Vector;
                         t : DoblDobl_Complex_Numbers.Complex_Number )
                       return DoblDobl_Complex_Vectors.Vector is

    use DoblDobl_Complex_Numbers;

    res : DoblDobl_Complex_Vectors.Vector(a'range);
    zero : constant Complex_Number := Create(integer(0));
    one : Complex_Number := Create(integer(1));

  begin
    if t = zero then
       res := a;
    elsif t = one then
       res := b;
    else
       one := one - t;
       for i in a'range loop
         res(i) := one*a(i) + t*b(i);  
       end loop;
    end if;
    return res;
  end Interpolate;

  function Interpolate ( a,b : QuadDobl_Complex_Vectors.Vector;
                         t : QuadDobl_Complex_Numbers.Complex_Number )
                       return QuadDobl_Complex_Vectors.Vector is

    use QuadDobl_Complex_Numbers;

    res : QuadDobl_Complex_Vectors.Vector(a'range);
    zero : constant Complex_Number := Create(integer(0));
    one : Complex_Number := Create(integer(1));

  begin
    if t = zero then
       res := a;
    elsif t = one then
       res := b;
    else
       one := one - t;
       for i in a'range loop
         res(i) := one*a(i) + t*b(i);  
       end loop;
    end if;
    return res;
  end Interpolate;

  function Circulate ( a,b : Standard_Complex_Vectors.Vector;
                       gamma,t : Standard_Complex_Numbers.Complex_Number )
                     return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Numbers;

    s : constant Complex_Number := t + gamma*t*(Create(1.0)-t);

  begin
    return Interpolate(a,b,s);
  end Circulate;

  function Circulate ( a,b : DoblDobl_Complex_Vectors.Vector;
                       gamma,t : DoblDobl_Complex_Numbers.Complex_Number )
                     return DoblDobl_Complex_Vectors.Vector is

    use DoblDobl_Complex_Numbers;

    s : constant Complex_Number := t + gamma*t*(Create(integer(1))-t);

  begin
    return Interpolate(a,b,s);
  end Circulate;

  function Circulate ( a,b : QuadDobl_Complex_Vectors.Vector;
                       gamma,t : QuadDobl_Complex_Numbers.Complex_Number )
                     return QuadDobl_Complex_Vectors.Vector is

    use QuadDobl_Complex_Numbers;

    s : constant Complex_Number := t + gamma*t*(Create(integer(1))-t);

  begin
    return Interpolate(a,b,s);
  end Circulate;

  function Differentiate ( a,b : Standard_Complex_Vectors.Vector )
                         return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Vectors;

    res : constant Vector(a'range) := b-a;

  begin
    return res;
  end Differentiate;

  function Differentiate ( a,b : DoblDobl_Complex_Vectors.Vector )
                         return DoblDobl_Complex_Vectors.Vector is

    use DoblDobl_Complex_Vectors;

    res : constant Vector(a'range) := b-a;

  begin
    return res;
  end Differentiate;

  function Differentiate ( a,b : QuadDobl_Complex_Vectors.Vector )
                         return QuadDobl_Complex_Vectors.Vector is

    use QuadDobl_Complex_Vectors;

    res : constant Vector(a'range) := b-a;

  begin
    return res;
  end Differentiate;

  procedure Standard_Reporting_Parameter_Continuation
              ( file : in file_type;
                n : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                pars : in Standard_Integer_Vectors.Vector;
                vars : in Standard_Integer_Vectors.Vector;
		sols : in out Standard_Complex_Solutions.Solution_List;
                output : in boolean ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Poly_Functions;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;
    use Standard_IncFix_Continuation;

    ep : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    ejm : Eval_Jaco_Mat(p'range,1..n) := Create(jm);
  
    function xvt ( x : Standard_Complex_Vectors.Vector;
                   t : Standard_Complex_Numbers.Complex_Number )
                 return Standard_Complex_Vectors.Vector is

    -- DESCRIPTION :
    --   Returns the value of all variables and parameters at
    --   the current value of the continuation parameter t.

      z : Standard_Complex_Vectors.Vector(1..n);
      v : constant Standard_Complex_Vectors.Vector(pars'range)
        := Evaluate_Parameters(t);

    begin
      for i in pars'range loop
        z(integer32(pars(i))) := v(i);
      end loop;
      for i in vars'range loop
        z(integer32(vars(i))) := x(i);
      end loop;
      return z;
    end xvt;

    function Eval ( x : Standard_Complex_Vectors.Vector;
                    t : Standard_Complex_Numbers.Complex_Number )
                  return Standard_Complex_Vectors.Vector is

      z : constant Standard_Complex_Vectors.Vector(1..n) := xvt(x,t);

    begin
      return Eval(ep,z);
    end Eval;

    function Diff ( x : Standard_Complex_Vectors.Vector;
                    t : Standard_Complex_Numbers.Complex_Number )
                  return Standard_Complex_Matrices.Matrix is

      res : Standard_Complex_Matrices.Matrix(p'range,x'range);
      z : constant Standard_Complex_Vectors.Vector(1..n) := xvt(x,t);

    begin
      for i in p'range loop
        for j in x'range loop
          res(i,j) := Eval(ejm(i,integer32(vars(j))),z);
        end loop;
      end loop;
      return res;
    end Diff;

    function Diff_t ( x : Standard_Complex_Vectors.Vector;
                      t : Standard_Complex_Numbers.Complex_Number )
                    return Standard_Complex_Vectors.Vector is

      res : Standard_Complex_Vectors.Vector(p'range)
          := (p'range => Create(0.0));
      z : constant Standard_Complex_Vectors.Vector(1..n) := xvt(x,t);
      d : constant Standard_Complex_Vectors.Vector
        := Differentiate_Parameters(t);

    begin
      for i in p'range loop
        for j in pars'range loop
          res(i) := res(i) + Eval(ejm(i,integer32(pars(j))),z)*d(j);
        end loop;
      end loop;
      return res;
    end Diff_t;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Eval,Diff_t,Diff);
    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,Eval,Diff_t,Diff);

  begin
    if output
     then Rep_Cont(file,sols,false,target=>Create(1.0));
     else Sil_Cont(sols,false,target=>Create(1.0));
    end if;
    Clear(ep); Clear(jm); Clear(ejm);
  end Standard_Reporting_Parameter_Continuation;

  procedure DoblDobl_Reporting_Parameter_Continuation
              ( file : in file_type;
                n : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pars : in Standard_Integer_Vectors.Vector;
                vars : in Standard_Integer_Vectors.Vector;
		sols : in out DoblDobl_Complex_Solutions.Solution_List;
                output : in boolean ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Poly_Functions;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_IncFix_Continuation;

    ep : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    ejm : Eval_Jaco_Mat(p'range,1..n) := Create(jm);
  
    function xvt ( x : DoblDobl_Complex_Vectors.Vector;
                   t : DoblDobl_Complex_Numbers.Complex_Number )
                 return DoblDobl_Complex_Vectors.Vector is

    -- DESCRIPTION :
    --   Returns the value of all variables and parameters at
    --   the current value of the continuation parameter t.

      z : DoblDobl_Complex_Vectors.Vector(1..n);
      v : constant DoblDobl_Complex_Vectors.Vector(pars'range)
        := Evaluate_Parameters(t);

    begin
      for i in pars'range loop
        z(integer32(pars(i))) := v(i);
      end loop;
      for i in vars'range loop
        z(integer32(vars(i))) := x(i);
      end loop;
      return z;
    end xvt;

    function Eval ( x : DoblDobl_Complex_Vectors.Vector;
                    t : DoblDobl_Complex_Numbers.Complex_Number )
                  return DoblDobl_Complex_Vectors.Vector is

      z : constant DoblDobl_Complex_Vectors.Vector(1..n) := xvt(x,t);

    begin
      return Eval(ep,z);
    end Eval;

    function Diff ( x : DoblDobl_Complex_Vectors.Vector;
                    t : DoblDobl_Complex_Numbers.Complex_Number )
                  return DoblDobl_Complex_Matrices.Matrix is

      res : DoblDobl_Complex_Matrices.Matrix(p'range,x'range);
      z : constant DoblDobl_Complex_Vectors.Vector(1..n) := xvt(x,t);

    begin
      for i in p'range loop
        for j in x'range loop
          res(i,j) := Eval(ejm(i,integer32(vars(j))),z);
        end loop;
      end loop;
      return res;
    end Diff;

    function Diff_t ( x : DoblDobl_Complex_Vectors.Vector;
                      t : DoblDobl_Complex_Numbers.Complex_Number )
                    return DoblDobl_Complex_Vectors.Vector is

      res : DoblDobl_Complex_Vectors.Vector(p'range)
          := (p'range => Create(integer(0)));
      z : constant DoblDobl_Complex_Vectors.Vector(1..n) := xvt(x,t);
      d : constant DoblDobl_Complex_Vectors.Vector
        := Differentiate_Parameters(t);

    begin
      for i in p'range loop
        for j in pars'range loop
          res(i) := res(i) + Eval(ejm(i,integer32(pars(j))),z)*d(j);
        end loop;
      end loop;
      return res;
    end Diff_t;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Eval,Diff_t,Diff);
    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,Eval,Diff_t,Diff);

  begin
    if output
     then Rep_Cont(file,sols,Create(integer(1)));
     else Sil_Cont(sols,Create(integer(1)));
    end if;
    Clear(ep); Clear(jm); Clear(ejm);
  end DoblDobl_Reporting_Parameter_Continuation;

  procedure QuadDobl_Reporting_Parameter_Continuation
              ( file : in file_type;
                n : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pars : in Standard_Integer_Vectors.Vector;
                vars : in Standard_Integer_Vectors.Vector;
		sols : in out QuadDobl_Complex_Solutions.Solution_List;
                output : in boolean ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Poly_Functions;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_IncFix_Continuation;

    ep : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    ejm : Eval_Jaco_Mat(p'range,1..n) := Create(jm);
  
    function xvt ( x : QuadDobl_Complex_Vectors.Vector;
                   t : QuadDobl_Complex_Numbers.Complex_Number )
                 return QuadDobl_Complex_Vectors.Vector is

    -- DESCRIPTION :
    --   Returns the value of all variables and parameters at
    --   the current value of the continuation parameter t.

      z : QuadDobl_Complex_Vectors.Vector(1..n);
      v : constant QuadDobl_Complex_Vectors.Vector(pars'range)
        := Evaluate_Parameters(t);

    begin
      for i in pars'range loop
        z(integer32(pars(i))) := v(i);
      end loop;
      for i in vars'range loop
        z(integer32(vars(i))) := x(i);
      end loop;
      return z;
    end xvt;

    function Eval ( x : QuadDobl_Complex_Vectors.Vector;
                    t : QuadDobl_Complex_Numbers.Complex_Number )
                  return QuadDobl_Complex_Vectors.Vector is

      z : constant QuadDobl_Complex_Vectors.Vector(1..n) := xvt(x,t);

    begin
      return Eval(ep,z);
    end Eval;

    function Diff ( x : QuadDobl_Complex_Vectors.Vector;
                    t : QuadDobl_Complex_Numbers.Complex_Number )
                  return QuadDobl_Complex_Matrices.Matrix is

      res : QuadDobl_Complex_Matrices.Matrix(p'range,x'range);
      z : constant QuadDobl_Complex_Vectors.Vector(1..n) := xvt(x,t);

    begin
      for i in p'range loop
        for j in x'range loop
          res(i,j) := Eval(ejm(i,integer32(vars(j))),z);
        end loop;
      end loop;
      return res;
    end Diff;

    function Diff_t ( x : QuadDobl_Complex_Vectors.Vector;
                      t : QuadDobl_Complex_Numbers.Complex_Number )
                    return QuadDobl_Complex_Vectors.Vector is

      res : QuadDobl_Complex_Vectors.Vector(p'range)
          := (p'range => Create(integer(0)));
      z : constant QuadDobl_Complex_Vectors.Vector(1..n) := xvt(x,t);
      d : constant QuadDobl_Complex_Vectors.Vector
        := Differentiate_Parameters(t);

    begin
      for i in p'range loop
        for j in pars'range loop
          res(i) := res(i) + Eval(ejm(i,integer32(pars(j))),z)*d(j);
        end loop;
      end loop;
      return res;
    end Diff_t;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Eval,Diff_t,Diff);
    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,Eval,Diff_t,Diff);

  begin
    if output
     then Rep_Cont(file,sols,Create(integer(1)));
     else Sil_Cont(sols,Create(integer(1)));
    end if;
    Clear(ep); Clear(jm); Clear(ejm);
  end QuadDobl_Reporting_Parameter_Continuation;

  procedure Standard_Silent_Parameter_Continuation
              ( n : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                pars : in Standard_Integer_Vectors.Vector;
                vars : in Standard_Integer_Vectors.Vector;
		sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Poly_Functions;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;
    use Standard_IncFix_Continuation;

    ep : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    ejm : Eval_Jaco_Mat(p'range,1..n) := Create(jm);
  
    function xvt ( x : Standard_Complex_Vectors.Vector;
                   t : Standard_Complex_Numbers.Complex_Number )
                 return Standard_Complex_Vectors.Vector is

    -- DESCRIPTION :
    --   Returns the value of all variables and parameters at
    --   the current value of the continuation parameter t.

      z : Standard_Complex_Vectors.Vector(1..n);
      v : constant Standard_Complex_Vectors.Vector(pars'range)
        := Evaluate_Parameters(t);

    begin
      for i in pars'range loop
        z(integer32(pars(i))) := v(i);
      end loop;
      for i in vars'range loop
        z(integer32(vars(i))) := x(i);
      end loop;
      return z;
    end xvt;

    function Eval ( x : Standard_Complex_Vectors.Vector;
                    t : Standard_Complex_Numbers.Complex_Number )
                  return Standard_Complex_Vectors.Vector is

      z : constant Standard_Complex_Vectors.Vector(1..n) := xvt(x,t);

    begin
      return Eval(ep,z);
    end Eval;

    function Diff ( x : Standard_Complex_Vectors.Vector;
                    t : Standard_Complex_Numbers.Complex_Number )
                  return Standard_Complex_Matrices.Matrix is

      res : Standard_Complex_Matrices.Matrix(p'range,x'range);
      z : constant Standard_Complex_Vectors.Vector(1..n) := xvt(x,t);

    begin
      for i in p'range loop
        for j in x'range loop
          res(i,j) := Eval(ejm(i,integer32(vars(j))),z);
        end loop;
      end loop;
      return res;
    end Diff;

    function Diff_t ( x : Standard_Complex_Vectors.Vector;
                      t : Standard_Complex_Numbers.Complex_Number )
                    return Standard_Complex_Vectors.Vector is

      res : Standard_Complex_Vectors.Vector(p'range)
          := (p'range => Create(0.0));
      z : constant Standard_Complex_Vectors.Vector(1..n) := xvt(x,t);
      d : constant Standard_Complex_Vectors.Vector
        := Differentiate_Parameters(t);

    begin
      for i in p'range loop
        for j in pars'range loop
          res(i) := res(i) + Eval(ejm(i,integer32(pars(j))),z)*d(j);
        end loop;
      end loop;
      return res;
    end Diff_t;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Eval,Diff_t,Diff);

  begin
    Sil_Cont(sols,false,target=>Create(1.0));
    Clear(ep); Clear(jm); Clear(ejm);
  end Standard_Silent_Parameter_Continuation;

  procedure DoblDobl_Silent_Parameter_Continuation
              ( n : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pars : in Standard_Integer_Vectors.Vector;
                vars : in Standard_Integer_Vectors.Vector;
		sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Poly_Functions;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_IncFix_Continuation;

    ep : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    ejm : Eval_Jaco_Mat(p'range,1..n) := Create(jm);
  
    function xvt ( x : DoblDobl_Complex_Vectors.Vector;
                   t : DoblDobl_Complex_Numbers.Complex_Number )
                 return DoblDobl_Complex_Vectors.Vector is

    -- DESCRIPTION :
    --   Returns the value of all variables and parameters at
    --   the current value of the continuation parameter t.

      z : DoblDobl_Complex_Vectors.Vector(1..n);
      v : constant DoblDobl_Complex_Vectors.Vector(pars'range)
        := Evaluate_Parameters(t);

    begin
      for i in pars'range loop
        z(integer32(pars(i))) := v(i);
      end loop;
      for i in vars'range loop
        z(integer32(vars(i))) := x(i);
      end loop;
      return z;
    end xvt;

    function Eval ( x : DoblDobl_Complex_Vectors.Vector;
                    t : DoblDobl_Complex_Numbers.Complex_Number )
                  return DoblDobl_Complex_Vectors.Vector is

      z : constant DoblDobl_Complex_Vectors.Vector(1..n) := xvt(x,t);

    begin
      return Eval(ep,z);
    end Eval;

    function Diff ( x : DoblDobl_Complex_Vectors.Vector;
                    t : DoblDobl_Complex_Numbers.Complex_Number )
                  return DoblDobl_Complex_Matrices.Matrix is

      res : DoblDobl_Complex_Matrices.Matrix(p'range,x'range);
      z : constant DoblDobl_Complex_Vectors.Vector(1..n) := xvt(x,t);

    begin
      for i in p'range loop
        for j in x'range loop
          res(i,j) := Eval(ejm(i,integer32(vars(j))),z);
        end loop;
      end loop;
      return res;
    end Diff;

    function Diff_t ( x : DoblDobl_Complex_Vectors.Vector;
                      t : DoblDobl_Complex_Numbers.Complex_Number )
                    return DoblDobl_Complex_Vectors.Vector is

      res : DoblDobl_Complex_Vectors.Vector(p'range)
          := (p'range => Create(integer(0)));
      z : constant DoblDobl_Complex_Vectors.Vector(1..n) := xvt(x,t);
      d : constant DoblDobl_Complex_Vectors.Vector
        := Differentiate_Parameters(t);

    begin
      for i in p'range loop
        for j in pars'range loop
          res(i) := res(i) + Eval(ejm(i,integer32(pars(j))),z)*d(j);
        end loop;
      end loop;
      return res;
    end Diff_t;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Eval,Diff_t,Diff);

  begin
    Sil_Cont(sols,Create(integer(1)));
    Clear(ep); Clear(jm); Clear(ejm);
  end DoblDobl_Silent_Parameter_Continuation;

  procedure QuadDobl_Silent_Parameter_Continuation
              ( n : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pars : in Standard_Integer_Vectors.Vector;
                vars : in Standard_Integer_Vectors.Vector;
		sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Poly_Functions;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_IncFix_Continuation;

    ep : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    ejm : Eval_Jaco_Mat(p'range,1..n) := Create(jm);
  
    function xvt ( x : QuadDobl_Complex_Vectors.Vector;
                   t : QuadDobl_Complex_Numbers.Complex_Number )
                 return QuadDobl_Complex_Vectors.Vector is

    -- DESCRIPTION :
    --   Returns the value of all variables and parameters at
    --   the current value of the continuation parameter t.

      z : QuadDobl_Complex_Vectors.Vector(1..n);
      v : constant QuadDobl_Complex_Vectors.Vector(pars'range)
        := Evaluate_Parameters(t);

    begin
      for i in pars'range loop
        z(integer32(pars(i))) := v(i);
      end loop;
      for i in vars'range loop
        z(integer32(vars(i))) := x(i);
      end loop;
      return z;
    end xvt;

    function Eval ( x : QuadDobl_Complex_Vectors.Vector;
                    t : QuadDobl_Complex_Numbers.Complex_Number )
                  return QuadDobl_Complex_Vectors.Vector is

      z : constant QuadDobl_Complex_Vectors.Vector(1..n) := xvt(x,t);

    begin
      return Eval(ep,z);
    end Eval;

    function Diff ( x : QuadDobl_Complex_Vectors.Vector;
                    t : QuadDobl_Complex_Numbers.Complex_Number )
                  return QuadDobl_Complex_Matrices.Matrix is

      res : QuadDobl_Complex_Matrices.Matrix(p'range,x'range);
      z : constant QuadDobl_Complex_Vectors.Vector(1..n) := xvt(x,t);

    begin
      for i in p'range loop
        for j in x'range loop
          res(i,j) := Eval(ejm(i,integer32(vars(j))),z);
        end loop;
      end loop;
      return res;
    end Diff;

    function Diff_t ( x : QuadDobl_Complex_Vectors.Vector;
                      t : QuadDobl_Complex_Numbers.Complex_Number )
                    return QuadDobl_Complex_Vectors.Vector is

      res : QuadDobl_Complex_Vectors.Vector(p'range)
          := (p'range => Create(integer(0)));
      z : constant QuadDobl_Complex_Vectors.Vector(1..n) := xvt(x,t);
      d : constant QuadDobl_Complex_Vectors.Vector
        := Differentiate_Parameters(t);

    begin
      for i in p'range loop
        for j in pars'range loop
          res(i) := res(i) + Eval(ejm(i,integer32(pars(j))),z)*d(j);
        end loop;
      end loop;
      return res;
    end Diff_t;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Eval,Diff_t,Diff);

  begin
    Sil_Cont(sols,Create(integer(1)));
    Clear(ep); Clear(jm); Clear(ejm);
  end QuadDobl_Silent_Parameter_Continuation;

end Complex_Convex_Continuation;
