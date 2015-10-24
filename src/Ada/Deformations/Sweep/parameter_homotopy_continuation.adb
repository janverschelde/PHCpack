with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Numbers_io;                        use Numbers_io;
with Standard_Natural_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Symbol_Table,Symbol_Table_io;      use Symbol_Table;
with Standard_Floating_Polynomials;
with Standard_Complex_to_Real_Poly;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Solution_Diagnostics;
with Standard_IncFix_Continuation;      use Standard_IncFix_Continuation;
with Drivers_for_Poly_Continuation;     use Drivers_for_Poly_Continuation;
with Standard_Root_Refiners;            use Standard_Root_Refiners;
with Standard_Parameter_Systems;        use Standard_Parameter_Systems;
with Standard_Parameter_Solutions;      use Standard_Parameter_Solutions;
with Standard_Quad_Parameters;
with Standard_Quad_Turn_Points_io;      use Standard_Quad_Turn_Points_io;
with Standard_Quad_Sweepers;            use Standard_Quad_Sweepers;

package body Parameter_Homotopy_Continuation is

-- TARGET ROUTINES :

  function Define_Start ( v : Standard_Complex_Vectors.Vector;
                          ip : Standard_Integer_Vectors.Vector )
                        return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(ip'range);

  begin
    for i in ip'range loop
      res(i) := v(integer32(ip(i)));
    end loop;
    return res;
  end Define_Start;

  function Define_Complex_Target
             ( ip : Standard_Integer_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(ip'range);
    rv,iv : double_float;

  begin
    put_line("Reading complex target values for the parameters...");
    for i in res'range loop
      put(" ");
      Symbol_Table_io.put(Symbol_Table.Get(natural32(ip(i))));
      put(" : ");
      Read_Double_Float(rv);
      Read_Double_Float(iv);
      res(i) := Standard_Complex_Numbers.Create(rv,iv);
    end loop;
    return res;
  end Define_Complex_Target;

  function Define_Real_Target
             ( ip : Standard_Integer_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(ip'range);
    f : double_float;

  begin
    put_line("Reading real target values for the parameters...");
    for i in res'range loop
      put(" ");
      Symbol_Table_io.put(Symbol_Table.Get(natural32(ip(i))));
      put(" : ");
      Read_Double_Float(f);
      res(i) := Standard_Complex_Numbers.Create(f);
    end loop;
    return res;
  end Define_Real_Target;

  procedure Write_Complex_Parameter_Values
             ( file : in file_type;
               labels : in Standard_Integer_Vectors.Vector;
               values : in Standard_Complex_Vectors.Vector ) is
  begin
    for i in labels'range loop
      put(file," ");
      Symbol_Table_io.put(file,Symbol_Table.Get(natural32(labels(i))));
      put(file," : ");
      put(file,values(i));
      new_line(file);
    end loop;
  end Write_Complex_Parameter_Values;

  procedure Write_Real_Parameter_Values
             ( file : in file_type;
               labels : in Standard_Integer_Vectors.Vector;
               values : in Standard_Complex_Vectors.Vector ) is
  begin
    for i in labels'range loop
      put(file," ");
      Symbol_Table_io.put(file,Symbol_Table.Get(natural32(labels(i))));
      put(file," : ");
      put(file,Standard_Complex_Numbers.REAL_PART(values(i)));
      new_line(file);
    end loop;
  end Write_Real_Parameter_Values;

  procedure Determine_Parameter_Values
              ( file : in file_type;
                sols : in Standard_Complex_Solutions.Solution_List;
                par : in Standard_Integer_Vectors.Vector;
                isreal : in out boolean;
                start,target : out Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

  begin
    start := Define_Start(Head_Of(sols).v,par);
    if isreal then
      for i in start'range loop
        isreal := IMAG_PART(start(i)) = 0.0;
        exit when not isreal;
      end loop;
    end if;
    if isreal then
      put_line("Real start values for the parameters :");
      Write_Real_Parameter_Values(Standard_Output,par,start);
      new_line(file);
      put_line(file,"REAL START VALUES FOR THE PARAMETERS :");
      Write_Real_Parameter_Values(file,par,start);
      target := Define_Real_Target(par);
      put_line("Real target values for the parameters :");
      Write_Real_Parameter_Values(Standard_Output,par,target);
      new_line(file);
      put_line(file,"REAL TARGET VALUES FOR THE PARAMETERS :");
      Write_Real_Parameter_Values(file,par,target);
    else
      put_line("Complex start values for the parameters :");
      Write_Complex_Parameter_Values(Standard_Output,par,start);
      new_line(file);
      put_line(file,"Complex START VALUES FOR THE PARAMETERS :");
      Write_Complex_Parameter_Values(file,par,start);
      target := Define_Complex_Target(par);
      put_line("Complex target values for the parameters :");
      Write_Complex_Parameter_Values(Standard_Output,par,target);
      new_line(file);
      put_line(file,"COMPLEX TARGET VALUES FOR THE PARAMETERS :");
      Write_Complex_Parameter_Values(file,par,target);
    end if;
  end Determine_Parameter_Values;

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

  function Circulate ( a,b : Standard_Complex_Vectors.Vector;
                       gamma,t : Standard_Complex_Numbers.Complex_Number )
                     return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Numbers;

    s : constant Complex_Number := t + gamma*t*(Create(1.0)-t);

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

  procedure Standard_Parameter_Continuation
              ( file : in file_type;
                n : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                pars : in Standard_Integer_Vectors.Vector;
                vars : in Standard_Integer_Vectors.Vector;
		sols : in out Standard_Complex_Solutions.Solution_List;
                output : in boolean ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Complex_Solutions;

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

--    procedure Simulate_Continuation is
--
--    -- DESCRIPTION :
--    --   Simulates ten steps in the continuation.
--
--      t : Complex_Number := Create(0.0);
--      p : Standard_Complex_Vectors.Vector(pars'range);
--      s : Standard_Complex_Vectors.Vector(vars'range)
--        := Head_Of(sols).v;
--    begin
--      for i in 0..10 loop
--        put("Parameters at t = ");
--        put(t); put_line(" :");
--        p := Evaluate_Parameters(t);
--        put_line(p);
--        put_line("value at current solution : ");
--        put_line(Eval(s,t));
--        put_line("derivatives at current solution : ");
--        put(Diff(s,t));
--        t := t + 0.1;
--      end loop;
--    end Simulate_Continuation;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Eval,Diff_t,Diff);
    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,Eval,Diff_t,Diff);

  begin
   -- Simulate_Continuation;
    if output
     then Rep_Cont(file,sols,false,Create(1.0));
     else Sil_Cont(sols,false,Create(1.0));
    end if;
    Clear(ep); Clear(jm); Clear(ejm);
  end Standard_Parameter_Continuation;

  procedure Call_Standard_Root_Refiner
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List ) is

    epsxa,epsfa,tolsing : double_float;
    numit : natural32 := 0;
    max : constant natural32 := 5;
    deflate : boolean := false;

  begin
    epsxa := 1.0E-14;
    epsfa := 1.0E-14;
    tolsing := 1.0E-08;
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,max,deflate,false);
  end Call_Standard_Root_Refiner;

  function Select_Symbols ( s : Array_of_Symbols;
                            v : Standard_Integer_Vectors.Vector ) 
                          return Array_of_Symbols is

    res : Array_of_Symbols(v'range);

  begin
    for i in v'range loop
      res(i) := s(integer32(v(i)));
    end loop;
    return res;
  end Select_Symbols;

  procedure Coefficient_Parameter_Homotopy_Continuation
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                nb_equ,nb_unk,nb_par : in integer32 ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    ind_par : constant Standard_Integer_Vectors.Vector
            := Define_Parameters(nb_equ,nb_unk,nb_par);
   -- nb_var : constant integer32 := nb_unk - nb_par;
    ind_var : constant Standard_Integer_Vectors.Vector
            := Complement(nb_unk,ind_par);
    start_pars : Standard_Complex_Vectors.Vector(ind_par'range);
    target_pars : Standard_Complex_Vectors.Vector(ind_par'range);
    var_sols,new_sols : Solution_List;
    gamma : constant Complex_Number := Random1;
    oc : natural32;
    output : boolean;
    s : constant Array_of_Symbols := Symbol_Table.Content;
    sv : constant Array_of_Symbols := Select_Symbols(s,ind_var);
    sp : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    isreal : boolean := Standard_Complex_to_Real_Poly.Is_Real(p);

    function Eval_Pars ( t : Complex_Number )
                       return Standard_Complex_Vectors.Vector is
    begin
     -- return Interpolate(start_pars,target_pars,t);
      return Circulate(start_pars,target_pars,gamma,t);
    end Eval_Pars;

    function Diff_Pars ( t : Complex_Number )
                       return Standard_Complex_Vectors.Vector is
    begin
      return Differentiate(start_pars,target_pars);
    end Diff_Pars;

    procedure Par_Con is
      new Standard_Parameter_Continuation(Eval_Pars,Diff_Pars);

  begin
    Determine_Parameter_Values
      (file,sols,ind_par,isreal,start_pars,target_pars);
    var_sols := Select_Variables(sols,nb_equ,ind_var);
    Symbol_Table.Clear;
    Symbol_Table.Init(sv);
    new_line(file);
    put_line(file,"THE SOLUTIONS : (after selecting variables)");
    put(file,Length_Of(var_sols),natural32(Head_Of(var_sols).n),var_sols);
    new_line;
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    output := (oc > 0);
    new_line;
    put_line("Starting parameter continuation...");
    put_line("See the output file for results.");
    new_line;
    Par_Con(file,nb_unk,p,ind_par,ind_var,var_sols,output);
    sp := Substitute(p,ind_par,target_pars);
    new_line(file);
    put_line(file,"THE TARGET SYSTEM :");
    put(file,natural32(sp'last),sp);
    Call_Standard_Root_Refiner(file,sp,var_sols);
    new_sols := Join_Variables(var_sols,nb_unk,ind_var,ind_par,target_pars);
    Symbol_Table.Clear;
    Symbol_Table.Init(s);
    new_line(file);
    put_line(file,"THE SOLUTIONS : (with target values of parameters)");
    put(file,Length_Of(new_sols),natural32(Head_Of(new_sols).n),new_sols);
  end Coefficient_Parameter_Homotopy_Continuation;

  function Create_Term
             ( n,k : natural32 ) return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    t : Term;
    p : Poly;

  begin
    t.cf := Create(1.0);
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.dg(integer32(k)) := 1;
    p := Create(t);
    return p;
  end Create_Term;

  function Complex_Sweep_Line
             ( n,k : integer32;
               start,target : Standard_Complex_Numbers.Complex_Number )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    t : Term;
    res : Poly;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n+1 => 0);
    t.cf := -start;
    res := Create(t);
    t.cf := Create(1.0);
    t.dg(k) := 1;
    Add(res,t);
    t.dg(k) := 0;
    t.dg(n+1) := 1;
    t.cf := start-target;
    Add(res,t);
    Clear(t);
    return res;
  end Complex_Sweep_Line;

  procedure Add_Symbol_for_Continuation_Parameter is

  -- DESCRIPTION :
  --   Adds the symbol 't' to the symbol table.

    sb : Symbol;
 
  begin
    sb := (sb'range => ' ');
    sb(1) := 't';
    Symbol_Table.Enlarge(1);
    Symbol_Table.Add(sb);
  end Add_Symbol_for_Continuation_Parameter;

  procedure Run_Complex_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Standard_Complex_Solutions.Link_to_Solution ) is

    ep : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(h'range)
       := Standard_Complex_Poly_SysFun.Create(h);
    jm : Standard_Complex_Jaco_Matrices.Jaco_Mat(h'range,1..integer32(nv))
       := Standard_Complex_Jaco_Matrices.Create(h);
    ejm : Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat
            (h'range,1..integer32(nv))
        := Standard_Complex_Jaco_Matrices.Create(jm);

  begin
    Run_Complex_Sweep(file,output,k,nq,nv,h,ep,ejm,s);
    Standard_Complex_Poly_SysFun.Clear(ep);
    Standard_Complex_Jaco_Matrices.Clear(jm);
    Standard_Complex_Jaco_Matrices.Clear(ejm);
  end Run_Complex_Sweep;

  procedure Run_Real_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in Standard_Floating_Poly_Systems.Poly_Sys;
                s : in Standard_Complex_Solutions.Link_to_Solution ) is

    ep : Standard_Floating_Poly_SysFun.Eval_Poly_Sys(h'range)
       := Standard_Floating_Poly_SysFun.Create(h);
    jm : Standard_Floating_Jaco_Matrices.Jaco_Mat(h'range,1..integer32(nv))
       := Standard_Floating_Jaco_Matrices.Create(h);
    ejm : Standard_Floating_Jaco_Matrices.Eval_Jaco_Mat
            (h'range,1..integer32(nv))
        := Standard_Floating_Jaco_Matrices.Create(jm);

  begin
    Run_Real_Sweep(file,output,k,nq,nv,h,ep,ejm,s);
    Standard_Floating_Poly_SysFun.Clear(ep);
    Standard_Floating_Jaco_Matrices.Clear(jm);
    Standard_Floating_Jaco_Matrices.Clear(ejm);
  end Run_Real_Sweep;

  function Hyperplane ( dx,x : Standard_Complex_Vectors.Vector )
                      return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;
    t : Term;
    y : Complex_Number := Create(0.0);

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..dx'last => 0);
    for i in 1..dx'last loop
      if dx(i) /= Create(0.0) then
        t.dg(i) := 1;
        t.cf := dx(i);
        Add(res,t);
        t.dg(i) := 0;
        y := x(i)*dx(i);
      end if;
    end loop;
    t.cf := -y;
    Add(res,t);
    Clear(t);
    return res;
  end Hyperplane;

  function Hyperplane ( dx,x : Standard_Floating_Vectors.Vector )
                      return Standard_Floating_Polynomials.Poly is

    use Standard_Complex_Numbers;
    use Standard_Floating_Polynomials;

    res : Poly := Null_Poly;
    t : Term;
    y : double_float := 0.0;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..dx'last => 0);
    for i in 1..dx'last loop
      if dx(i) /= 0.0 then
        t.dg(i) := 1;
        t.cf := dx(i);
        Add(res,t);
        t.dg(i) := 0;
        y := x(i)*dx(i);
      end if;
    end loop;
    t.cf := -y;
    Add(res,t);
    Clear(t);
    return res;
  end Hyperplane;

  procedure Run_Complex_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                f : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                s : in Standard_Complex_Solutions.Link_to_Solution ) is

    use Standard_Complex_Numbers;

    x,dx : Standard_Complex_Vectors.Vector(1..integer32(nv));
   -- st : Standard_Complex_Poly_Systems.Poly_Sys(1..h'last+1);
   -- sols : Solution_List;
   -- nb : natural := 0;
   -- ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    for i in s.v'range loop
      x(i) := s.v(i);
    end loop;
    x(integer32(nv)) := Create(0.0);
    new_line(file);
    put(file,"Starting the complex sweep at solution ");
    put(file,k,1); put_line(file," :");
    Write_Vector(file,x);
    new_line(file);
   -- if output
   --  then
    Start_Complex_Sweep(file,output,nq,nv,1.0,f,jf,x,dx);
   --  else Start_Complex_Sweep(nq,nv,1.0,f,jf,x,dx);
   -- end if;
    put_line(file,"The solution and its tangent at the end");
    Write_Vector_and_its_Tangent(file,x,dx);
    s.t := x(x'last);
    for i in s.v'range loop
      s.v(i) := x(i);
    end loop;
   -- st(h'range) := h;
   -- st(st'last) := Hyperplane(dx,x);
   -- put_line(file,"The system used to refine : "); put(file,st);
   -- ls := new Solution(nv);
   -- ls.t := s.t;
   -- ls.m := s.m;
   -- ls.v := x;
   -- Construct(ls,sols);
   -- Reporting_Root_Refiner
   --   (file,st,sols,1.0E-12,1.0E-12,1.0E-8,nb,5,true,false);
   -- Standard_Complex_Polynomials.Clear(st(st'last));
   -- ls := Head_Of(sols);
   -- s.v := ls.v(s.v'range);
   -- s.t := ls.v(ls.v'last);
  end Run_Complex_Sweep;

  procedure Run_Real_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in Standard_Floating_Poly_Systems.Poly_Sys;
                f : in Standard_Floating_Poly_SysFun.Eval_Poly_Sys;
                jf : in Standard_Floating_Jaco_Matrices.Eval_Jaco_Mat;
                s : in Standard_Complex_Solutions.Link_to_Solution ) is

    use Standard_Complex_Numbers;

    x,dx : Standard_Floating_Vectors.Vector(1..integer32(nv));
    evl : constant boolean := false;
   -- st : Standard_Floating_Poly_Systems.Poly_Sys(1..h'last+1);
   -- cst : Standard_Complex_Poly_Systems.Poly_Sys(1..h'last+1);
   -- sols : Solution_List;
   -- nb : natural := 0;
   -- ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    for i in s.v'range loop
      x(i) := REAL_PART(s.v(i));
    end loop;
    x(integer32(nv)) := 0.0;
    new_line(file);
    put(file,"Starting the real sweep at solution ");
    put(file,k,1); put_line(file," :");
    Write_Vector(file,x);
    new_line(file);
   -- if output
   --  then
    Start_Real_Sweep(file,output,evl,nq,nv,1.0,f,jf,x,dx);
   --  else Start_Real_Sweep(evl,nq,nv,1.0,f,jf,x,dx);
   -- end if;
    put_line(file,"The solution and its tangent at the end");
    Write_Vector_and_its_Tangent(file,x,dx);
    s.t := Create(x(x'last));
    for i in s.v'range loop
      s.v(i) := Create(x(i));
    end loop;
   -- st(h'range) := h;
   -- st(st'last) := Hyperplane(dx,x);   
   -- put_line(file,"the system used to refine : "); put(file,st);
   -- cst := Standard_Complex_to_Real_Poly.Convert_Real_to_Complex(st);
   -- ls := new Solution(nv);
   -- ls.t := s.t;
   -- ls.m := s.m;
   -- for i in x'range loop
   --   ls.v(i) := Create(x(i));
   -- end loop;
   -- Construct(ls,sols);
   -- Reporting_Root_Refiner
   --   (file,cst,sols,1.0E-12,1.0E-12,1.0E-8,nb,5,true,false);
   -- Standard_Floating_Polynomials.Clear(st(st'last));
   -- Standard_Complex_Poly_Systems.Clear(cst);
   -- ls := Head_Of(sols);
   -- s.v := ls.v(s.v'range);
   -- s.t := ls.v(ls.v'last);
  end Run_Real_Sweep;

  procedure Sweep ( file : in file_type; isreal : in out boolean;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    sols : in Standard_Complex_Solutions.Solution_List;
                    nb_equ,nb_unk,nb_par : in integer32 ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    par : Standard_Integer_Vectors.Vector(1..nb_par);
   -- nb_var : constant natural := nb_unk-nb_par;
   -- var : Standard_Integer_Vectors.Vector(1..nb_var);
    start : Standard_Complex_Vectors.Vector(par'range);
    target : Standard_Complex_Vectors.Vector(par'range);
    n : constant natural32 := natural32(nb_unk);
    spf : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(p'range);
    sjm : Standard_Complex_Jaco_Matrices.Jaco_Mat(p'range,1..integer32(n));
    sjf : Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat
            (p'range,1..integer32(n));
    rp : Standard_Floating_Poly_Systems.Poly_Sys(p'range);
    rpf : Standard_Floating_Poly_SysFun.Eval_Poly_Sys(p'range);
    rjm : Standard_Floating_Jaco_Matrices.Jaco_Mat(p'range,1..integer32(n));
    rjf : Standard_Floating_Jaco_Matrices.Eval_Jaco_Mat
            (p'range,1..integer32(n));
    tmp : Solution_List;
    ls : Link_to_Solution;
    len : constant natural32 := Length_Of(sols);
    ans : character;
    output : boolean;

  begin
    par := Define_Parameters(nb_equ,nb_unk,nb_par);
   -- var := Complement(nb_unk,par);
    Determine_Parameter_Values(file,sols,par,isreal,start,target);
    new_line;
    put("Starting sweep at "); put(len,1); put_line(" solutions.");
    new_line;
    put("Do you want intermediate output along the paths ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    new_line;
    Standard_Quad_Parameters.Tune;
    new_line;
    put_line("See the output file for results...");
    new_line;
    tstart(timer);
   -- Add_Symbol_for_Continuation_Parameter;
   -- for i in p'range loop
   --   sp(i) := Standard_Embed_Polynomials.Add_Variables(p(i),1); 
   -- end loop;
   -- for i in 1..nb_par loop
   --   sp(p'last+i) := Complex_Sweep_Line(nb_unk,par(i),start(i),target(i));
   -- end loop;
    spf := Standard_Complex_Poly_SysFun.Create(p);
    sjm := Standard_Complex_Jaco_Matrices.Create(p);
    sjf := Standard_Complex_Jaco_Matrices.Create(sjm);
    if isreal then
      rp := Standard_Complex_to_Real_Poly.Convert_Complex_to_Real(p);
      rpf := Standard_Floating_Poly_SysFun.Create(rp);
      rjm := Standard_Floating_Jaco_Matrices.Create(rp);
      rjf := Standard_Floating_Jaco_Matrices.Create(rjm);
    end if;
    new_line(file);
    put_line(file,"THE SWEEP HOMOTOPY :");
    put(file,natural32(p'last),n,p);
    tmp := sols;
    for i in 1..len loop
      ls := Head_Of(tmp);
      ls.t := Create(0.0);
      if Standard_Solution_Diagnostics.Is_Real(ls.all,1.0E-14) then
        if isreal then
          Run_Real_Sweep(file,output,i,natural32(rp'last),n,rp,rpf,rjf,ls);
        else
          Run_Complex_Sweep(file,output,i,natural32(p'last),n,p,spf,sjf,ls);
        end if;
      else
        Run_Complex_Sweep(file,output,i,natural32(p'last),n,p,spf,sjf,ls);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    if isreal then
      Standard_Floating_Poly_Systems.Clear(rp);
      Standard_Floating_Poly_SysFun.Clear(rpf);
      Standard_Floating_Jaco_Matrices.Clear(rjm);
      Standard_Floating_Jaco_Matrices.Clear(rjf);
    end if;
    tstop(timer);
    new_line(file);
    put_line(file,"THE SOLUTIONS at the end of the sweep :");
    put(file,len,natural32(nb_unk),sols);
    Write_Sweep_Summary(file,sols);
    new_line(file);
    print_times(file,timer,"running a sweep");
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Sweep;

  function Show_Menu return character is

    ans : character;

  begin
    put_line("MENU for coefficient parameter polynomial continuation :");
    put_line("  1. complex homotopy with random gamma constant;");
    put_line("  2. sweep in parameter space with singularity detection.");
    put("Type 1 or 2 to choose : "); Ask_Alternative(ans,"12");
    return ans;
  end Show_Menu;

end Parameter_Homotopy_Continuation;
