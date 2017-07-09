with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Multprec_Complex_Matrices;          use Multprec_Complex_Matrices;
with Multprec_Complex_Linear_Solvers;    use Multprec_Complex_Linear_Solvers;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;
with Multprec_Solution_Diagnostics;      use Multprec_Solution_Diagnostics;
with Multprec_Condition_Tables;          use Multprec_Condition_Tables;
with Multprec_Complex_Newton_Steps;      use Multprec_Complex_Newton_Steps;

package body Multprec_Root_Refiners is

-- AUXILIARIES :

  procedure Write_Bar ( file : in file_type ) is
  begin
    for i in 1..75 loop
      put(file,'=');
    end loop;
    new_line(file);
  end Write_Bar;

  procedure Write_Info ( file : in file_type; zero : in Solution;
                         initres : double_float;
                         i,numb : in natural32; fail : in boolean ) is

  -- DESCRIPTION :
  --   The information concerning the zero is written

  begin
    put(file,"solution "); put(file,i,1); put(file," : ");
    put(file,"   start residual : "); put(file,initres,2,3,3);
    put(file,"   #iterations : "); put(file,numb,1);
    if fail
     then put_line(file,"   failure");
     else put_line(file,"   success");
    end if;
    put(file,zero);
  end Write_Info;
 
  procedure Root_Accounting
               ( file : in file_type; ls : in out Link_to_Solution;
                 nb : in natural32; sa : in out Solution_Array;
                 fail : in boolean; tolsing,tolclus : in Floating_Number;
                 nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus 
                   : in out natural32 ) is

  -- DESCRIPTION :
  --   This procedure does root accounting of the solution sol, w.r.t. the
  --   solution list sols.  Information will be provided concerning the type
  --   of solution.

    tolreal : Floating_Number := Create(1.0E-13);

  begin
    if fail then
      put_line(file," no solution ==");
      nbfail := nbfail + 1;
      ls.m := 0;
    else 
      if Is_Real(ls.all,tolreal) then
        put(file," real ");
        nbreal := nbreal + 1;
      else
        put(file," complex ");
        nbcomp := nbcomp + 1;
      end if;
      if sa(integer32(nb)).rco < tolsing then
        declare
          m : natural32 := Multiplicity(ls.all,sa,tolclus);
        begin
          if m = 1
           then m := 0;
          end if;
          ls.m := integer32(m);
        end;
        put_line(file,"singular ==");
        nbsing := nbsing + 1;
      else
        declare
          nb2 : constant natural32 := Is_Clustered(ls.all,nb,sa,tolclus);
        begin
          if nb2 = nb then
            put_line(file,"regular ==");
            nbreg := nbreg + 1;
            -- ls.m := 1;
          else 
            put(file,"clustered : ");
            put(file,nb2,1);
            put_line(file," ==");
            nbclus := nbclus + 1;
            ls.m := -integer32(nb2);
            sa(integer32(nb2)).m := -integer32(nb);
          end if;
        end;	   
      end if;  
    end if;
    Clear(tolreal);
  end Root_Accounting;

  procedure Root_Accounting 
                ( ls : in out Link_to_Solution; nb : in natural32;
                  sa : in out Solution_Array; fail : in boolean;
                  tolsing,tolclus : in Floating_Number ) is

  -- DESCRIPTION :
  --   This procedure does root accounting of the solution sol, w.r.t. the
  --   solution list sols.  Information will be provided concerning the type
  --   of solution.

  begin
    if fail then
      ls.m := 0;
    elsif sa(integer32(nb)).rco < tolsing then
      declare
        m : constant natural32 := Multiplicity(ls.all,sa,tolclus);
      begin
        if m = 1
         then ls.m := 0;
         else ls.m := integer32(m);
        end if;
      end;
    else -- sa(nb).rco > tolsing, check for clustering
      declare
        nb2 : constant natural32 := Is_Clustered(ls.all,nb,sa,tolclus);
      begin
        if nb2 /= nb
         then ls.m := -integer32(nb2);
              sa(integer32(nb2)).m := -integer32(nb);
        end if;
       -- if nb2 = nb
       --  then ls.m := 1;
       --  else ls.m := -nb2;
       -- end if;
      end;	   
    end if;
  end Root_Accounting;

  procedure Write_Global_Info
             ( file : in file_type;
               tot,nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus : in natural32 ) is

  begin
    Write_Bar(file);
    put(file,"A list of "); put(file,tot,1);
    put_line(file," solutions has been refined :");
    put(file,"Number of regular solutions   : "); put(file,nbreg,1);
    put_line(file,".");
    put(file,"Number of singular solutions  : "); put(file,nbsing,1);
    put_line(file,".");
    put(file,"Number of real solutions      : "); put(file,nbreal,1);
    put_line(file,".");
    put(file,"Number of complex solutions   : "); put(file,nbcomp,1);
    put_line(file,".");
    put(file,"Number of clustered solutions : "); put(file,nbclus,1);
    put_line(file,".");
    put(file,"Number of failures            : "); put(file,nbfail,1);
    put_line(file,".");
    Write_Bar(file);
  end Write_Global_Info;

-- ONE NEWTON STEP :

  procedure Multprec_Newton_Step
              ( f : in Multprec_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out Multprec_Complex_Vectors.Vector;
                err,rco,res : out Floating_Number ) is

    y : Vector(f'range) := eval(f,x);
    A : Matrix(f'range,x'range) := eval(jf,x);
    ipvt : Standard_Integer_Vectors.Vector(A'range(1));

  begin
    Min(y);
    Clear(rco);
    lufco(A,A'last(1),ipvt,rco);
    lusolve(A,A'last(1),ipvt,y);
    Add(x,y);
    Clear(err); err := Max_Norm(y);
    Clear(y);
    y := eval(f,x);
    Clear(res); res := Max_Norm(y);
    Clear(y); Clear(A);
  end Multprec_Newton_Step;

  procedure Multprec_Newton_Step
              ( f : in Multprec_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in Multprec_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out Multprec_Complex_Vectors.Vector;
                err,rco,res : out Floating_Number ) is

    y : Vector(f'range) := Multprec_Complex_Laur_SysFun.Eval(f,x);
    A : Matrix(f'range,x'range)
      := Multprec_Complex_Laur_JacoMats.Eval(jf,x);
    ipvt : Standard_Integer_Vectors.Vector(A'range(1));

  begin
    Min(y);
    Clear(rco);
    lufco(A,A'last(1),ipvt,rco);
    lusolve(A,A'last(1),ipvt,y);
    Add(x,y);
    Clear(err); err := Max_Norm(y);
    Clear(y);
    y := Multprec_Complex_Laur_SysFun.Eval(f,x);
    Clear(res); res := Max_Norm(y);
    Clear(y); Clear(A);
  end Multprec_Newton_Step;

-- TARGET ROUTINES, NEWTON's METHOD :

  procedure Silent_Newton
               ( p_eval : in Eval_Poly_Sys;
                 j_eval : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                 zero : in out Solution; epsxa,epsfa : in Floating_Number; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    n : constant integer32 := p_eval'length;
    jacobian : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    y,deltax : Vector(1..n);
    tmp : Floating_Number;

  begin
    y := Eval(p_eval,zero.v);               -- y = f(zero)
    for i in 1..max loop
      jacobian := Eval(j_eval,zero.v);      -- solve jacobian*deltax = -f(zero)
     -- lufco(jacobian,n,ipvt,zero.rco);
      lufac(jacobian,n,ipvt,info);
     -- if Equal(zero.rco,0.0)
      if info /= 0 then
        tmp := Sum_Norm(y);
        fail := (tmp > epsfa);
        Clear(tmp);
        return;                             -- singular Jacobian matrix
      end if;
      deltax := -y; Clear(y);
      lusolve(jacobian,n,ipvt,deltax);  
      Add(zero.v,deltax);                   -- make the updates
      y := Eval(p_eval,zero.v);
      Clear(zero.err);              Clear(zero.res);
      zero.err := Sum_Norm(deltax); zero.res := Sum_Norm(y);
      Clear(jacobian); Clear(deltax); 
     -- Clear(y);
      numit := numit + 1;
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) then
        fail := false; exit; -- stopping criteria
      elsif numit >= max then
        fail := true; exit;
      end if;
    end loop;
    jacobian := Eval(j_eval,zero.v);        -- compute condition number
    lufco(jacobian,n,ipvt,zero.rco);
    Clear(jacobian);
  exception
    when constraint_error => fail := true; return;
  end Silent_Newton;

  procedure Silent_Newton
               ( p_eval : in Eval_Laur_Sys;
                 j_eval : in Multprec_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                 zero : in out Solution; epsxa,epsfa : in Floating_Number; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    n : constant integer32 := p_eval'length;
    jacobian : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    y,deltax : Vector(1..n);
    tmp : Floating_Number;

    use Multprec_Complex_Laur_JacoMats;

  begin
    y := Eval(p_eval,zero.v);               -- y = f(zero)
    for i in 1..max loop
      jacobian := Eval(j_eval,zero.v); -- solve jacobian*deltax = -f(zero)
     -- lufco(jacobian,n,ipvt,zero.rco);
      lufac(jacobian,n,ipvt,info);
     -- if Equal(zero.rco,0.0)
      if info /= 0 then
        tmp := Sum_Norm(y);
        fail := (tmp > epsfa);
        Clear(tmp);
        return;                             -- singular Jacobian matrix
      end if;
      deltax := -y; Clear(y);
      lusolve(jacobian,n,ipvt,deltax);  
      Add(zero.v,deltax);                   -- make the updates
      y := Eval(p_eval,zero.v);
      Clear(zero.err);              Clear(zero.res);
      zero.err := Sum_Norm(deltax); zero.res := Sum_Norm(y);
      Clear(jacobian); Clear(deltax); 
     -- Clear(y);
      numit := numit + 1;
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) then
        fail := false; exit; -- stopping criteria
      elsif numit >= max then
        fail := true; exit;
      end if;
    end loop;
    jacobian := Eval(j_eval,zero.v);        -- compute condition number
    lufco(jacobian,n,ipvt,zero.rco);
    Clear(jacobian);
  exception
    when constraint_error => fail := true; return;
  end Silent_Newton;

  procedure Reporting_Newton
               ( file : in file_type; p_eval : in Eval_Poly_Sys;
                 j_eval : in Multprec_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                 zero : in out Solution; epsxa,epsfa : in Floating_Number; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    n : constant integer32 := p_eval'length;
    jacobian : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    y,deltax : Vector(1..n);
    tmp : Floating_Number;

  begin
    y := Eval(p_eval,zero.v);              -- y = f(zero)
    for i in 1..max loop
      jacobian := Eval(j_eval,zero.v);     -- solve jacobian*deltax = -f(zero)
      lufac(jacobian,n,ipvt,info);
      if info /= 0 then
        tmp := Sum_Norm(y);
        fail := (tmp > epsfa);             -- accuracy not reached yet
        Clear(tmp);
        Clear(jacobian); Clear(y);
        return;
      end if;
      deltax := -y; Clear(y);
      lusolve(jacobian,n,ipvt,deltax);  
      Add(zero.v,deltax);                  -- make the updates
      y := Eval(p_eval,zero.v);
      Clear(zero.err);              Clear(zero.res);
      zero.err := Sum_Norm(deltax); zero.res := Sum_Norm(y);
      Clear(deltax); Clear(jacobian);
      numit := numit + 1;
      put(file,"Step "); put(file,numit,4); new_line(file);      -- output
      put(file," |errxa| : "); put(file,zero.err); new_line(file);
      put(file," |errfa| : "); put(file,zero.res); new_line(file);
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) then
        fail := false; exit;                          -- stopping criteria
      elsif numit >= max then
        fail := true; exit;
      end if;
    end loop;
    Clear(y);
    jacobian := Eval(j_eval,zero.v);              -- compute condition number
    lufco(jacobian,n,ipvt,zero.rco);
    Clear(jacobian);
  exception
    when constraint_error => fail := true; return;
  end Reporting_Newton;

  procedure Reporting_Newton
               ( file : in file_type; p_eval : in Eval_Laur_Sys;
                 j_eval : in Multprec_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                 zero : in out Solution; epsxa,epsfa : in Floating_Number; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    n : constant integer32 := p_eval'length;
    jacobian : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    y,deltax : Vector(1..n);
    tmp : Floating_Number;

    use Multprec_Complex_Laur_JacoMats;

  begin
    y := Eval(p_eval,zero.v);              -- y = f(zero)
    for i in 1..max loop
      jacobian := Eval(j_eval,zero.v);     -- solve jacobian*deltax = -f(zero)
      lufac(jacobian,n,ipvt,info);
      if info /= 0 then
        tmp := Sum_Norm(y);
        fail := (tmp > epsfa);             -- accuracy not reached yet
        Clear(tmp);
        Clear(jacobian); Clear(y);
        return;
      end if;
      deltax := -y; Clear(y);
      lusolve(jacobian,n,ipvt,deltax);  
      Add(zero.v,deltax);                  -- make the updates
      y := Eval(p_eval,zero.v);
      Clear(zero.err);              Clear(zero.res);
      zero.err := Sum_Norm(deltax); zero.res := Sum_Norm(y);
      Clear(deltax); Clear(jacobian);
      numit := numit + 1;
      put(file,"Step "); put(file,numit,4); new_line(file);      -- output
      put(file," |errxa| : "); put(file,zero.err); new_line(file);
      put(file," |errfa| : "); put(file,zero.res); new_line(file);
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) then
        fail := false; exit;                          -- stopping criteria
      elsif numit >= max then
        fail := true; exit;
      end if;
    end loop;
    Clear(y);
    jacobian := Eval(j_eval,zero.v);              -- compute condition number
    lufco(jacobian,n,ipvt,zero.rco);
    Clear(jacobian);
  exception
    when constraint_error => fail := true; return;
  end Reporting_Newton;

-- FOR OVERDETERMINED SYSTEMS :

  procedure Silent_Gauss_Newton
               ( p_eval : in Eval_Poly_Sys; j_eval : in Eval_Jaco_Mat;
                 zero : in out Solution; tol : in double_float;
                 epsxa,epsfa : in Floating_Number; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    rank : natural32;

    function f ( x : Vector ) return Vector is
    begin
      return Eval(p_eval,x);
    end f;

    function jf ( x : Vector ) return Matrix is
    begin
      return Eval(j_eval,x);
    end jf;

    procedure Newton is new Silent_Newton_Step(f,jf);

  begin
    for i in 1..max loop
      Newton(natural32(p_eval'last),zero.v,tol,
             zero.err,zero.rco,zero.res,rank);
      numit := i;
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) then
        fail := false; exit;
      elsif numit >= max then
        fail := true; exit;
      end if;
    end loop;
  end Silent_Gauss_Newton;

  procedure Reporting_Gauss_Newton
               ( file : in file_type;
                 p_eval : in Eval_Poly_Sys; j_eval : in Eval_Jaco_Mat;
                 zero : in out Solution; tol : in double_float;
                 epsxa,epsfa : in Floating_Number; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    rank : natural32;

    function f ( x : Vector ) return Vector is
    begin
      return Eval(p_eval,x);
    end f;

    function jf ( x : Vector ) return Matrix is
    begin
      return Eval(j_eval,x);
    end jf;

    procedure Newton is new Reporting_Newton_Step(f,jf);

  begin
    for i in 1..max loop
      Newton(file,natural32(p_eval'last),zero.v,tol,
             zero.err,zero.rco,zero.res,rank);
      numit := i;
      put(file,"Step "); put(file,numit,4); new_line(file);      -- output
      put(file," |errxa| : "); put(file,zero.err); new_line(file);
      put(file," |errfa| : "); put(file,zero.res); new_line(file);
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) then
        fail := false; exit;
      elsif numit >= max then
        fail := true; exit;
      end if;
    end loop;
    null;
  end Reporting_Gauss_Newton;

-- APPLICATION OF NEWTON's METHOD TO LISTS OF SOLUTIONS :

  procedure Silent_Root_Refiner
               ( p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate : in boolean  ) is

    n : constant integer32 := p'length;
    p_eval : Eval_Poly_Sys(1..n) := Create(p);
    jac : Jaco_Mat(1..n,1..n) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..n,1..n) := Create(jac);
    numb : natural32;
    fail : boolean;
    sa : Solution_Array(1..integer32(Length_Of(sols))) := Create(sols);
    eva : Floating_Number;
    eval_acc : Vector(1..n);

  begin
    for i in sa'range loop
      numb := 0;
      eval_acc := Eval(p_eval,sa(i).v);
      eva := Sum_Norm(eval_acc);
      Clear(eval_acc);
      if eva < 1.0 then
        Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
      else
        fail := true;
      end if;
      Root_Accounting(sa(i),natural32(i),sa(sa'first..i),fail,tolsing,epsxa);
      numit := numit + numb;
      Clear(eva);
    end loop;
    Clear(p_eval); Clear(jac); Clear(jac_eval);
    Shallow_Clear(sols); sols := Create(sa); -- Shallow_Clear(sa);
  end Silent_Root_Refiner;

  procedure Silent_Root_Refiner
               ( p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate : in boolean ) is

    n : constant integer32 := p'length;
    p_eval : Eval_Poly_Sys(1..n) := Create(p);
    jac : Jaco_Mat(1..n,1..n) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..n,1..n) := Create(jac);
    numb : natural32;
    fail : boolean;
    sa : Solution_Array(1..integer32(Length_Of(sols))) := Create(sols);
    eva : Floating_Number;
    refsols_last : Solution_List;
    eval_acc : Vector(1..n);

  begin
    for i in sa'range loop
      numb := 0;
      eval_acc := Eval(p_eval,sa(i).v);
      eva := Sum_Norm(eval_acc);
      Clear(eval_acc);
      if eva < 1.0 then
        Silent_Newton
          (p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
      else
        fail := true;
      end if;
      Root_Accounting(sa(i),natural32(i),sa(sa'first..i),fail,tolsing,epsxa);
      if not fail
       then Append(refsols,refsols_last,sa(i).all);
      end if;
      numit := numit + numb;
      Clear(eva);
    end loop;
    Clear(p_eval); Clear(jac); Clear(jac_eval);
    Shallow_Clear(sols); sols := Create(sa); -- Shallow_Clear(sa);
  end Silent_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate,wout : in boolean ) is

    n : constant integer32 := p'length;
    p_eval : Eval_Poly_Sys(1..n) := Create(p);
    jac : Jaco_Mat(1..n,1..n) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..n,1..n) := Create(jac);
    numb : natural32;
    nbfail,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    nbtot : constant natural32 := Length_Of(sols);
    fail : boolean;
    sa : Solution_Array(1..integer32(nbtot)) := Create(sols);
    initres : Standard_Floating_Vectors.Vector(sa'range);
    eval_acc : Vector(1..n);
    eva : Floating_Number;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..20)
                      := Multprec_Condition_Tables.Create(20); 

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); new_line(file);
    put(file,nbtot,1); put(file," "); put(file,n,1); new_line(file);
    Write_Bar(file);
    for i in sa'range loop 
      numb := 0;
      eval_acc := Eval(p_eval,sa(i).v);
      eva := Sum_Norm(eval_acc);
      Clear(eval_acc);
      initres(i) := Trunc(eva);
      if eva < 1.0 then
        if wout then
          Reporting_Newton
            (file,p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
        else
          Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
        end if;
      else
        fail := true;
      end if;
      Write_Info(file,sa(i).all,initres(i),natural32(i),numb,fail);
      Root_Accounting(file,sa(i),natural32(i),sa(sa'first..i),fail,tolsing,
                      epsxa,nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus);
      Update_Corrector(t_err,sa(i).all);
      Update_Condition(t_rco,sa(i).all);
      Update_Residuals(t_res,sa(i).all);
      numit := numit + numb;
      Clear(eva);
    end loop;
    Write_Global_Info(file,nbtot,nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus);
    Clear(p_eval); Clear(jac); Clear(jac_eval);
    Shallow_Clear(sols); sols := Create(sa); -- Shallow_Clear(sa);
    Write_Tables(file,t_err,t_res,t_rco);
  end Reporting_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate,wout : in boolean ) is

    n : constant integer32 := p'length;
    p_eval : Eval_Poly_Sys(1..n) := Create(p);
    jac : Jaco_Mat(1..n,1..n) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..n,1..n) := Create(jac);
    numb : natural32;
    nbfail,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    nbtot : constant natural32 := Length_Of(sols);
    fail : boolean;
    sa : Solution_Array(1..integer32(nbtot)) := Create(sols);
    initres : Standard_Floating_Vectors.Vector(sa'range);
    refsols_last : Solution_List;
    eval_acc : Vector(1..n);
    eva : Floating_Number;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..20)
                      := Multprec_Condition_Tables.Create(20); 

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); new_line(file);
    put(file,nbtot,1); put(file," "); put(file,n,1); new_line(file);
    Write_Bar(file);
    for i in sa'range loop 
      numb := 0;
      eval_acc := Eval(p_eval,sa(i).v);
      eva := Sum_Norm(eval_acc);
      Clear(eval_acc);
      initres(i) := Trunc(eva);
      if eva < 1.0 then
        if wout then
          Reporting_Newton(file,p_eval,jac_eval,sa(i).all,epsxa,epsfa,
                           numb,max,fail);
        else
          Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
        end if;
      else
        fail := true;
      end if;
      Write_Info(file,sa(i).all,initres(i),natural32(i),numb,fail);
      Root_Accounting(file,sa(i),natural32(i),sa(sa'first..i),fail,tolsing,
                      epsxa,nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus);
      if not fail
       then Append(refsols,refsols_last,sa(i).all);
      end if;
      Update_Corrector(t_err,sa(i).all);
      Update_Condition(t_rco,sa(i).all);
      Update_Residuals(t_res,sa(i).all);
      numit := numit + numb;
      Clear(eva);
    end loop;
    Write_Global_Info(file,nbtot,nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus);
    Clear(p_eval); Clear(jac); Clear(jac_eval);
    Shallow_Clear(sols); sols := Create(sa); -- Shallow_Clear(sa);
    Write_Tables(file,t_err,t_res,t_rco);
  end Reporting_Root_Refiner;

-- APPLICATION OF GAUSS-NEWTON's METHOD TO LISTS OF SOLUTIONS :

  procedure Silent_Root_Sharpener
               ( p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate : in boolean  ) is

    nq : constant integer32 := p'last;
    nv : constant integer32 := Head_Of(sols).n;
    p_eval : Eval_Poly_Sys(1..nq) := Create(p);
    jac : Jaco_Mat(1..nq,1..nv) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..nq,1..nv) := Create(jac);
    numb : natural32;
    tolrnk : constant double_float := 1.0E-6;
    fail : boolean;
    sa : Solution_Array(1..integer32(Length_Of(sols))) := Create(sols);
    eva : Floating_Number;
    eval_acc : Vector(1..nq);

  begin
    for i in sa'range loop
      numb := 0;
      eval_acc := Eval(p_eval,sa(i).v);
      eva := Sum_Norm(eval_acc);
      Clear(eval_acc);
      if eva < 1.0 then
        Silent_Gauss_Newton
          (p_eval,jac_eval,sa(i).all,tolrnk,epsxa,epsfa,numb,max,fail);
      else
        fail := true;
      end if;
      Root_Accounting(sa(i),natural32(i),sa(sa'first..i),fail,tolsing,epsxa);
      numit := numit + numb;
      Clear(eva);
    end loop;
    Clear(p_eval); Clear(jac); Clear(jac_eval);
    Shallow_Clear(sols); sols := Create(sa); -- Shallow_Clear(sa);
  end Silent_Root_Sharpener;

  procedure Silent_Root_Sharpener
               ( p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate : in boolean ) is

    nq : constant integer32 := p'last;
    nv : constant integer32 := Head_Of(sols).n;
    p_eval : Eval_Poly_Sys(1..nq) := Create(p);
    jac : Jaco_Mat(1..nq,1..nv) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..nq,1..nv) := Create(jac);
    numb : natural32;
    tolrnk : constant double_float := 1.0E-6;
    fail : boolean;
    sa : Solution_Array(1..integer32(Length_Of(sols))) := Create(sols);
    eva : Floating_Number;
    refsols_last : Solution_List;
    eval_acc : Vector(1..nq);

  begin
    for i in sa'range loop
      numb := 0;
      eval_acc := Eval(p_eval,sa(i).v);
      eva := Sum_Norm(eval_acc);
      Clear(eval_acc);
      if eva < 1.0 then
        Silent_Gauss_Newton
          (p_eval,jac_eval,sa(i).all,tolrnk,epsxa,epsfa,numb,max,fail);
      else
        fail := true;
      end if;
      Root_Accounting(sa(i),natural32(i),sa(sa'first..i),fail,tolsing,epsxa);
      if not fail
       then Append(refsols,refsols_last,sa(i).all);
      end if;
      numit := numit + numb;
      Clear(eva);
    end loop;
    Clear(p_eval); Clear(jac); Clear(jac_eval);
    Shallow_Clear(sols); sols := Create(sa); -- Shallow_Clear(sa);
  end Silent_Root_Sharpener;

  procedure Reporting_Root_Sharpener
               ( file : in file_type;
                 p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate,wout : in boolean ) is

    nq : constant integer32 := p'last;
    nv : constant integer32 := Head_Of(sols).n;
    p_eval : Eval_Poly_Sys(1..nq) := Create(p);
    jac : Jaco_Mat(1..nq,1..nv) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..nq,1..nv) := Create(jac);
    numb : natural32;
    tolrnk : constant double_float := 1.0E-6;
    nbfail,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    nbtot : constant natural32 := Length_Of(sols);
    fail : boolean;
    sa : Solution_Array(1..integer32(nbtot)) := Create(sols);
    initres : Standard_Floating_Vectors.Vector(sa'range);
    eval_acc : Vector(1..nq);
    eva : Floating_Number;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..20)
                      := Multprec_Condition_Tables.Create(20); 

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); new_line(file);
    put(file,nbtot,1); put(file," "); put(file,nv,1); new_line(file);
    Write_Bar(file);
    for i in sa'range loop 
      numb := 0;
      eval_acc := Eval(p_eval,sa(i).v);
      eva := Sum_Norm(eval_acc);
      Clear(eval_acc);
      initres(i) := Trunc(eva);
      if eva < 1.0 then
        if wout then
          Reporting_Gauss_Newton
            (file,p_eval,jac_eval,sa(i).all,tolrnk,epsxa,epsfa,numb,max,fail);
        else
          Silent_Gauss_Newton
            (p_eval,jac_eval,sa(i).all,tolrnk,epsxa,epsfa,numb,max,fail);
        end if;
      else
        fail := true;
      end if;
      Write_Info(file,sa(i).all,initres(i),natural32(i),numb,fail);
      Root_Accounting(file,sa(i),natural32(i),sa(sa'first..i),fail,tolsing,
                      epsxa,nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus);
      Update_Corrector(t_err,sa(i).all);
      Update_Condition(t_rco,sa(i).all);
      Update_Residuals(t_res,sa(i).all);
      numit := numit + numb;
      Clear(eva);
    end loop;
    Write_Global_Info(file,nbtot,nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus);
    Clear(p_eval); Clear(jac); Clear(jac_eval);
    Shallow_Clear(sols); sols := Create(sa); -- Shallow_Clear(sa);
    Write_Tables(file,t_err,t_res,t_rco);
  end Reporting_Root_Sharpener;

  procedure Reporting_Root_Sharpener
               ( file : in file_type;
                 p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in Floating_Number;
                 numit : in out natural32; max : in natural32;
                 deflate,wout : in boolean ) is

    nq : constant integer32 := p'last;
    nv : constant integer32 := Head_Of(sols).n;
    p_eval : Eval_Poly_Sys(1..nq) := Create(p);
    jac : Jaco_Mat(1..nq,1..nv) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..nq,1..nv) := Create(jac);
    numb : natural32;
    tolrnk : constant double_float := 1.0E-6;
    nbfail,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    nbtot : constant natural32 := Length_Of(sols);
    fail : boolean;
    sa : Solution_Array(1..integer32(nbtot)) := Create(sols);
    initres : Standard_Floating_Vectors.Vector(sa'range);
    refsols_last : Solution_List;
    eval_acc : Vector(1..nq);
    eva : Floating_Number;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..20)
                      := Multprec_Condition_Tables.Create(20); 

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); new_line(file);
    put(file,nbtot,1); put(file," "); put(file,nv,1); new_line(file);
    Write_Bar(file);
    for i in sa'range loop 
      numb := 0;
      eval_acc := Eval(p_eval,sa(i).v);
      eva := Sum_Norm(eval_acc);
      Clear(eval_acc);
      initres(i) := Trunc(eva);
      if eva < 1.0 then
        if wout then
          Reporting_Gauss_Newton
            (file,p_eval,jac_eval,sa(i).all,tolrnk,epsxa,epsfa,numb,max,fail);
        else
          Silent_Gauss_Newton
            (p_eval,jac_eval,sa(i).all,tolrnk,epsxa,epsfa,numb,max,fail);
        end if;
      else
        fail := true;
      end if;
      Write_Info(file,sa(i).all,initres(i),natural32(i),numb,fail);
      Root_Accounting(file,sa(i),natural32(i),sa(sa'first..i),fail,tolsing,
                      epsxa,nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus);
      if not fail
       then Append(refsols,refsols_last,sa(i).all);
      end if;
      Update_Corrector(t_err,sa(i).all);
      Update_Condition(t_rco,sa(i).all);
      Update_Residuals(t_res,sa(i).all);
      numit := numit + numb;
      Clear(eva);
    end loop;
    Write_Global_Info(file,nbtot,nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus);
    Clear(p_eval); Clear(jac); Clear(jac_eval);
    Shallow_Clear(sols); sols := Create(sa); -- Shallow_Clear(sa);
    Write_Tables(file,t_err,t_res,t_rco);
  end Reporting_Root_Sharpener;

end Multprec_Root_Refiners;
