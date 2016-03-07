with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Natural64_VecVecs;
with Standard_Floating_Vectors;
with Standard_Random_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Standard_Complex_Singular_Values;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Solution_Diagnostics;      use Standard_Solution_Diagnostics;
with Standard_Condition_Tables;
with Standard_Condition_Report;
with Handle_Underflow_Gracefully;
with Standard_Jacobian_Trees;            use Standard_Jacobian_Trees;
with Monomial_Hashing;                   use Monomial_Hashing;
with Standard_Deflation_Trees_io;        use Standard_Deflation_Trees_io;
with Standard_Deflation_Methods;         use Standard_Deflation_Methods;
with Standard_Multiple_Solutions;        use Standard_Multiple_Solutions;
with Standard_Complex_Newton_Steps;      use Standard_Complex_Newton_Steps;

package body Standard_Root_Refiners is

-- ROOT ACCOUNTING :

  procedure Write_Info ( file : in file_type; zero : in Solution;
                         initres : in double_float;
                         i,numb,nbdef : in natural32;
                         fail,infty : in boolean ) is
  begin
    put(file,"solution "); put(file,i,1); put(file," : ");
    put(file,"   start residual : "); put(file,initres,2,3,3);
    if nbdef = 0
     then put(file,"   #iterations : "); put(file,numb,1);
     else put(file,"   #deflations : "); put(file,nbdef,1);
    end if;
    if infty then
      put_line(file,"   at infinity");
    elsif fail then
      put_line(file,"   failure");
    else
      put_line(file,"   success");
    end if;
    put(file,zero);
  end Write_Info;

  procedure Root_Accounting
               ( file : in file_type;
                 h1,h2 : in Standard_Complex_Vectors.Vector;
                 pl : in out Point_List; ls : in Link_to_Solution;
                 nb : in natural32; sa : in out Solution_Array;
                 fail,infty,deflate : in boolean;
                 tolsing,tolclus : in double_float; nbfail,nbinfty,
                 nbreal,nbcomp,nbreg,nbsing,nbclus : in out natural32 ) is
  begin
    if infty then
      put_line(file," at infinity ==");
      nbinfty := nbinfty + 1;
    elsif fail then
      put_line(file," no solution ==");
      nbfail := nbfail + 1;
      ls.m := 0;
    else
      if Is_Real(ls.all,1.0E-13)
       then put(file," real ");    nbreal := nbreal + 1;
       else put(file," complex "); nbcomp := nbcomp + 1;
      end if;
      if deflate or sa(integer32(nb)).rco < tolsing then
        if ls.m <= 1 then  -- do not change input multiplicity field
          declare       -- to determine multiplicity, count clusters
            m : natural32; -- := Multiplicity(ls.all,sa,tolclus);
          begin
            Standard_Condition_Report.Multiplicity
              (ls.all,nb,sa,tolclus,h1,h2,pl,m);
            if ((m = 1) and (not deflate))
             then m := 0;
            end if;
            ls.m := integer32(m);
          end;
        end if;
        if deflate then
          if ls.m = 1
           then put_line(file,"single ==");
           else put_line(file,"multiple ==");
                nbsing := nbsing + 1;
          end if;
          nbreg := nbreg + 1;
        else
          put_line(file,"singular ==");
          nbsing := nbsing + 1;
        end if;
      elsif sa(integer32(nb)).rco > tolsing then
        declare
         -- nb2 : constant natural32 := Is_Clustered(ls.all,nb,sa,tolclus);
          nb2 : natural32;
        begin
          Standard_Condition_Report.Is_Clustered
            (ls.all,nb,sa,tolclus,h1,h2,pl,nb2);
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
  end Root_Accounting;

  procedure Write_Type
               ( file : in file_type; 
                 h1,h2 : in Standard_Complex_Vectors.Vector;
                 pl : in out Point_List; ls : in Link_to_Solution;
                 nb : in natural32; sa : in out Solution_Array;
                 fail,infty,deflate : in boolean;
                 tolsing,tolclus : in double_float; nbfail,nbinfty,
                 nbreal,nbcomp,nbreg,nbsing,nbclus : in out natural32 ) is

    nb2 : natural32;

  begin
    if infty then
      put_line(file," at infinity ==");
      nbinfty := nbinfty + 1;
    elsif fail then
      put_line(file," no solution ==");
      nbfail := nbfail + 1;
      ls.m := 0;
    else
      if Is_Real(ls.all,1.0E-13)
       then put(file," real ");    nbreal := nbreal + 1;
       else put(file," complex "); nbcomp := nbcomp + 1;
      end if;
      if deflate or sa(integer32(nb)).rco < tolsing then
        if deflate then
          if ls.m = 1 then
            put_line(file,"single ==");
          else
            nb2 := Is_Clustered(ls.all,nb,sa,tolclus);
           -- Standard_Condition_Report.Is_Clustered
           --   (ls.all,nb,sa,tolclus,h1,h2,pl,nb2);
            if nb2 = nb then
              put_line(file,"multiple ==");
            else
              put(file,"multiple: ");
              put(file,nb2,1); put_line(file," ==");
              nbsing := nbsing + 1;
            end if;
          end if;
          nbreg := nbreg + 1;
        else
          put_line(file,"singular ==");
          nbsing := nbsing + 1;
        end if;
      elsif sa(integer32(nb)).rco > tolsing then
        if ls.m > 0 then
          put_line(file,"regular ==");
          nbreg := nbreg + 1;
        else
          put(file,"clustered : ");
          put(file,ls.m,1);
          put_line(file," ==");
          nbclus := nbclus + 1;
        end if;
      end if;  
    end if;
  end Write_Type;

  procedure Write_Type
               ( file : in file_type; ls : in Link_to_Solution;
                 fail,infty : in boolean;
                 tolsing : in double_float;
                 nbfail,nbinfty : in out natural32;
                 nbreal,nbcomp,nbreg,nbsing : in out natural32 ) is
  begin
    if infty then
      put_line(file," at infinity =="); nbinfty := nbinfty + 1;
    elsif fail then
      put_line(file," no solution =="); nbfail := nbfail + 1;
      ls.m := 0;
    else
      if Is_Real(ls.all,1.0E-13)
       then put(file," real ");    nbreal := nbreal + 1;
       else put(file," complex "); nbcomp := nbcomp + 1;
      end if;
      if ls.rco < tolsing
       then put_line(file,"singular =="); nbsing := nbsing + 1;
       else put_line(file,"regular ==");  nbreg := nbreg + 1;
      end if;  
    end if;
  end Write_Type;

  procedure Multiplicity
                ( h1,h2 : in Standard_Complex_Vectors.Vector;
                  pl : in out Point_List; ls : in Link_to_Solution;
                  nb : in natural32; sa : in out Solution_Array;
                  fail,infty,deflate : in boolean;
                  tolsing,tolclus : in double_float ) is
  begin
    if infty then
      null;
    elsif fail then
      ls.m := 0;
    elsif sa(integer32(nb)).rco < tolsing or deflate then
      if ls.m <= 1 then -- do not change input multiplicity field
        declare   -- to determine multiplicity count clustered solutions
         -- m : constant natural32 := Multiplicity(ls.all,sa,tolclus);
          m : natural32;
        begin
          Standard_Condition_Report.Multiplicity
            (ls.all,nb,sa,tolclus,h1,h2,pl,m);
          if ((m = 1) and (not deflate))
           then ls.m := 0;
           else ls.m := integer32(m);
          end if;
        end;
      end if;
    else  -- sa(nb).rco > tolsing, check for clustering
      declare
        nb2 : constant natural32 := Is_Clustered(ls.all,nb,sa,tolclus);
      begin
        if nb2 /= nb
         then ls.m := -integer32(nb2); sa(integer32(nb2)).m := -integer32(nb);
        end if;
      end;	   
    end if;
  end Multiplicity;

  procedure Multiplicity
                ( h1,h2 : in Standard_Complex_Vectors.Vector;
                  pl : in out Point_List; ls : in Link_to_Solution;
                  nb : in natural32; sols : in out Solution_List;
                  fail,infty,deflate : in boolean;
                  tolsing,tolclus : in double_float ) is
  begin
    if infty then
      null;
    elsif fail then
      ls.m := 0;
    elsif ls.rco < tolsing or deflate then
      if ls.m <= 1 then -- do not change input multiplicity field
        declare   -- to determine multiplicity count clustered solutions
         -- m : constant natural32 := Multiplicity(ls.all,sols,tolclus);
          m : natural32;
        begin
          Standard_Condition_Report.Multiplicity
            (ls.all,nb,sols,tolclus,h1,h2,pl,m);
          if ((m = 1) and (not deflate))
           then ls.m := 0;
           else ls.m := integer32(m);
          end if;
        end;
      end if;
    else  -- ls.rco > tolsing, check for clustering
      declare
       -- nb2 : constant natural32 := Is_Clustered(ls.all,nb,sols,tolclus);
        nb2 : natural32;
      begin
        Standard_Condition_Report.Is_Clustered
          (ls.all,nb,sols,tolclus,h1,h2,pl,nb2);
        if nb2 /= nb then
          ls.m := -integer32(nb2);
          Change_Multiplicity(sols,nb2,-integer32(nb));
        end if;
      end;	   
    end if;
  end Multiplicity;

  procedure Write_Global_Info
             ( file : in file_type; tot,nbfail,nbinfty,
               nbreal,nbcomp,nbreg,nbsing,nbclus : in natural32 ) is
  begin
    put_bar(file);
    put(file,"A list of "); put(file,tot,1);
    put_line(file," solutions has been refined :");
    put(file,"Number of regular solutions     : "); put(file,nbreg,1);
    put_line(file,".");
    put(file,"Number of singular solutions    : "); put(file,nbsing,1);
    put_line(file,".");
    put(file,"Number of real solutions        : "); put(file,nbreal,1);
    put_line(file,".");
    put(file,"Number of complex solutions     : "); put(file,nbcomp,1);
    put_line(file,".");
    put(file,"Number of clustered solutions   : "); put(file,nbclus,1);
    put_line(file,".");
    put(file,"Number of solutions at infinity : "); put(file,nbinfty,1);
    put_line(file,".");
    put(file,"Number of failures              : "); put(file,nbfail,1);
    put_line(file,".");
    put_bar(file);
  end Write_Global_Info;

-- ONE NEWTON STEP :

  procedure Standard_SVD_Newton_Step
              ( f : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float ) is

    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;

    y : Standard_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,x'range) := eval(jf,x);
    n : constant integer32 := f'last;
    p : constant integer32 := x'last;
    dx : Standard_Complex_Vectors.Vector(1..p);
    mm : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    sv : Standard_Complex_Vectors.Vector(1..mm);
    e : Standard_Complex_Vectors.Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    Standard_Complex_Singular_Values.SVD(A,n,p,sv,e,u,v,job,info);
    rco := Standard_Complex_Singular_Values.Inverse_Condition_Number(sv);
    Standard_Complex_Vectors.Min(y);
    dx := Standard_Complex_Singular_Values.Solve(u,v,sv,y);
    Standard_Complex_Vectors.Add(x,dx);
    err := Max_Norm(dx);
    y := eval(f,x);
    res := Max_Norm(y);
  end Standard_SVD_Newton_Step;

  procedure Standard_SVD_Newton_Step
              ( f : in Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float ) is

    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;

    y : Standard_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,x'range) := eval(jf,x);
    n : constant integer32 := f'last;
    p : constant integer32 := x'last;
    dx : Standard_Complex_Vectors.Vector(1..p);
    mm : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    sv : Standard_Complex_Vectors.Vector(1..mm);
    e : Standard_Complex_Vectors.Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    Standard_Complex_Singular_Values.SVD(A,n,p,sv,e,u,v,job,info);
    rco := Standard_Complex_Singular_Values.Inverse_Condition_Number(sv);
    Standard_Complex_Vectors.Min(y);
    dx := Standard_Complex_Singular_Values.Solve(u,v,sv,y);
    Standard_Complex_Vectors.Add(x,dx);
    err := Max_Norm(dx);
    y := eval(f,x);
    res := Max_Norm(y);
  end Standard_SVD_Newton_Step;

  procedure Standard_LU_Newton_Step
              ( f : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float ) is

    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;

    y : Standard_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,f'range) := eval(jf,x);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : constant double_float := Norm1(A);

  begin
    Standard_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,rco);
    lusolve(A,A'last(1),ipvt,y);
    Standard_Complex_Vectors.Add(x,y);
    err := Max_Norm(y);
    y := eval(f,x);
    res := Max_Norm(y);
  end Standard_LU_Newton_Step;

  procedure Standard_LU_Newton_Step
              ( f : in Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float ) is

    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;

    y : Standard_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,f'range) := eval(jf,x);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : constant double_float := Norm1(A);

  begin
    Standard_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,rco);
    lusolve(A,A'last(1),ipvt,y);
    Standard_Complex_Vectors.Add(x,y);
    err := Max_Norm(y);
    y := eval(f,x);
    res := Max_Norm(y);
  end Standard_LU_Newton_Step;

-- WRAPPING ONE NEWTON STEP :

  procedure Standard_Newton_Step
              ( f : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float ) is
  begin
    if f'last > x'last
     then Standard_SVD_Newton_Step(f,jf,x,err,rco,res);
     else Standard_LU_Newton_Step(f,jf,x,err,rco,res);
    end if;
  end Standard_Newton_Step;

  procedure Standard_Newton_Step
              ( f : in Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out Standard_Complex_Vectors.Vector;
                err,rco,res : out double_float ) is
  begin
    if f'last > x'last
     then Standard_SVD_Newton_Step(f,jf,x,err,rco,res);
     else Standard_LU_Newton_Step(f,jf,x,err,rco,res);
    end if;
  end Standard_Newton_Step;

-- TARGET ROUTINES, NEWTON's METHOD :

  procedure Silent_Newton
               ( p_eval : in Eval_Poly_Sys;
                 j_eval : in Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                 zero : in out Solution; epsxa,epsfa : in double_float; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    use Standard_Complex_Jaco_Matrices;

    n : constant integer32 := zero.n;
    jacobian : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    y,deltax : Vector(1..n);

  begin
    y := Eval(p_eval,zero.v);               -- y = f(zero)
    for i in 1..max loop
      jacobian := Eval(j_eval,zero.v);      -- solve jacobian*deltax = -f(zero)
      lufco(jacobian,n,ipvt,zero.rco);
      if zero.rco < 1.0E-14 then
        zero.res := Sum_Norm(y);
        fail := (zero.res > epsfa);
        return;                             -- singular Jacobian matrix
      end if;
      deltax := -y;
      lusolve(jacobian,n,ipvt,deltax);  
      Add(zero.v,deltax);                   -- make the updates
      y := Eval(p_eval,zero.v);
      zero.err := Sum_Norm(deltax); zero.res := Sum_Norm(y);
      numit := numit + 1;
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) then -- stop ?
        fail := false; exit;
      elsif numit >= max then
        fail := true; exit;
      end if;
    end loop;
    jacobian := Eval(j_eval,zero.v);        -- compute condition number
    lufco(jacobian,n,ipvt,zero.rco);
  exception
   -- when constraint_error => fail := (zero.res > epsfa); return;
    when others =>
      --put_line("exception in silent Newton 1");
      fail := (zero.res > epsfa);
      --if not fail
      -- then
      Handle_Underflow_Gracefully.Underflow_to_Zero(zero.v);
      --put_line("setting to underflow done");
      -- else put_line("fail is true");
      --end if;
      return;
  end Silent_Newton;

  procedure Silent_Newton
               ( p_eval : in Eval_Laur_Sys;
                 j_eval : in Standard_Complex_Laur_Jacomats.Eval_Jaco_Mat;
                 zero : in out Solution; epsxa,epsfa : in double_float; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    use Standard_Complex_Laur_Jacomats;

    n : constant integer32 := zero.n;
    jacobian : Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    y,deltax : Vector(1..n);

  begin
    y := Eval(p_eval,zero.v);               -- y = f(zero)
    for i in 1..max loop
      jacobian := Eval(j_eval,zero.v);      -- solve jacobian*deltax = -f(zero)
      lufco(jacobian,n,ipvt,zero.rco);
      if zero.rco < 1.0E-14 then
        zero.res := Sum_Norm(y);
        fail := (zero.res > epsfa);
        return;                             -- singular Jacobian matrix
      end if;
      deltax := -y;
      lusolve(jacobian,n,ipvt,deltax);  
      Add(zero.v,deltax);                   -- make the updates
      y := Eval(p_eval,zero.v);
      zero.err := Sum_Norm(deltax); zero.res := Sum_Norm(y);
      numit := numit + 1;
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) then -- stop ?
        fail := false; exit;
      elsif numit >= max then
        fail := true; exit;
      end if;
    end loop;
    jacobian := Eval(j_eval,zero.v);        -- compute condition number
    lufco(jacobian,n,ipvt,zero.rco);
  exception
   -- when constraint_error => fail := (zero.res > epsfa); return;
    when others =>
      --put_line("exception in silent Newton 1");
      fail := (zero.res > epsfa);
      --if not fail
      -- then
      Handle_Underflow_Gracefully.Underflow_to_Zero(zero.v);
      --put_line("setting to underflow done");
      -- else put_line("fail is true");
      --end if;
      return;
  end Silent_Newton;

  procedure Silent_Newton
               ( p_eval : in Standard_Complex_Poly_SysFun.Evaluator;
                 j_eval : in Standard_Complex_Jaco_Matrices.Evaluator;
                 zero : in out Solution; epsxa,epsfa : in double_float;
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    n : constant integer32 := zero.n;
    jacobian : matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    y,deltax : Vector(1..n);

  begin
    y := p_eval(zero.v);                    -- y = f(zero)
    for i in 1..max loop
      jacobian := j_eval(zero.v);           -- solve jacobian*deltax = -f(zero)
      lufco(jacobian,n,ipvt,zero.rco);
      if zero.rco < 1.0E-14 then
        zero.res := Sum_Norm(y); 
        fail := (zero.res > epsfa);
        return;                             -- singular Jacobian matrix
      end if;
      deltax := -y;
      lusolve(jacobian,n,ipvt,deltax);
      Add(zero.v,deltax);                   -- make the updates
      y := p_eval(zero.v);
      zero.err := Sum_Norm(deltax); zero.res := Sum_Norm(y);
      numit := numit + 1;
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) then -- stop ?
        fail := false; exit;
      elsif numit >= max then
        fail := true; exit;
      end if;
    end loop;
    jacobian := j_eval(zero.v);            -- compute condition number
    lufco(jacobian,n,ipvt,zero.rco);
  exception
   -- when constraint_error => fail := (zero.res > epsfa); return;
    when others =>
      --put_line("exception in silent Newton 2");
      fail := (zero.res > epsfa);
      --if not fail
      -- then
      Handle_Underflow_Gracefully.Underflow_to_Zero(zero.v);
      --put_line("setting to underflow done");
      -- else put_line("fail is true");
      --end if;
      return;
  end Silent_Newton;

  procedure Reporting_Newton
               ( file : in file_type; p_eval : in Eval_Poly_Sys;
                 j_eval : in Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                 zero : in out Solution; epsxa,epsfa : in double_float; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    use Standard_Complex_Jaco_Matrices;

    n : constant integer32 := zero.n;
    jacobian : matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    y,deltax : Vector(1..n);

  begin
    y := Eval(p_eval,zero.v);              -- y = f(zero)
    for i in 1..max loop
      jacobian := Eval(j_eval,zero.v);     -- solve jacobian*deltax = -f(zero)
      lufco(jacobian,n,ipvt,zero.rco);
      put(file,"zero.rco : "); put(file,zero.rco,3); new_line(file);
      if zero.rco < 1.0E-14                -- singular jacobian matrix
       then zero.res := Sum_Norm(y); 
            fail := (zero.res > epsfa);    -- accuracy not reached yet
            return;
      end if;
      deltax := -y;
      lusolve(jacobian,n,ipvt,deltax);  
      Add(zero.v,deltax);                  -- make the updates
      y := Eval(p_eval,zero.v);
      zero.err := Sum_Norm(deltax); zero.res := Sum_Norm(y);
      numit := numit + 1;
      put(file,"Step "); put(file,numit,4); new_line(file);      -- output
      put(file," |errxa| : "); put(file,zero.err); new_line(file);
      put(file," |errfa| : "); put(file,zero.res); new_line(file);
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) then -- stop ?
        fail := false; exit;
      elsif numit >= max then
        fail := true; exit;
      end if;
    end loop;
    jacobian := Eval(j_eval,zero.v);              -- compute condition number
    lufco(jacobian,n,ipvt,zero.rco);
  exception
   -- when constraint_error => fail := (zero.res > epsfa); return;
    when others => 
      --put_line("exception in reporting Newton 1");
      fail := (zero.res > epsfa);
      --if not fail
      -- then
      Handle_Underflow_Gracefully.Underflow_to_Zero(zero.v);
      --put_line("underflow handling done");
      -- else put_line("fail is true");
      --end if;
      return;
  end Reporting_Newton;

  procedure Reporting_Newton
               ( file : in file_type; p_eval : in Eval_Laur_Sys;
                 j_eval : in Standard_Complex_Laur_Jacomats.Eval_Jaco_Mat;
                 zero : in out Solution; epsxa,epsfa : in double_float; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    use Standard_Complex_Laur_Jacomats;

    n : constant integer32 := zero.n;
    jacobian : matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    y,deltax : Vector(1..n);

  begin
    y := Eval(p_eval,zero.v);              -- y = f(zero)
    for i in 1..max loop
      jacobian := Eval(j_eval,zero.v);     -- solve jacobian*deltax = -f(zero)
      lufco(jacobian,n,ipvt,zero.rco);
      put(file,"zero.rco : "); put(file,zero.rco,3); new_line(file);
      if zero.rco < 1.0E-14                -- singular jacobian matrix
       then zero.res := Sum_Norm(y); 
            fail := (zero.res > epsfa);    -- accuracy not reached yet
            return;
      end if;
      deltax := -y;
      lusolve(jacobian,n,ipvt,deltax);  
      Add(zero.v,deltax);                  -- make the updates
      y := Eval(p_eval,zero.v);
      zero.err := Sum_Norm(deltax); zero.res := Sum_Norm(y);
      numit := numit + 1;
      put(file,"Step "); put(file,numit,4); new_line(file);      -- output
      put(file," |errxa| : "); put(file,zero.err); new_line(file);
      put(file," |errfa| : "); put(file,zero.res); new_line(file);
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) then -- stop ?
        fail := false; exit;
      elsif numit >= max then
        fail := true; exit;
      end if;
    end loop;
    jacobian := Eval(j_eval,zero.v);              -- compute condition number
    lufco(jacobian,n,ipvt,zero.rco);
  exception
   -- when constraint_error => fail := (zero.res > epsfa); return;
    when others => 
      --put_line("exception in reporting Newton 1");
      fail := (zero.res > epsfa);
      --if not fail
      -- then
      Handle_Underflow_Gracefully.Underflow_to_Zero(zero.v);
      --put_line("underflow handling done");
      -- else put_line("fail is true");
      --end if;
      return;
  end Reporting_Newton;

  procedure Reporting_Newton
               ( file : in file_type;
                 p_eval : in Standard_Complex_Poly_SysFun.Evaluator;
                 j_eval : in Standard_Complex_Jaco_Matrices.Evaluator;
                 zero : in out Solution; epsxa,epsfa : in double_float;
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    n : constant integer32 := zero.n;
    jacobian : matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    y,deltax : Vector(1..n);

  begin
    y := p_eval(zero.v);                   -- y = f(zero)
    for i in 1..max loop
      jacobian := j_eval(zero.v);          -- solve jacobian*deltax = -f(zero)
      lufco(jacobian,n,ipvt,zero.rco);
      if zero.rco < 1.0E-14                -- singular jacobian matrix
       then zero.res := Sum_Norm(y); 
            fail := (zero.res > epsfa); -- accuracy not reached yet
            return;
      end if;
      deltax := -y;
      lusolve(jacobian,n,ipvt,deltax);
      Add(zero.v,deltax);                  -- make the updates
      y := p_eval(zero.v);
      zero.err := Sum_Norm(deltax); zero.res := Sum_Norm(y);
      numit := numit + 1;
      put(file,"Step "); put(file,numit,4); new_line(file);      -- output
      put(file," |errxa| : "); put(file,zero.err); new_line(file);
      put(file," |errfa| : "); put(file,zero.res); new_line(file);
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) then
        fail := false; exit;                  -- stopping criteria
      elsif numit >= max then
        fail := true; exit;
      end if;
    end loop;
    jacobian := j_eval(zero.v);                   -- compute condition number
    lufco(jacobian,n,ipvt,zero.rco);
  exception
   -- when constraint_error => fail := (zero.res > epsfa); return;
    when others =>
      --put_line("exception in reporting Newton 2");
      fail := (zero.res > epsfa);
      --if not fail
      -- then
      Handle_Underflow_Gracefully.Underflow_to_Zero(zero.v);
      --put_line("underflow setting is done");
      -- else put_line("fail is true");
      --end if;
      return;
  end Reporting_Newton;

-- DEFLATION PROCEDURES :

  procedure Silent_Deflate
              ( max : in natural32; f : Eval_Poly_Sys;
                jf : Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                ls : in out Link_to_Solution; order : in integer32;
                tolrnk : in double_float; nd : in Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                nv,nq,R1 : in out Standard_Natural_Vectors.Vector;
                nbitr,nbdef : out natural32; fail : out boolean ) is

  -- DECRIPTION :
  --   Performs deflation on a system without intermediate output.

    use Standard_Complex_Jaco_Matrices;

    k : integer32 := 0;
    rank,m : natural32;
    done : boolean := false;
    B : Standard_Complex_VecMats.VecMat(1..order);
    h : Standard_Complex_VecVecs.VecVec(1..order);
    x : Standard_Complex_VecVecs.VecVec(0..order);
    rsd : double_float;

  begin
    nbitr := 0;
    nbdef := 0;
    loop
      if k = 0 then
        Apply_Newton(max,f,jf,ls,tolrnk,rank);
        done := (integer32(rank) = ls.n);
        if not done then
          x(0) := new Standard_Complex_Vectors.Vector'(ls.v);
        end if;
        nbitr := 1;
      else
        Apply_Newton(max,f,jf,nd,monkeys,k,nv,nq,R1,B,h,x,
                     ls.err,ls.rco,ls.res,tolrnk,rank);
        ls.v := x(0).all;
        done := (rank = nv(k));
      end if;
      exit when done or (k >= order);
      k := k + 1; m := rank + 1;
      Add_Deflation_Data(k,m,nv,nq,R1,B,h);
      Add_Multipliers(f,jf,nd,monkeys,k,nv(0..k),
                      nq(0..k),R1(1..k),B(1..k),h(1..k),x(0..k),rsd);
      exit when (rsd > 0.1); -- deflation fails
    end loop;
    Standard_Complex_Vectors.Clear(x(0));
    for i in 1..k loop
      Standard_Complex_Matrices.Clear(B(i)); 
      Standard_Complex_Vectors.Clear(h(i));
      Standard_Complex_Vectors.Clear(x(i));
    end loop;
    fail := not done;
    nbdef := natural32(k);
  exception
   -- when constraint_error -- same underflow as in Silent_Newton may occur
    when others 
       => --put_line("exception raised in silent deflate");
          Handle_Underflow_Gracefully.Underflow_to_Zero(ls.v);
          ls.rco := Handle_Underflow_Gracefully.Underflow_to_Zero(ls.rco);
          fail := not done;
  end Silent_Deflate;

  procedure Reporting_Deflate
              ( file : in file_type; output : in boolean; 
                max : in natural32; f : Eval_Poly_Sys;
                jf : Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                ls : in out Link_to_Solution; order : in integer32;
                tolrnk : in double_float; nd : in Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                nv,nq,R1 : in out Standard_Natural_Vectors.Vector;
                nbitr,nbdef : out natural32; fail : out boolean ) is

  -- DECRIPTION :
  --   Performs deflation on a system with intermediate output.

    use Standard_Complex_Jaco_Matrices;

    k : integer32 := 0;
    rank,m : natural32;
    done : boolean := false;
    B : Standard_Complex_VecMats.VecMat(1..order);
    h : Standard_Complex_VecVecs.VecVec(1..order);
    x : Standard_Complex_VecVecs.VecVec(0..order);
    rsd : double_float;

  begin
    nbitr := 0;
    nbdef := 0;
    loop
      if k = 0 then
        Apply_Newton(file,output,max,f,jf,ls,tolrnk,rank);
        done := (integer32(rank) = ls.n);
        if not done then
          x(0) := new Standard_Complex_Vectors.Vector'(ls.v);
        end if;
        nbitr := 1;
      else
        Apply_Newton(file,output,max,f,jf,nd,monkeys,k,nv,nq,
                     R1,B,h,x,ls.err,ls.rco,ls.res,tolrnk,rank);
        ls.v := x(0).all;
        done := (rank = nv(k));
      end if;
      exit when done or (k >= order);
      k := k + 1; m := rank + 1;
      Add_Multiplier_Symbols(natural32(k),m);
      Add_Deflation_Data(k,m,nv,nq,R1,B,h);
      Add_Multipliers(file,output,f,jf,nd,monkeys,k,nv(0..k),
                      nq(0..k),R1(1..k),B(1..k),h(1..k),x(0..k),rsd);
      exit when (rsd > 0.1); -- deflation fails
    end loop;
    Standard_Complex_Vectors.Clear(x(0));
    for i in 1..k loop
      Standard_Complex_Matrices.Clear(B(i)); 
      Standard_Complex_Vectors.Clear(h(i));
      Standard_Complex_Vectors.Clear(x(i));
    end loop;
    fail := not done;
    nbdef := natural32(k);
  exception
   -- when constraint_error -- same underflow as in Reporting_Newton may occur
    when others -- 
       => -- put_line("exception raised in reporting deflate");
          Handle_Underflow_Gracefully.Underflow_to_Zero(ls.v);
          ls.rco := Handle_Underflow_Gracefully.Underflow_to_Zero(ls.rco);
          fail := not done;
  end Reporting_Deflate;

-- APPLICATION OF NEWTON's METHOD ON LIST OF SOLUTIONS :

  procedure Silent_Root_Refiner
               ( p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in out boolean ) is

    use Standard_Complex_Jaco_Matrices;

    n : constant integer32 := p'length;
    p_eval : Eval_Poly_Sys(1..n) := Create(p);
    jac : Jaco_Mat(1..n,1..n) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..n,1..n) := Create(jac);
    numb,nbdef,nit : natural32 := 0;
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(Length_Of(sols))) := Create(sols);
    tolrnk : constant double_float := 100.0*tolsing;
    order : constant integer32 := 3;
    nv : Standard_Natural_Vectors.Vector(0..order);
    nq : Standard_Natural_Vectors.Vector(0..order);
    R1 : Standard_Natural_Vectors.Vector(1..order);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..order);
    nd : Link_to_Eval_Node;
    backup : Solution(n);
    merge : boolean := false; -- to merge clustered solutions
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;

  begin
    if deflate then
      declare
      begin
        nv(0) := natural32(n); nq(0) := natural32(p'last); R1(1) := 0;
        Create_Remember_Derivatives(jac,order,nd);
        monkeys := Monomial_Keys(natural32(order),nv(0));
      exception
        when others =>
          put_line("The system is too large to apply deflation.");
          deflate := false;
      end;
    end if;
    for i in sa'range loop
      numb := 0; nbdef := 0;
      sa(i).res := Sum_Norm(Eval(p_eval,sa(i).v));
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        if deflate then
          backup := sa(i).all;
          Silent_Deflate
            (max,p_eval,jac_eval,sa(i),order,tolrnk,nd,monkeys,nv,nq,R1,
             numb,nbdef,fail);
          Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,nit,max,fail);
          if fail and backup.res < sa(i).res then
            sa(i).all := backup;
            Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,nit,max,fail);
          end if;
        else
          Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
        end if;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,deflate,tolsing,epsxa);
      if not fail and then deflate
       then merge := merge or (sa(i).m > 1);
      end if;
      numit := numit + numb;
    end loop;
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
    if deflate then
      Standard_Natural64_VecVecs.Clear(monkeys);
      if merge
       then Merge_Multiple_Solutions(sols,tolsing);
      end if;
      Standard_Jacobian_Trees.Clear(nd);
    else
      Clear(jac); -- otherwise crash after Clear(nd)
    end if;
    Clear(p_eval); Clear(jac_eval);
    Clear(pl);
  --exception
  --  when others => put_line("exception raised in silent root refiner 1");
  --                 raise;
  end Silent_Root_Refiner;

  procedure Silent_Root_Refiner
               ( p : in Standard_Complex_Poly_SysFun.Evaluator;
                 j : in Standard_Complex_Jaco_Matrices.Evaluator;
                 sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32 ) is

    numb : natural32;
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(Length_Of(sols))) := Create(sols);
    n : constant integer32 := Head_Of(sols).n;
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;

  begin
    for i in sa'range loop
      numb := 0;
      sa(i).res := Sum_Norm(p(sa(i).v));
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1
       then Silent_Newton(p,j,sa(i).all,epsxa,epsfa,numb,max,fail);
       else fail := true;
      end if;
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,false,tolsing,epsxa);
      numit := numit + numb;
    end loop;
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
    Clear(pl);
  end Silent_Root_Refiner;

  procedure Silent_Root_Refiner
               ( p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in out boolean ) is

    use Standard_Complex_Jaco_Matrices;

    n : constant integer32 := p'last;
    p_eval : Eval_Poly_Sys(1..n) := Create(p);
    jac : Jaco_Mat(1..n,1..n) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..n,1..n) := Create(jac);
    numb,nbdef,nit : natural32 := 0;
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(Length_Of(sols))) := Create(sols);
    refsols_last : Solution_List;
    tolrnk : constant double_float := 100.0*tolsing;
    order : constant integer32 := 3;
    nv : Standard_Natural_Vectors.Vector(0..order);
    nq : Standard_Natural_Vectors.Vector(0..order);
    R1 : Standard_Natural_Vectors.Vector(1..order);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..order);
    nd : Link_to_Eval_Node;
    backup : Solution(n);
    merge : boolean := false;  -- flag to merge clustered solutions
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;

  begin
    if deflate then
      declare
      begin
        nv(0) := natural32(n); nq(0) := natural32(p'last); R1(1) := 0;
        Create_Remember_Derivatives(jac,order,nd);
        monkeys := Monomial_Keys(natural32(order),nv(0));
      exception
        when others =>
          put_line("The system is too large to apply deflation.");
          deflate := false;
      end;
    end if;
    for i in sa'range loop
      numb := 0; nbdef := 0;
      sa(i).res := Sum_Norm(Eval(p_eval,sa(i).v));
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        if deflate then
          backup := sa(i).all;
          Silent_Deflate
            (max,p_eval,jac_eval,sa(i),order,tolrnk,nd,monkeys,nv,nq,R1,
             numb,nbdef,fail);
          Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,nit,max,fail);
          if fail and backup.res < sa(i).res then
            sa(i).all := backup;
            Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,nit,max,fail);
          end if;
        else
          Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
        end if;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,deflate,tolsing,epsxa);
      if not fail then
        Append(refsols,refsols_last,sa(i).all);
        if deflate
         then merge := merge or (sa(i).m > 1);
        end if;
      end if;
      numit := numit + numb;
    end loop;
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
    if deflate then
      Standard_Natural64_VecVecs.Clear(monkeys);
      Standard_Jacobian_Trees.Clear(nd);
      if merge then
        Merge_Multiple_Solutions(refsols,tolsing);
        Merge_Multiple_Solutions(sols,tolsing);
      end if;
    else
      Clear(jac); -- otherwise crash after clear(nd)
    end if;
    Clear(p_eval); Clear(jac_eval);
    Clear(pl);
  --exception
  --  when others => put_line("exception raised in silent root refiner 2");
  --                 raise;
  end Silent_Root_Refiner;

  procedure Silent_Root_Refiner
               ( p : in Laur_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32 ) is

    use Standard_Complex_Laur_Jacomats;

    n : constant integer32 := p'last;
    p_eval : Eval_Laur_Sys(1..n) := Create(p);
    jac : Jaco_Mat(1..n,1..n) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..n,1..n) := Create(jac);
    numb : natural32 := 0;
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(Length_Of(sols))) := Create(sols);
    refsols_last : Solution_List;
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;

  begin
    for i in sa'range loop
      numb := 0;
      sa(i).res := Sum_Norm(Eval(p_eval,sa(i).v));
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,false,tolsing,epsxa);
      if not fail
       then Append(refsols,refsols_last,sa(i).all);
      end if;
      numit := numit + numb;
    end loop;
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
    Clear(jac); Clear(p_eval); Clear(jac_eval);
    Clear(pl);
  end Silent_Root_Refiner;

  procedure Silent_Root_Refiner
               ( p : in Standard_Complex_Poly_SysFun.Evaluator;
                 j : in Standard_Complex_Jaco_Matrices.Evaluator;
                 sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32 ) is

    numb : natural32;
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(Length_Of(sols))) := Create(sols);
    refsols_last : Solution_List;
    n : constant integer32 := Head_Of(sols).n;
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;

  begin
    for i in sa'range loop
      numb := 0;
      sa(i).res := Sum_Norm(p(sa(i).v));
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if infty and sa(i).res < 0.1
       then Silent_Newton(p,j,sa(i).all,epsxa,epsfa,numb,max,fail);
       else fail := true;
      end if;
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,false,tolsing,epsxa);
      if not fail
       then Append(refsols,refsols_last,sa(i).all);
      end if;
      numit := numit + numb;
    end loop;
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
    Clear(pl);
  end Silent_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in out boolean; wout : in boolean ) is

    use Standard_Complex_Jaco_Matrices;

    n : constant integer32 := p'length;
    p_eval : Eval_Poly_Sys(1..n) := Create(p);
    jac : Jaco_Mat(1..n,1..n) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..n,1..n) := Create(jac);
    numb,nbdef,nit : natural32 := 0;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    nbtot : constant natural32 := Length_Of(sols);
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(nbtot)) := Create(sols);
    initres : Standard_Floating_Vectors.Vector(sa'range);
   -- t_err,t_rco,t_dis,t_res : Standard_Natural_Vectors.Vector(0..15)
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15)
                      := Standard_Condition_Tables.Create(15); 
    tolrnk : constant double_float := 100.0*tolsing;
    order : constant integer32 := 3;
    nv : Standard_Natural_Vectors.Vector(0..order);
    nq : Standard_Natural_Vectors.Vector(0..order);
    R1 : Standard_Natural_Vectors.Vector(1..order);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..order);
    nd : Link_to_Eval_Node;
    backup : Solution(n);
    merge : boolean := false;  -- flag to merge clustered solutions
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;

  begin
   -- put_line("in this reporting root refiner ...");
    if deflate then
      declare
      begin
        nv(0) := natural32(n); nq(0) := natural32(p'last); R1(1) := 0;
        Create_Remember_Derivatives(jac,order,nd);
        monkeys := Monomial_Keys(natural32(order),nv(0));
      exception
        when others =>
          put_line("The system is too large to apply deflation.");
          deflate := false;
      end;
    end if;
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); new_line(file);
    put(file,nbtot,1); put(file," "); put(file,n,1); new_line(file);
    put_bar(file);
    for i in sa'range loop 
      numb := 0; nbdef := 0;
      sa(i).res := Sum_Norm(Eval(p_eval,sa(i).v));
      initres(i) := sa(i).res;
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        declare
        begin
          if deflate then
            backup := sa(i).all;
            Reporting_Deflate(file,wout,max,p_eval,jac_eval,sa(i),order,
                              tolrnk,nd,monkeys,nv,nq,R1,numb,nbdef,fail);
            if wout then
              Reporting_Newton(file,p_eval,jac_eval,sa(i).all,epsxa,epsfa,
                               nit,max,fail);
            else
              Silent_Newton
                (p_eval,jac_eval,sa(i).all,epsxa,epsfa,nit,max,fail);
            end if;
            if fail and backup.res < sa(i).res then
              sa(i).all := backup;
              Silent_Newton
                (p_eval,jac_eval,sa(i).all,epsxa,epsfa,nit,max,fail);
            end if;
          elsif wout then
            Reporting_Newton(file,p_eval,jac_eval,sa(i).all,epsxa,epsfa,
                             numb,max,fail);
          else
            Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
          end if;
       -- exception
       --   when others => fail := (sa(i).res > epsfa);
       --      put("exception raised at solution "); put(i,1); new_line;
       --      flush(file); put_line("flushed buffers..."); raise;
        end;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,deflate,tolsing,epsxa);
      Write_Info(file,sa(i).all,initres(i),natural32(i),numb,nbdef,fail,infty);
      Write_Type
        (file,h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,infty,deflate,
         tolsing,epsxa,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
      if not fail and then deflate 
       then merge := merge or (sa(i).m > 1);
      end if;
      Standard_Condition_Tables.Update_Corrector(t_err,sa(i).all);
      Standard_Condition_Tables.Update_Condition(t_rco,sa(i).all);
      Standard_Condition_Tables.Update_Residuals(t_res,sa(i).all);
      numit := numit + numb;
    end loop;
   -- put_line("done with the loop, writing global info ...");
    Write_Global_Info
      (file,nbtot,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
   -- put_line("cleaning up ...");
    if deflate then
      Standard_Natural64_VecVecs.Clear(monkeys);
      Standard_Jacobian_Trees.Clear(nd);
      if merge
       then Merge_Multiple_Solutions(sols,tolsing);
      end if;
    else
      Clear(jac); -- otherwise crash after clear(nd)
    end if;
    Clear(p_eval); Clear(jac_eval);
   -- Distances_Table(t_dis,sols);
   -- put_line("writing the condition tables ...");
    Standard_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco); --,t_dis);
    Clear(pl);
  --exception
  --  when others => put_line("Exception raised in reporting root refiner 1");
  --                 raise;
  end Reporting_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in Standard_Complex_Poly_SysFun.Evaluator;
                 j : in Standard_Complex_Jaco_Matrices.Evaluator;
                 sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean ) is

    n : constant integer32 := Head_Of(sols).n;
    numb : natural32;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    nbtot : constant natural32 := Length_Of(sols);
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(nbtot)) := Create(sols);
    initres : Standard_Floating_Vectors.Vector(sa'range);
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); new_line(file);
    put(file,nbtot,1); put(file," "); put(file,n,1); new_line(file);
    put_bar(file);
    for i in sa'range loop 
      numb := 0;
      sa(i).res := Sum_Norm(p(sa(i).v));
      initres(i) := sa(i).res;
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        if wout
         then Reporting_Newton
                (file,p,j,sa(i).all,epsxa,epsfa,numb,max,fail);
         else Silent_Newton(p,j,sa(i).all,epsxa,epsfa,numb,max,fail);
        end if;
      else
        fail := true;
      end if;
      Write_Info(file,sa(i).all,initres(i),natural32(i),numb,0,fail,infty);
      Root_Accounting
        (file,h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,infty,false,
         tolsing,epsxa,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
      numit := numit + numb;
    end loop;
    Write_Global_Info
      (file,nbtot,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
    Clear(pl);
  end Reporting_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in out boolean; wout : in boolean ) is

    use Standard_Complex_Jaco_Matrices;

    n : constant integer32 := p'length;
    p_eval : Eval_Poly_Sys(1..n) := Create(p);
    eyp : Standard_Complex_Vectors.Vector(p'range);
    jac : Jaco_Mat(1..n,1..n) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..n,1..n) := Create(jac);
    numb,nbdef,nit : natural32 := 0;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    nbtot : constant natural32 := Length_Of(sols);
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(nbtot)) := Create(sols);
    initres : Standard_Floating_Vectors.Vector(sa'range);
    refsols_last : Solution_List;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15)
                      := Standard_Condition_Tables.Create(15); 
    tolrnk : constant double_float := tolsing*100.0;
    order : constant integer32 := 3;
    nv : Standard_Natural_Vectors.Vector(0..order);
    nq : Standard_Natural_Vectors.Vector(0..order);
    R1 : Standard_Natural_Vectors.Vector(1..order);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..order);
    nd : Link_to_Eval_Node;
    backup : Solution(n);
    merge : boolean := false;  -- to merge clustered solutions
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;

  begin
    if deflate then
      declare
      begin
        nv(0) := natural32(n); nq(0) := natural32(p'last); R1(1) := 0;
        Create_Remember_Derivatives(jac,order,nd);
        monkeys := Monomial_Keys(natural32(order),nv(0));
      exception
        when others =>
          put_line("The system is too large to apply deflation.");
          deflate := false;
      end;
    end if;
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); new_line(file);
    put(file,nbtot,1); put(file," "); put(file,n,1); new_line(file);
    put_bar(file);
    for i in sa'range loop 
      numb := 0; nbdef := 0;
      declare
      begin
        eyp := Eval(p_eval,sa(i).v);
      exception
        when others => put_line("exception with eval..."); raise;
      end;
      sa(i).res := Sum_Norm(eyp);
      initres(i) := sa(i).res;
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        if deflate then
          backup := sa(i).all;
          Reporting_Deflate(file,wout,max,p_eval,jac_eval,sa(i),order,
                            tolrnk,nd,monkeys,nv,nq,R1,numb,nbdef,fail);
          if wout then
            Reporting_Newton(file,p_eval,jac_eval,sa(i).all,epsxa,epsfa,
                             nit,max,fail);
          else 
            Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,nit,max,fail);
          end if;
          if fail and backup.res < sa(i).res then
            sa(i).all := backup;
            Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,nit,max,fail);
          end if;
        elsif wout then
          Reporting_Newton(file,p_eval,jac_eval,sa(i).all,epsxa,epsfa,
                           numb,max,fail);
        else 
          Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
        end if;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,deflate,tolsing,epsxa);
      Write_Info(file,sa(i).all,initres(i),natural32(i),numb,nbdef,fail,infty);
      Write_Type
        (file,h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,infty,deflate,
         tolsing,epsxa,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
      if not fail then
        Append(refsols,refsols_last,sa(i).all);
        if deflate
         then merge := merge and (sa(i).m > 1);
        end if;
      end if;
      Standard_Condition_Tables.Update_Corrector(t_err,sa(i).all);
      Standard_Condition_Tables.Update_Condition(t_rco,sa(i).all);
      Standard_Condition_Tables.Update_Residuals(t_res,sa(i).all);
      numit := numit + numb;
    end loop;
    Write_Global_Info
      (file,nbtot,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
    if deflate then
      Standard_Natural64_VecVecs.Clear(monkeys);
      Standard_Jacobian_Trees.Clear(nd);
      if merge then
        Merge_Multiple_Solutions(sols,tolsing);
        Merge_Multiple_Solutions(refsols,tolsing);
      end if;
    else
      Clear(jac); -- otherwise crash after Clear(nd)
    end if;
    Clear(p_eval); Clear(jac_eval);
    Clear(pl);
   -- Distances_Table(t_dis,sols);
    Standard_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco); -- ,t_dis);
  --exception
  --  when others => put_line("exception raised in root refiner 2."); raise;
  end Reporting_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in Laur_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean ) is

    use Standard_Complex_Laur_Jacomats;

    n : constant integer32 := p'length;
    p_eval : Eval_Laur_Sys(1..n) := Create(p);
    eyp : Standard_Complex_Vectors.Vector(p'range);
    jac : Jaco_Mat(1..n,1..n) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..n,1..n) := Create(jac);
    numb,nbdef : natural32 := 0;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    nbtot : constant natural32 := Length_Of(sols);
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(nbtot)) := Create(sols);
    initres : Standard_Floating_Vectors.Vector(sa'range);
    refsols_last : Solution_List;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15)
                      := Standard_Condition_Tables.Create(15); 
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); new_line(file);
    put(file,nbtot,1); put(file," "); put(file,n,1); new_line(file);
    put_bar(file);
    for i in sa'range loop 
      numb := 0; nbdef := 0;
      eyp := Eval(p_eval,sa(i).v);
      sa(i).res := Sum_Norm(eyp);
      initres(i) := sa(i).res;
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        if wout then
          Reporting_Newton(file,p_eval,jac_eval,sa(i).all,epsxa,epsfa,
                           numb,max,fail);
        else 
          Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
        end if;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,false,tolsing,epsxa);
      Write_Info(file,sa(i).all,initres(i),natural32(i),numb,nbdef,fail,infty);
      Write_Type
        (file,h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,infty,false,
         tolsing,epsxa,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
      if not fail
       then Append(refsols,refsols_last,sa(i).all);
      end if;
      Standard_Condition_Tables.Update_Corrector(t_err,sa(i).all);
      Standard_Condition_Tables.Update_Condition(t_rco,sa(i).all);
      Standard_Condition_Tables.Update_Residuals(t_res,sa(i).all);
      numit := numit + numb;
    end loop;
    Write_Global_Info
      (file,nbtot,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
    Clear(jac); Clear(p_eval); Clear(jac_eval);
   -- Distances_Table(t_dis,sols);
    Standard_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco); -- ,t_dis);
    Clear(pl);
  end Reporting_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in Standard_Complex_Poly_SysFun.Evaluator;
                 j : in Standard_Complex_Jaco_Matrices.Evaluator;
                 sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean ) is

    n : constant integer32 := Head_Of(sols).n;
    numb : natural32;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    nbtot : constant natural32 := Length_Of(sols);
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(nbtot)) := Create(sols);
    initres : Standard_Floating_Vectors.Vector(sa'range);
    refsols_last : Solution_List;
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,n);
    pl : Point_List;

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); new_line(file);
    put(file,nbtot,1); put(file," "); put(file,n,1); new_line(file);
    put_bar(file);
    for i in sa'range loop 
      numb := 0;
      sa(i).res := Sum_Norm(p(sa(i).v));
      initres(i) := sa(i).res;
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        if wout
         then Reporting_Newton(file,p,j,sa(i).all,epsxa,epsfa,numb,max,fail);
         else Silent_Newton(p,j,sa(i).all,epsxa,epsfa,numb,max,fail);
        end if;
      else 
        fail := true;
      end if;
      Write_Info(file,sa(i).all,initres(i),natural32(i),numb,0,fail,infty);
      Root_Accounting
        (file,h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,infty,false,
         tolsing,epsxa,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
      if not fail
       then Append(refsols,refsols_last,sa(i).all);
      end if;
      numit := numit + numb;
    end loop;
    Write_Global_Info
      (file,nbtot,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
    Clear(pl);
  end Reporting_Root_Refiner;

-- APPLICATION of Gauss-Newton for overdetermined systems

  procedure Silent_Gauss_Newton
               ( f : in Eval_Poly_Sys;
                 jf : in Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat; 
                 zero : in out Solution; tol,epsxa,epsfa : in double_float; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    rank : natural32;

    function ef ( x : Vector ) return Vector is
    begin
      return Eval(f,x);
    end ef;

    function ejf ( x : Vector ) return Matrix is
    begin
      return Standard_Complex_Jaco_Matrices.Eval(jf,x);
    end ejf;

    procedure Newton is new Silent_Newton_Step(ef,ejf);

  begin
    for i in 1..max loop
      Newton(natural32(f'last),zero.v,tol,zero.err,zero.rco,zero.res,rank);
      numit := i;
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) then
        fail := false; exit;
      elsif numit >= max then
        fail := true; exit;
      end if;
    end loop;
  end Silent_Gauss_Newton;

  procedure Silent_Gauss_Newton
               ( f : in Eval_Laur_Sys;
                 jf : in Standard_Complex_Laur_Jacomats.Eval_Jaco_Mat; 
                 zero : in out Solution; tol,epsxa,epsfa : in double_float; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    rank : natural32;

    function ef ( x : Vector ) return Vector is
    begin
      return Eval(f,x);
    end ef;

    function ejf ( x : Vector ) return Matrix is
    begin
      return Standard_Complex_Laur_Jacomats.Eval(jf,x);
    end ejf;

    procedure Newton is new Silent_Newton_Step(ef,ejf);

  begin
    for i in 1..max loop
      Newton(natural32(f'last),zero.v,tol,zero.err,zero.rco,zero.res,rank);
      numit := i;
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) then
        fail := false; exit;
      elsif numit >= max then
        fail := true; exit;
      end if;
    end loop;
  end Silent_Gauss_Newton;

  procedure Reporting_Gauss_Newton
               ( file : in file_type; f : in Eval_Poly_Sys;
                 jf : in Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat; 
                 zero : in out Solution; tol,epsxa,epsfa : in double_float; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    rank : natural32;

    function ef ( x : Vector ) return Vector is
    begin
      return Eval(f,x);
    end ef;

    function ejf ( x : Vector ) return Matrix is
    begin
      return Standard_Complex_Jaco_Matrices.Eval(jf,x);
    end ejf;

    procedure Newton is new Reporting_Newton_Step(ef,ejf);

  begin
    for i in 1..max loop
      Newton
        (file,natural32(f'last),zero.v,tol,zero.err,zero.rco,zero.res,rank);
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
  end Reporting_Gauss_Newton;

  procedure Reporting_Gauss_Newton
               ( file : in file_type; f : in Eval_Laur_Sys;
                 jf : in Standard_Complex_Laur_Jacomats.Eval_Jaco_Mat; 
                 zero : in out Solution; tol,epsxa,epsfa : in double_float; 
                 numit : in out natural32; max : in natural32;
                 fail : out boolean ) is

    rank : natural32;

    function ef ( x : Vector ) return Vector is
    begin
      return Eval(f,x);
    end ef;

    function ejf ( x : Vector ) return Matrix is
    begin
      return Standard_Complex_Laur_Jacomats.Eval(jf,x);
    end ejf;

    procedure Newton is new Reporting_Newton_Step(ef,ejf);

  begin
    for i in 1..max loop
      Newton
        (file,natural32(f'last),zero.v,tol,zero.err,zero.rco,zero.res,rank);
      put(file,"Step "); put(file,numit,4); new_line(file);      -- output
      put(file," |errxa| : "); put(file,zero.err); new_line(file);
      put(file," |errfa| : "); put(file,zero.res); new_line(file);
      numit := i;
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) then
        fail := false; exit;
      elsif numit >= max then
        fail := true; exit;
      end if;
    end loop;
  end Reporting_Gauss_Newton;
  
  procedure Silent_Root_Sharpener
               ( p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in boolean ) is

    use Standard_Complex_Jaco_Matrices;

    nq : constant integer32 := p'last;
    nv : constant integer32 := Head_Of(sols).n;
    p_eval : Eval_Poly_Sys(1..nq) := Create(p);
    jac : Jaco_Mat(1..nq,1..nv) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..nq,1..nv) := Create(jac);
    numb : natural32 := 0;
    tol_rnk : constant double_float := 1.0E-6;
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(Length_Of(sols))) := Create(sols);
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,nv);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,nv);
    pl : Point_List;

  begin
    if deflate
     then put_line("Deflation not yet provided for overdetermined systems.");
    end if;
    for i in sa'range loop
      numb := 0;
      sa(i).res := Sum_Norm(Eval(p_eval,sa(i).v));
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        Silent_Gauss_Newton
          (p_eval,jac_eval,sa(i).all,tol_rnk,epsxa,epsfa,numb,max,fail);
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,false,tolsing,epsxa);
      numit := numit + numb;
    end loop;
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
    Clear(jac); Clear(p_eval); Clear(jac_eval);
    Clear(pl);
  end Silent_Root_Sharpener;

  procedure Silent_Root_Sharpener
               ( p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in boolean ) is

    use Standard_Complex_Jaco_Matrices;

    nq : constant integer32 := p'last;
    nv : constant integer32 := Head_Of(sols).n;
    p_eval : Eval_Poly_Sys(1..nq) := Create(p);
    jac : Jaco_Mat(1..nq,1..nv) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..nq,1..nv) := Create(jac);
    numb : natural32 := 0;
    tol_rnk : constant double_float := 1.0E-6;
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(Length_Of(sols))) := Create(sols);
    refsols_last : Solution_List;
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,nv);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,nv);
    pl : Point_List;

  begin
    if deflate
     then put_line("Deflation not yet provided for overdetermined systems.");
    end if;
    for i in sa'range loop
      numb := 0;
      sa(i).res := Sum_Norm(Eval(p_eval,sa(i).v));
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        Silent_Gauss_Newton
          (p_eval,jac_eval,sa(i).all,tol_rnk,epsxa,epsfa,numb,max,fail);
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,false,tolsing,epsxa);
      if not fail
       then Append(refsols,refsols_last,sa(i).all);
      end if;
      numit := numit + numb;
    end loop;
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
    Clear(jac); Clear(p_eval); Clear(jac_eval);
    Clear(pl);
  end Silent_Root_Sharpener;

  procedure Silent_Root_Sharpener
               ( p : in Laur_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32 ) is

    use Standard_Complex_Laur_Jacomats;

    nq : constant integer32 := p'last;
    nv : constant integer32 := Head_Of(sols).n;
    p_eval : Eval_Laur_Sys(1..nq) := Create(p);
    jac : Jaco_Mat(1..nq,1..nv) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..nq,1..nv) := Create(jac);
    numb : natural32 := 0;
    tol_rnk : constant double_float := 1.0E-6;
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(Length_Of(sols))) := Create(sols);
    refsols_last : Solution_List;
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,nv);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,nv);
    pl : Point_List;

  begin
    for i in sa'range loop
      numb := 0;
      sa(i).res := Sum_Norm(Eval(p_eval,sa(i).v));
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
        Silent_Gauss_Newton
          (p_eval,jac_eval,sa(i).all,tol_rnk,epsxa,epsfa,numb,max,fail);
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,false,tolsing,epsxa);
      if not fail
       then Append(refsols,refsols_last,sa(i).all);
      end if;
      numit := numit + numb;
    end loop;
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
    Clear(jac); Clear(p_eval); Clear(jac_eval);
    Clear(pl);
  end Silent_Root_Sharpener;

  procedure Reporting_Root_Sharpener
               ( file : in file_type;
                 p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate,wout : in boolean ) is

    use Standard_Complex_Jaco_Matrices;

    nq : constant integer32 := p'length;
    nv : constant integer32 := Head_Of(sols).n;
    p_eval : Eval_Poly_Sys(1..nq) := Create(p);
    eyp : Standard_Complex_Vectors.Vector(p'range);
    jac : Jaco_Mat(1..nq,1..nv) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..nq,1..nv) := Create(jac);
    numb,nbdef : natural32 := 0;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    nbtot : constant natural32 := Length_Of(sols);
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(nbtot)) := Create(sols);
    initres : Standard_Floating_Vectors.Vector(sa'range);
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15)
                      := Standard_Condition_Tables.Create(15); 
    tolrnk : constant double_float := tolsing*100.0;
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,nv);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,nv);
    pl : Point_List;

  begin
    if deflate
     then put_line("Deflation not yet provided for overdetermined systems.");
    end if;
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); new_line(file);
    put(file,nbtot,1); put(file," ");
    put(file,nv,1); new_line(file);
    put_bar(file);
    for i in sa'range loop 
      numb := 0; nbdef := 0;
      eyp := Eval(p_eval,sa(i).v);
      sa(i).res := Sum_Norm(eyp);
      initres(i) := sa(i).res;
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
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
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,false,tolsing,epsxa);
      Write_Info(file,sa(i).all,initres(i),natural32(i),numb,nbdef,fail,infty);
      Write_Type
        (file,h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,infty,false,
         tolsing,epsxa,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
      Standard_Condition_Tables.Update_Corrector(t_err,sa(i).all);
      Standard_Condition_Tables.Update_Condition(t_rco,sa(i).all);
      Standard_Condition_Tables.Update_Residuals(t_res,sa(i).all);
      numit := numit + numb;
    end loop;
    Write_Global_Info
      (file,nbtot,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
    Clear(jac); Clear(p_eval); Clear(jac_eval);
    Standard_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    Clear(pl);
  end Reporting_Root_Sharpener;

  procedure Reporting_Root_Sharpener
               ( file : in file_type;
                 p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate,wout : in boolean ) is

    use Standard_Complex_Jaco_Matrices;

    nq : constant integer32 := p'length;
    nv : constant integer32 := Head_Of(sols).n;
    p_eval : Eval_Poly_Sys(1..nq) := Create(p);
    eyp : Standard_Complex_Vectors.Vector(p'range);
    jac : Jaco_Mat(1..nq,1..nv) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..nq,1..nv) := Create(jac);
    numb,nbdef : natural32 := 0;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    nbtot : constant natural32 := Length_Of(sols);
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(nbtot)) := Create(sols);
    initres : Standard_Floating_Vectors.Vector(sa'range);
    refsols_last : Solution_List;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15)
                      := Standard_Condition_Tables.Create(15); 
    tolrnk : constant double_float := tolsing*100.0;
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,nv);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,nv);
    pl : Point_List;

  begin
    if deflate
     then put_line("Deflation not yet provided for overdetermined systems.");
    end if;
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); new_line(file);
    put(file,nbtot,1); put(file," "); put(file,nv,1); new_line(file);
    put_bar(file);
    for i in sa'range loop 
      numb := 0; nbdef := 0;
      eyp := Eval(p_eval,sa(i).v);
      sa(i).res := Sum_Norm(eyp);
      initres(i) := sa(i).res;
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
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
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,false,tolsing,epsxa);
      Write_Info(file,sa(i).all,initres(i),natural32(i),numb,nbdef,fail,infty);
      Write_Type
        (file,h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,infty,false,
         tolsing,epsxa,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
      if not fail
       then Append(refsols,refsols_last,sa(i).all);
      end if;
      Standard_Condition_Tables.Update_Corrector(t_err,sa(i).all);
      Standard_Condition_Tables.Update_Condition(t_rco,sa(i).all);
      Standard_Condition_Tables.Update_Residuals(t_res,sa(i).all);
      numit := numit + numb;
    end loop;
    Write_Global_Info
      (file,nbtot,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
    Clear(jac); Clear(p_eval); Clear(jac_eval);
    Standard_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    Clear(pl);
  end Reporting_Root_Sharpener;

  procedure Reporting_Root_Sharpener
               ( file : in file_type;
                 p : in Laur_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean ) is

    use Standard_Complex_Laur_Jacomats;

    nq : constant integer32 := p'length;
    nv : constant integer32 := Head_Of(sols).n;
    p_eval : Eval_Laur_Sys(1..nq) := Create(p);
    eyp : Standard_Complex_Vectors.Vector(p'range);
    jac : Jaco_Mat(1..nq,1..nv) := Create(p);
    jac_eval : Eval_Jaco_Mat(1..nq,1..nv) := Create(jac);
    numb,nbdef : natural32 := 0;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    nbtot : constant natural32 := Length_Of(sols);
    fail,infty : boolean;
    sa : Solution_Array(1..integer32(nbtot)) := Create(sols);
    initres : Standard_Floating_Vectors.Vector(sa'range);
    refsols_last : Solution_List;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..15)
                      := Standard_Condition_Tables.Create(15); 
    tolrnk : constant double_float := tolsing*100.0;
    h1 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,nv);
    h2 : constant Standard_Complex_Vectors.Vector
       := Standard_Random_Vectors.Random_Vector(1,nv); 
    pl : Point_List;

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); new_line(file);
    put(file,nbtot,1); put(file," "); put(file,nv,1); new_line(file);
    put_bar(file);
    for i in sa'range loop 
      numb := 0; nbdef := 0;
      eyp := Eval(p_eval,sa(i).v);
      sa(i).res := Sum_Norm(eyp);
      initres(i) := sa(i).res;
      infty := At_Infinity(sa(i).all,false,1.0E+8);
      if not infty and sa(i).res < 0.1 and sa(i).err < 0.1 then
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
      Multiplicity(h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,
                   infty,false,tolsing,epsxa);
      Write_Info(file,sa(i).all,initres(i),natural32(i),numb,nbdef,fail,infty);
      Write_Type
       (file,h1,h2,pl,sa(i),natural32(i),sa(sa'first..i),fail,infty,false,
        tolsing,epsxa,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
      if not fail
       then Append(refsols,refsols_last,sa(i).all);
      end if;
      Standard_Condition_Tables.Update_Corrector(t_err,sa(i).all);
      Standard_Condition_Tables.Update_Condition(t_rco,sa(i).all);
      Standard_Condition_Tables.Update_Residuals(t_res,sa(i).all);
      numit := numit + numb;
    end loop;
    Write_Global_Info
      (file,nbtot,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
    Clear(jac); Clear(p_eval); Clear(jac_eval);
    Standard_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    Clear(pl);
  end Reporting_Root_Sharpener;

end Standard_Root_Refiners;
